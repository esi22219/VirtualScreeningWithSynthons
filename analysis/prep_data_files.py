import contextlib
import os
import subprocess
import tempfile
from pathlib import Path
import argparse
import multiprocessing as mp
import hashlib
from tqdm import tqdm

import AutoDockTools
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Chem.inchi import MolToInchiKey
from meeko import MoleculePreparation
from openbabel import pybel


def supress_stdout(func):
    def wrapper(*a, **ka):
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                return func(*a, **ka)
    return wrapper


# -----------------------------
# Ligand ID helpers
# -----------------------------
def get_ligand_id(ob_mol):
    """
    Determine a stable ligand ID.
    Priority:
      1) Existing title
      2) InChIKey
      3) Hash of canonical SMILES
    """
    title = ob_mol.title.strip()
    if title:
        return title.replace(" ", "_")

    # Convert to RDKit mol for InChIKey
    sdf_block = ob_mol.write("sdf")
    rd_mol = Chem.MolFromMolBlock(sdf_block, sanitize=True)

    if rd_mol:
        try:
            return MolToInchiKey(rd_mol)
        except Exception:
            pass

        smiles = Chem.MolToSmiles(rd_mol, canonical=True)
        return hashlib.sha1(smiles.encode()).hexdigest()[:16]

    return hashlib.sha1(sdf_block.encode()).hexdigest()[:16]


def write_sdf_with_id(ob_mol, out_sdf, ligand_id):
    ob_mol.title = ligand_id
    with open(out_sdf, "w") as f:
        f.write(ob_mol.write("sdf"))


# -----------------------------
# Ligand preparation
# -----------------------------
class PrepLig(object):
    def __init__(self, ob_mol):
        self.ob_mol = ob_mol

    def addH(self):
        self.ob_mol.OBMol.AddHydrogens(True, True, 7)

    def gen_conf(self):
        if not self.ob_mol.OBMol.Has3D():
            self.ob_mol.make3D()

    @supress_stdout
    def get_pdbqt(self, out_pdbqt):
        preparator = MoleculePreparation()
        preparator.prepare(self.ob_mol.OBMol)
        preparator.write_pdbqt_file(out_pdbqt)
    
    def insert_id_title(self, pdbqt_file):

        with open(pdbqt_file, "r") as f:
            lines = f.readlines()

        # Avoid duplicating remark if rerun
        for line in lines:
            if line.startswith("REMARK") and "Name =" in line:
                return

        new_lines = [f"REMARK Name = {self.ob_mol.title}\n"] + lines

        with open(pdbqt_file, "w") as f:
            f.writelines(new_lines)


def _prepare_single_ligand(args):
    source, mol_index, lig_index, sdf_dir, pdbqt_dir, index_records = args

    # Reload molecule
    mols = list(pybel.readfile("sdf", str(source)))
    ob_mol = mols[mol_index]
    
    ligand_id = get_ligand_id(ob_mol)
    
    if not ob_mol.title:
        ob_mol.title = ligand_id
    
    out_sdf = sdf_dir / f"ligand{lig_index}.sdf"
    out_pdbqt = pdbqt_dir / f"ligand{lig_index}.pdbqt"

    if out_pdbqt.exists():
        return f"Skipping {source.stem}{i}_ligand (exists)"

    write_sdf_with_id(ob_mol, out_sdf, ligand_id)

    with tempfile.TemporaryDirectory():

        lig = PrepLig(ob_mol)
        lig.addH()
        lig.gen_conf()
        lig.get_pdbqt(out_pdbqt)
        lig.insert_id_title(out_pdbqt) 

    index_records.append((f"ligand{lig_index}", ligand_id))
    
    return f"Saved ligand{lig_index}"



# -----------------------------
# Protein preparation (unchanged)
# -----------------------------
class PrepProt(object):
    def __init__(self, pdb_file):
        self.prot = pdb_file

    def del_water(self, dry_pdb_file):
        with open(self.prot) as f:
            lines = [
                l for l in f.readlines()
                if l.startswith("ATOM") or l.startswith("HETATM")
            ]
        dry_lines = [l for l in lines if "HOH" not in l]
        with open(dry_pdb_file, "w") as f:
            f.write("".join(dry_lines))
        self.prot = dry_pdb_file

    def addH(self, prot_pqr):
        self.prot_pqr = prot_pqr
        subprocess.Popen(
            ["pdb2pqr30", "--ff=AMBER", self.prot, self.prot_pqr],
            stderr=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL
        ).communicate()

    @supress_stdout
    def get_pdbqt(self, prot_pdbqt):
        prepare_receptor = os.path.join(
            AutoDockTools.__path__[0],
            "Utilities24/prepare_receptor4.py"
        )
        subprocess.Popen(
            ["python3", prepare_receptor, "-r", self.prot_pqr, "-o", prot_pdbqt],
            stderr=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL
        ).communicate()


# -----------------------------
# Main
# -----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=str)
    parser.add_argument("-o", "--outdir", type=str, default=None, help="Directory to save converted files")
    parser.add_argument("--del_water", action="store_true", help="Remove Water from Protein Structure")
    parser.add_argument("--nproc", type=int, default=1)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    outdir = Path(args.outdir) if args.outdir else input_dir

    sdf_dir = outdir / "sdf_ligands"
    pdbqt_dir = outdir / "pdbqt_ligands"
    sdf_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_dir.mkdir(parents=True, exist_ok=True)

    ligand_tasks = []
    lig_counter = 0
    index_records = mp.Manager().list()
    for file in input_dir.iterdir():
        ext = file.suffix.lower()

        if ext == ".pdb":
            print(f"Converting pdb receptor {file}")
            out_pdbqt = outdir / f"{file.stem}_rec.pdbqt"
            if out_pdbqt.exists():
                print(f"Skipping file (exists)")
                continue

            with tempfile.TemporaryDirectory() as tmpdir:
                tmp = Path(tmpdir)
                prep = PrepProt(file)
                if args.del_water:
                    dry_pdb = tmp / f"{file.stem}_dry.pdb"
                    prep.del_water(dry_pdb)
                pqr = tmp / f"{file.stem}.pqr"
                prep.addH(pqr)
                prep.get_pdbqt(out_pdbqt)

        elif ext == ".sdf":
            mol_count = sum(1 for _ in pybel.readfile("sdf", str(file)))
            for i in range(mol_count):
                ligand_tasks.append((file, i, lig_counter, sdf_dir, pdbqt_dir, index_records))
                lig_counter += 1

    if ligand_tasks:
        with mp.Pool(processes=args.nproc) as pool:
            for msg in tqdm(
                pool.imap_unordered(_prepare_single_ligand, ligand_tasks),
                total=len(ligand_tasks),
                desc="Ligand preparation"
            ):
                if msg:
                    print(msg)
    index_file = outdir / "ligand_index.tsv"
    with open(index_file, "w") as f:
        f.write("ligand_file\tligand_id\n")
        for lig_file, lig_id in index_records:
            f.write(f"{lig_file}\t{lig_id}\n")
