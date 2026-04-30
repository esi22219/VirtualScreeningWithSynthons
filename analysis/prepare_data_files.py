import contextlib
import os
import subprocess
import tempfile
from pathlib import Path
import argparse
import multiprocessing as mp
from tqdm import tqdm
import haslib

import AutoDockTools
import rdkit.Chem as Chem
from rdkit.Chem import AllChem, inchi
from meeko import MoleculePreparation
from openbabel import pybel


def supress_stdout(func):
    def wrapper(*a, **ka):
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                return func(*a, **ka)
    return wrapper


# -----------------------------
# Ligand preparation
# -----------------------------


def generate_ligand_id(rdkit_mol):
    """
    Generate a deterministic structure-based ID.
    Uses InChIKey if possible, otherwise SMILES hash.
    """
    try:
        ik = inchi.MolToInchiKey(rdkit_mol)
        if ik:
            return f"LIG_{ik.replace('-', '_')}"
    except Exception:
        pass

    # Fallback: canonical SMILES hash
    smiles = Chem.MolToSmiles(rdkit_mol, canonical=True)
    h = hashlib.sha1(smiles.encode()).hexdigest()[:10].upper()
    return f"LIG_SMILES_{h}"


class PrepLig(object):
    def __init__(self, input_mol, mol_format):
        self.input_path = input_mol
        self.mol_format = mol_format

        if mol_format == "smi":
            self.ob_mol = pybel.readstring("smi", input_mol)
            rdkit_mol = Chem.MolFromSmiles(input_mol)
            self.ligand_id = generate_ligand_id(rdkit_mol)

        elif mol_format == "sdf":
            self.ob_mol = next(pybel.readfile("sdf", input_mol))

            # RDKit copy for identity generation
            rdkit_mol = Chem.MolFromMolFile(
                input_mol,
                sanitize=True,
                removeHs=False
            )

            # 1️⃣ Prefer existing SDF title
            title = (self.ob_mol.title or "").strip()
            if title:
                self.ligand_id = title
            else:
                self.ligand_id = generate_ligand_id(rdkit_mol)

                # Insert generated ID back into SDF title
                self.ob_mol.title = self.ligand_id

        else:
            raise ValueError(f"Unsupported ligand format {mol_format}")

    def addH(self):
        self.ob_mol.OBMol.AddHydrogens(True, True, 7)

    def gen_conf(self):
        if not self.ob_mol.OBMol.Has3D():
            self.ob_mol.make3D()

    
    def write_identified_sdf(self, out_sdf):
        self.ob_mol.write("sdf", str(out_sdf), overwrite=True)

    @supress_stdout
    def get_pdbqt(self, out_pdbqt):
        preparator = MoleculePreparation()
        preparator.prepare(self.ob_mol.OBMol)
        preparator.write_pdbqt_file(out_pdbqt)

        # Inject LigandID into PDBQT as REMARK
        with open(out_pdbqt, "r") as f:
            lines = f.readlines()

        with open(out_pdbqt, "w") as f:
            f.write(f"REMARK LigandID: {self.ligand_id}\n")
            f.writelines(lines)

# this is not handling many ligands super well, good for small amounts but def not 500k
def _prepare_single_ligand(args):
    """
    Worker for multiprocessing ligand preparation.
    Re-loads molecule inside the process (pickle-safe).
    """
    source, mol_index, out_pdbqt = args

    if out_pdbqt.exists():
        return f"Skipping ligand {out_pdbqt.name} (exists)"

    # Reconstruct the molecule INSIDE the worker
    if source.suffix.lower() == ".sdf":
        mols = list(pybel.readfile("sdf", str(source)))
        ob_mol = mols[mol_index]
    elif source.suffix.lower() == ".smi":
        smi = source.read_text().strip()
        ob_mol = pybel.readstring("smi", smi)
    else:
        return f"Unsupported ligand source {source}"

    with tempfile.TemporaryDirectory():
        lig = PrepLig(ob_mol)
        lig.addH()
        lig.gen_conf()
        lig.get_pdbqt(out_pdbqt)

    return f"Saved ligand pdbqt -> {out_pdbqt}"


# -----------------------------
# Protein preparation
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
    parser.add_argument("input_dir", type=str, help="Directory with PDB / SDF / SMI files")
    parser.add_argument("--outdir", type=str, default=None)
    parser.add_argument("--nproc", type=int, default=1)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    outdir = Path(args.outdir) if args.outdir else input_dir
    outdir.mkdir(parents=True, exist_ok=True)
    lig_sdf_dir = outdir / "ligands_sdf"
    lig_pdbqt_dir = outdir / "ligands_pdbqt"

    lig_sdf_dir.mkdir(exist_ok=True)
    lig_pdbqt_dir.mkdir(exist_ok=True)

    ligand_tasks = []
    
    for file in input_dir.iterdir():
        ext = file.suffix.lower()

        # -------------------------
        # Protein (PDB)
        # -------------------------
        if ext == ".pdb":
            if "protein" in file.name:
                out_pdbqt = outdir / f"{file.stem}.pdbqt"
            else:
                out_pdbqt = outdir / f"{file.stem}_protein.pdbqt"

            if out_pdbqt.exists():
                print(f"Skipping protein {file.name} (exists)")
                continue

            print(f"Processing protein {file.name}")

            with tempfile.TemporaryDirectory() as tmpdir:
                tmp = Path(tmpdir)
                prep = PrepProt(file)

                dry_pdb = tmp / f"{file.stem}_dry.pdb"
                prep.del_water(dry_pdb)

                pqr = tmp / f"{file.stem}.pqr"
                prep.addH(pqr)

                prep.get_pdbqt(out_pdbqt)

            print(f"Saved protein pdbqt -> {out_pdbqt}")

        # -------------------------
        # Ligand (SDF / SMI)
        # -------------------------
                   
        elif ext == ".sdf":
            print(f"Processing multi-molecule SDF {file.name}")
            mol_count = sum(1 for _ in pybel.readfile("sdf", str(file)))

            if mol_count == 0:
                print(f"No molecules found in {file.name}")
                continue
            
        
            for i in range(mol_count):
                out_pdbqt = outdir / f"{file.stem}{i}_ligand.pdbqt"
                ligand_tasks.append((file, i, out_pdbqt))

            
        elif ext == ".smi":
            print(f"Processing SMILES {file.name}")

            ob_mol = pybel.readstring("smi", file.read_text().strip())
        
            out_pdbqt = outdir / f"{file.stem}_ligand.pdbqt"

            if out_pdbqt.exists():
                print(f"Skipping ligand {file.name} (exists)")
                continue

            with tempfile.TemporaryDirectory():
                lig = PrepLig(ob_mol)
                lig.addH()
                lig.gen_conf()
                lig.get_pdbqt(out_pdbqt)

            print(f"Saved ligand pdbqt -> {out_pdbqt}")

        # -------------------------
        # Ignore everything else
        # -------------------------
        else:
            continue

    if ligand_tasks:
        n_procs = args.nproc
        print(f"\nPreparing {len(ligand_tasks)} ligands using {n_procs} processes")

        with mp.Pool(processes=n_procs) as pool:
            for msg in tqdm(
                pool.imap_unordered(_prepare_single_ligand, ligand_tasks),
                total=len(ligand_tasks),
                desc="Ligand preparation",
            ):
             #   if msg:
                   # print(msg)

