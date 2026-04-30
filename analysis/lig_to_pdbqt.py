import contextlib
import os
import tempfile
from pathlib import Path
import argparse

import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, obutils
from openbabel import pybel

def supress_stdout(func):
    def wrapper(*a, **ka):
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                return func(*a, **ka)
    return wrapper

class PrepLig(object):
    def __init__(self, input_mol, mol_format):
        if mol_format == 'smi':
            self.ob_mol = pybel.readstring('smi', input_mol)
        elif mol_format == 'sdf':
            self.ob_mol = next(pybel.readfile(mol_format, input_mol))
        else:
            raise ValueError(f'mol_format {mol_format} not supported')

    def addH(self):
        self.ob_mol.OBMol.AddHydrogens(True, True, 7)

    def gen_conf(self):
        sdf_block = self.ob_mol.write('sdf')
        rdkit_mol = Chem.MolFromMolBlock(sdf_block, removeHs=False)
        AllChem.EmbedMolecule(rdkit_mol, Chem.rdDistGeom.ETKDGv3())
        self.ob_mol = pybel.readstring(
            'sdf',
            Chem.MolToMolBlock(rdkit_mol)
        )

    @supress_stdout
    def get_pdbqt(self, lig_pdbqt):
        preparator = MoleculePreparation()
        preparator.prepare(self.ob_mol.OBMol)
        preparator.write_pdbqt_file(lig_pdbqt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("lig_dir", type=str)
    parser.add_argument("--outdir", type=str, default=None)
    args = parser.parse_args()

    lig_dir = Path(args.lig_dir)
    outdir = Path(args.outdir) if args.outdir else lig_dir
    outdir.mkdir(parents=True, exist_ok=True)

    for lig_file in lig_dir.iterdir():
        mol_format = lig_file.suffix.lower().lstrip('.')
        if mol_format not in ('sdf', 'smi'):
            continue

        if "ligand" in lig_file.name:
            out_pdbqt = outdir / f"{lig_file.stem}.pdbqt"
        else:
            out_pdbqt = outdir / f"{lig_file.stem}_ligand.pdbqt"

        if out_pdbqt.exists():
            print(f"Skipping {lig_file.name} (exists)")
            continue

        print(f"Processing ligand {lig_file.name}")

        with tempfile.TemporaryDirectory():
            lig_prep = PrepLig(str(lig_file), mol_format)
            lig_prep.addH()
            lig_prep.gen_conf()
            lig_prep.get_pdbqt(out_pdbqt)

        print(f"Saved pdbqt to {out_pdbqt}")
