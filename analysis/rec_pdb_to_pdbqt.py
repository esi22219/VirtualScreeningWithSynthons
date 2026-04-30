import contextlib
import os
import subprocess
import tempfile
from pathlib import Path
import argparse
import AutoDockTools

def supress_stdout(func):
    def wrapper(*a, **ka):
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                return func(*a, **ka)
    return wrapper

class PrepProt(object):
    def __init__(self, pdb_file):
        self.prot = pdb_file

    def del_water(self, dry_pdb_file):
        with open(self.prot) as f:
            lines = [
                l for l in f.readlines()
                if l.startswith('ATOM') or l.startswith('HETATM')
            ]
        dry_lines = [l for l in lines if 'HOH' not in l]
        with open(dry_pdb_file, 'w') as f:
            f.write(''.join(dry_lines))
        self.prot = dry_pdb_file

    def addH(self, prot_pqr):
        self.prot_pqr = prot_pqr
        subprocess.Popen(
            ['pdb2pqr30', '--ff=AMBER', self.prot, self.prot_pqr],
            stderr=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL
        ).communicate()

    @supress_stdout
    def get_pdbqt(self, prot_pdbqt):
        prepare_receptor = os.path.join(
            AutoDockTools.__path__[0],
            'Utilities24/prepare_receptor4.py'
        )
        subprocess.Popen(
            ['python3', prepare_receptor, '-r', self.prot_pqr, '-o', prot_pdbqt],
            stderr=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL
        ).communicate()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_dir", type=str)
    parser.add_argument("--outdir", type=str, default=None)
    args = parser.parse_args()

    pdb_dir = Path(args.pdb_dir)
    outdir = Path(args.outdir) if args.outdir else pdb_dir
    outdir.mkdir(parents=True, exist_ok=True)

    for pdb_file in pdb_dir.iterdir():
        if pdb_file.suffix.lower() != ".pdb":
            continue

        if "protein" in pdb_file.name:
            out_pdbqt = outdir / f"{pdb_file.stem}.pdbqt"
        else:
            out_pdbqt = outdir / f"{pdb_file.stem}_protein.pdbqt"

        if out_pdbqt.exists():
            print(f"Skipping {pdb_file.name} (exists)")
            continue

        print(f"Processing protein {pdb_file.name}")

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            prep = PrepProt(pdb_file)

            dry_pdb = tmp / f"{pdb_file.stem}_dry.pdb"
            prep.del_water(dry_pdb)

            pqr_file = tmp / f"{pdb_file.stem}.pqr"
            prep.addH(pqr_file)

            prep.get_pdbqt(out_pdbqt)

        print(f"Saved pdbqt to {out_pdbqt}")
