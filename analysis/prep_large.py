#!/usr/bin/env python3
"""
High-throughput ligand preparation for AutoDock / Vina-family docking.

Design goals:
- Stream SDF records once (no full-file rereads per ligand)
- RDKit -> Meeko -> direct PDBQT output
- No per-ligand intermediate SDF files
- Shard outputs into docking-friendly directories
- Write a global manifest and an error log
- Resume safely if output PDBQT already exists

Tested conceptually for large libraries; tune --nproc, --submit-batch, and --shard-size
for your storage and CPU environment.

Requires:
  rdkit
  meeko
Optional for receptor prep:
  pdb2pqr30
  AutoDockTools
"""

import argparse
import csv
import hashlib
import os
import re
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from itertools import islice
from pathlib import Path
import multiprocessing as mp

from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy


# -----------------------------
# General helpers
# -----------------------------
SAFE_CHARS = re.compile(r"[^A-Za-z0-9._-]+")


def safe_filename(text: str, max_len: int = 80) -> str:
    text = (text or "").strip()
    text = SAFE_CHARS.sub("_", text)
    text = text.strip("._")
    if not text:
        text = "lig"
    return text[:max_len]


def batched(iterable, n):
    """Yield lists of length <= n."""
    it = iter(iterable)
    while True:
        batch = list(islice(it, n))
        if not batch:
            break
        yield batch


def iter_sdf_blocks(sdf_path: Path):
    """
    Stream one mol block at a time from an SDF file.
    Avoids loading the whole file and avoids reparsing the full file per ligand.
    """
    buf = []
    with open(sdf_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            buf.append(line)
            if line.startswith("$$$$"):
                yield "".join(buf)
                buf = []
        # handle trailing block if file is missing final "$$$$"
        if buf and "".join(buf).strip():
            yield "".join(buf)


# -----------------------------
# Ligand ID logic
# -----------------------------
def ligand_id_from_mol(mol, mol_block: str, fallback_prefix: str = "lig") -> str:
    """
    Stable, relatively cheap ligand ID policy:
      1) existing name/title if present
      2) hash of canonical SMILES
      3) hash of mol block
    """
    if mol.HasProp("_Name"):
        name = mol.GetProp("_Name").strip()
        if name:
            return safe_filename(name)

    try:
        smiles = Chem.MolToSmiles(mol, canonical=True)
        if smiles:
            return f"{fallback_prefix}_{hashlib.sha1(smiles.encode()).hexdigest()[:16]}"
    except Exception:
        pass

    return f"{fallback_prefix}_{hashlib.sha1(mol_block.encode()).hexdigest()[:16]}"


# -----------------------------
# 3D handling
# -----------------------------
def has_3d(mol) -> bool:
    if mol.GetNumConformers() == 0:
        return False
    conf = mol.GetConformer()
    return conf.Is3D()


def ensure_3d(mol):
    """
    Preserve existing 3D when present.
    If absent, add Hs, embed, and do a short UFF optimization.
    """
    mol = Chem.AddHs(mol, addCoords=True)
    if has_3d(mol):
        return mol

    params = AllChem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        raise RuntimeError("RDKit ETKDG embedding failed")

    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        # UFF can fail for some chemistries; keep embedded geometry if so
        pass

    return mol


# -----------------------------
# Worker
# -----------------------------
def prepare_one(task):
    """
    Worker payload:
      task = {
        "source_file": str,
        "source_idx": int,
        "global_idx": int,
        "mol_block": str,
        "out_pdbqt": str,
        "shard_name": str,
      }
    """
    source_file = task["source_file"]
    source_idx = task["source_idx"]
    global_idx = task["global_idx"]
    mol_block = task["mol_block"]
    out_pdbqt = Path(task["out_pdbqt"])
    shard_name = task["shard_name"]

    if out_pdbqt.exists():
        ligand_file = out_pdbqt.name
        ligand_id = out_pdbqt.stem
        return {
            "status": "skipped",
            "source_file": source_file,
            "source_idx": source_idx,
            "global_idx": global_idx,
            "ligand_id": ligand_id,
            "ligand_file": ligand_file,
            "shard": shard_name,
            "out_pdbqt": str(out_pdbqt),
            "error": "",
        }

    try:
        mol = Chem.MolFromMolBlock(mol_block, sanitize=True, removeHs=False)
        if mol is None:
            raise ValueError("RDKit failed to parse mol block")

        ligand_id = ligand_id_from_mol(mol, mol_block)
        mol.SetProp("_Name", ligand_id)

        mol = ensure_3d(mol)

        preparator = MoleculePreparation()
        setups = preparator.prepare(mol)
        if not setups:
            raise ValueError("Meeko returned no setups")

        # Most workflows write the first setup; if you want all tautomers/protonation
        # variants written, that can be expanded later.
        pdbqt_string, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            raise ValueError(err if err else "Meeko PDBQTWriterLegacy failed")

        out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
        # Add a name remark like your current script does
        with open(out_pdbqt, "w") as fh:
            fh.write(f"REMARK Name = {ligand_id}\n")
            fh.write(pdbqt_string)

        return {
            "status": "ok",
            "source_file": source_file,
            "source_idx": source_idx,
            "global_idx": global_idx,
            "ligand_id": ligand_id,
            "ligand_file": out_pdbqt.name,
            "shard": shard_name,
            "out_pdbqt": str(out_pdbqt),
            "error": "",
        }

    except Exception as e:
        return {
            "status": "error",
            "source_file": source_file,
            "source_idx": source_idx,
            "global_idx": global_idx,
            "ligand_id": "",
            "ligand_file": "",
            "shard": shard_name,
            "out_pdbqt": str(out_pdbqt),
            "error": str(e),
        }


# -----------------------------
# Optional receptor prep
# -----------------------------
def strip_waters_from_pdb(in_pdb: Path, out_pdb: Path):
    with open(in_pdb) as f:
        lines = [l for l in f if l.startswith("ATOM") or l.startswith("HETATM")]
    dry_lines = [l for l in lines if "HOH" not in l]
    with open(out_pdb, "w") as f:
        f.write("".join(dry_lines))


def prepare_receptor_legacy(pdb_path: Path, out_pdbqt: Path, del_water: bool = False):
    """
    Optional receptor preparation using your current-style legacy approach:
      pdb -> (optional dry pdb) -> pqr -> pdbqt
    """
    try:
        import AutoDockTools
    except ImportError as e:
        raise RuntimeError(
            "AutoDockTools is required for --receptor prep in this legacy mode"
        ) from e

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)

        work_pdb = pdb_path
        if del_water:
            dry_pdb = td / f"{pdb_path.stem}_dry.pdb"
            strip_waters_from_pdb(pdb_path, dry_pdb)
            work_pdb = dry_pdb

        pqr_path = td / f"{pdb_path.stem}.pqr"
        subprocess.run(
            ["pdb2pqr30", "--ff=AMBER", str(work_pdb), str(pqr_path)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        prepare_receptor = (
            Path(AutoDockTools.__path__[0]) / "Utilities24" / "prepare_receptor4.py"
        )
        subprocess.run(
            ["python3", str(prepare_receptor), "-r", str(pqr_path), "-o", str(out_pdbqt)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


# -----------------------------
# Main
# -----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Scalable ligand preparation for large virtual screening libraries."
    )
    parser.add_argument(
        "input_path",
        type=str,
        help="Path to either a sinlge .sdf file or a directory containing one or more .sdf files.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        required=True,
        help="Output directory for prepared ligands and manifests.",
    )
    parser.add_argument(
        "--nproc",
        type=int,
        default=max(1, os.cpu_count() // 2),
        help="Number of worker processes for ligand prep.",
    )
    parser.add_argument(
        "--shard-size",
        type=int,
        default=5000,
        help="Number of ligands per shard directory.",
    )
    parser.add_argument(
        "--submit-batch",
        type=int,
        default=256,
        help="Number of ligands submitted to the process pool at a time.",
    )
    parser.add_argument(
        "--receptor",
        type=str,
        default=None,
        help="Optional receptor PDB file to prepare once.",
    )
    parser.add_argument(
        "--del-water",
        action="store_true",
        help="If receptor prep is requested, strip HOH waters first.",
    )

    args = parser.parse_args()


    input_path = Path(args.input_path)

    input_path = Path(args.input_path)

    if input_path.is_file():
        if input_path.suffix.lower() != ".sdf":
            parser.error(f"If input_path is a file, it must be an .sdf file: {input_path}")
        sdf_files = [input_path]

    elif input_path.is_dir():
        sdf_files = sorted(input_path.glob("*.sdf"))
        if not sdf_files:
            parser.error(f"No .sdf files found in directory: {input_path}")

    else:
        parser.error(f"Input path does not exist: {input_path}")


    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    shards_root = outdir / "shards"
    shards_root.mkdir(parents=True, exist_ok=True)

    # Optional receptor prep
    if args.receptor:
        receptor_path = Path(args.receptor)
        receptor_out = outdir / f"{receptor_path.stem}_rec.pdbqt"
        if receptor_out.exists():
            print(f"[receptor] exists: {receptor_out}")
        else:
            print(f"[receptor] preparing: {receptor_path.name} -> {receptor_out.name}")
            prepare_receptor_legacy(receptor_path, receptor_out, del_water=args.del_water)
            print(f"[receptor] done: {receptor_out}")

    if not sdf_files:
        raise SystemExit("No .sdf files found in input_dir")

    manifest_path = outdir / "ligand_manifest.tsv"
    errors_path = outdir / "ligand_errors.tsv"

    manifest_fh = open(manifest_path, "w", newline="")
    errors_fh = open(errors_path, "w", newline="")

    manifest_writer = csv.writer(manifest_fh, delimiter="\t")
    error_writer = csv.writer(errors_fh, delimiter="\t")

    manifest_writer.writerow(
        [
            "global_idx",
            "source_file",
            "source_idx",
            "ligand_id",
            "ligand_file",
            "shard",
            "out_pdbqt",
            "status",
        ]
    )
    error_writer.writerow(
        [
            "global_idx",
            "source_file",
            "source_idx",
            "shard",
            "out_pdbqt",
            "error",
        ]
    )

    total_seen = 0
    total_ok = 0
    total_skipped = 0
    total_error = 0

    # spawn is safer with many cheminformatics libs than fork on some systems
    mp_ctx = mp.get_context("spawn")

    with ProcessPoolExecutor(max_workers=args.nproc, mp_context=mp_ctx) as ex:
        for sdf_file in sdf_files:
            print(f"[input] streaming {sdf_file.name}")

            # enumerate molecules in this file once, in streaming mode
            stream = enumerate(iter_sdf_blocks(sdf_file))

            for batch in batched(stream, args.submit_batch):
                tasks = []
                for source_idx, mol_block in batch:
                    global_idx = total_seen
                    shard_idx = global_idx // args.shard_size
                    shard_name = f"shard_{shard_idx:05d}"
                    shard_dir = shards_root / shard_name / "ligands"
                    ligand_stub = f"lig_{global_idx:09d}"
                    out_pdbqt = shard_dir / f"{ligand_stub}.pdbqt"

                    tasks.append(
                        {
                            "source_file": sdf_file.name,
                            "source_idx": source_idx,
                            "global_idx": global_idx,
                            "mol_block": mol_block,
                            "out_pdbqt": str(out_pdbqt),
                            "shard_name": shard_name,
                        }
                    )
                    total_seen += 1

                for result in ex.map(prepare_one, tasks, chunksize=16):
                    if result["status"] in ("ok", "skipped"):
                        manifest_writer.writerow(
                            [
                                result["global_idx"],
                                result["source_file"],
                                result["source_idx"],
                                result["ligand_id"],
                                result["ligand_file"],
                                result["shard"],
                                result["out_pdbqt"],
                                result["status"],
                            ]
                        )
                        if result["status"] == "ok":
                            total_ok += 1
                        else:
                            total_skipped += 1
                    else:
                        total_error += 1
                        error_writer.writerow(
                            [
                                result["global_idx"],
                                result["source_file"],
                                result["source_idx"],
                                result["shard"],
                                result["out_pdbqt"],
                                result["error"],
                            ]
                        )

                manifest_fh.flush()
                errors_fh.flush()
                print(
                    f"[progress] seen={total_seen:,} ok={total_ok:,} "
                    f"skipped={total_skipped:,} error={total_error:,}"
                )

    manifest_fh.close()
    errors_fh.close()

    print("\nDone.")
    print(f"Prepared/kept ligands: {total_ok + total_skipped:,}")
    print(f"  Newly prepared:      {total_ok:,}")
    print(f"  Skipped (existing):  {total_skipped:,}")
    print(f"  Errors:              {total_error:,}")
    print(f"Manifest:              {manifest_path}")
    print(f"Errors:                {errors_path}")
    print(f"Shards root:           {shards_root}")


if __name__ == "__main__":
    main()

