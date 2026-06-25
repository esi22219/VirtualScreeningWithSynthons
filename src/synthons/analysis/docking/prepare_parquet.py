#!/usr/bin/env python3
import argparse
import csv
import hashlib
import os
import shutil
import re
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from itertools import islice
from pathlib import Path
import multiprocessing as mp
import AutoDockTools

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from tqdm import tqdm
import pyarrow.parquet as pq


# General helpers

SAFE_CHARS = re.compile(r"[^A-Za-z0-9._-]+")

def safe_filename(text: str, max_len: int = 120) -> str:
    text = (text or "").strip()
    text = SAFE_CHARS.sub("_", text)
    text = text.strip("._")
    if not text:
        text = "lig"
    return text[:max_len]


def normalize_text(value):
    if value is None:
        return ""
    # pandas / pyarrow missing values can come through as NaN/NA
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return str(value)


def normalize_bool(value) -> bool:
    if value is None:
        return False
    try:
        if pd.isna(value):
            return False
    except Exception:
        pass
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    s = str(value).strip().lower()
    return s in {"1", "true", "t", "yes", "y"}


def batched(iterable, n):
    """Yield lists of length <= n."""
    it = iter(iterable)
    while True:
        batch = list(islice(it, n))
        if not batch:
            break
        yield batch

def clean_marked_product_smiles(product_smiles: str, placeholder_symbols=("W","V","Rf", "Db", "Sg", "Rh5",), placeholder_atomic_nums=(0,)):
    mol = Chem.MolFromSmiles(product_smiles, sanitize=True)
    if mol is None:
        raise ValueError("RDKit failed to parse product_smiles")

    # Track original atom indices across atom deletions so that we can repair
    # placeholder-bearing neighbors after editing.
    for atom in mol.GetAtoms():
        atom.SetIntProp("orig_idx", atom.GetIdx())

    rw = Chem.RWMol(mol)
    touched_orig_idxs = set()
    remove_idxs = []

    for atom in rw.GetAtoms():
        is_placeholder = (
            atom.GetSymbol() in placeholder_symbols
            or atom.GetAtomicNum() in placeholder_atomic_nums
        )
        if is_placeholder:
            remove_idxs.append(atom.GetIdx())
            for nbr in atom.GetNeighbors():
                if nbr.HasProp("orig_idx"):
                    touched_orig_idxs.add(nbr.GetIntProp("orig_idx"))

    for idx in sorted(remove_idxs, reverse=True):
        rw.RemoveAtom(idx)

    mol = rw.GetMol()

    # Remove atom-map annotations from the chemistry we will dock.
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)

    # Repair common aromatic N mark case: [n:20]([W]) -> [nH]
    # Replace each removed placeholder with one hydrogen on the neighboring atom.
    # Examples:
    #   [n]([W])  -> [nH]
    #   [NH][W]   -> NH2
    #   [CH][W]   -> CH2
    for atom in mol.GetAtoms():
        if atom.HasProp("orig_idx") and atom.GetIntProp("orig_idx") in touched_orig_idxs:
            # clear any radical left behind by placeholder removal
            if atom.GetNumRadicalElectrons() > 0:
                atom.SetNumRadicalElectrons(0)

            # add one H to replace the removed placeholder attachment
            atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)

            # keep that H explicit so the replacement is deterministic
            atom.SetNoImplicit(True)


    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)

    for atom in mol.GetAtoms():
        if atom.HasProp("orig_idx"):
            atom.ClearProp("orig_idx")

    clean_smiles = Chem.MolToSmiles(mol, canonical=True)
    dedupe_key = hashlib.sha1(clean_smiles.encode("utf-8")).hexdigest()

    return mol, clean_smiles, len(remove_idxs), dedupe_key

# Parquet streaming

def count_parquet_rows(parquet_files):
    total = 0
    if pq is not None:
        for pf in parquet_files:
            total += pq.ParquetFile(pf).metadata.num_rows
    else:
        for pf in parquet_files:
            total += len(pd.read_parquet(pf))
    return total

def iter_parquet_records(parquet_path: Path, batch_rows: int = 5000):
    """
    Yield (source_idx, record_dict) from a parquet file.

    Uses pyarrow batch streaming when available. Falls back to pandas.read_parquet,
    which may read the full file into memory.
    """
    if pq is not None:
        pf = pq.ParquetFile(parquet_path)
        source_idx = 0
        for rb in pf.iter_batches(batch_size=batch_rows):
            batch = rb.to_pylist()
            for rec in batch:
                yield source_idx, rec
                source_idx += 1
        return

    # Fallback: full read
    df = pd.read_parquet(parquet_path)
    for source_idx, rec in enumerate(df.to_dict(orient="records")):
        yield source_idx, rec



# Ligand ID / provenance logic
def ligand_stub_from_record(source_file: str, run_id: str, route_id: str) -> str:
    source_stem = Path(source_file).stem
    run_id = normalize_text(run_id) or "run"
    route_id = normalize_text(route_id) or "route"
    base = f"{source_stem}__{run_id}__{route_id}"
    digest = hashlib.sha1(base.encode("utf-8")).hexdigest()[:10]
    base_safe = safe_filename(base, max_len=90)
    return f"{base_safe}__{digest}"


# 3D handling
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

# Worker
def prepare_one(task):
    """
    Worker payload:
      task = {
        "source_file": str,
        "source_idx": int,
        "global_idx": int,
        "record": dict,
        "ligand_id": str,
        "out_pdbqt": str,
        "shard_name": str,
      }
    """
    source_file = task["source_file"]
    source_idx = task["source_idx"]
    global_idx = task["global_idx"]
    rec = task["record"]
    ligand_id = task["ligand_id"]
    out_pdbqt = Path(task["out_pdbqt"])
    shard_name = task["shard_name"]
    mol = task["mol"]    
    canonical_prod_smiles = task["canonical_prod_smiles"]

    route_id = normalize_text(rec.get("route_id"))
    run_id = normalize_text(rec.get("run_id"))
    seed_id = normalize_text(rec.get("seed_id"))
    synthon_id = normalize_text(rec.get("synthon_id"))
    reaction_id = normalize_text(rec.get("reaction_id"))
    reaction_name = normalize_text(rec.get("reaction_name"))
    product_smiles = normalize_text(rec.get("product_smiles"))
    synthon_smiles = normalize_text(rec.get("synthon_smiles"))
    failure_reason = normalize_text(rec.get("failure_reason"))
    product_valid = normalize_bool(rec.get("product_valid"))
    

    invalid_atom_types = {"B", "Si"}
   # print(product_smiles)
    base_result = {
        "source_file": source_file,
        "source_idx": source_idx,
        "global_idx": global_idx,
        "route_id": route_id,
        "run_id": run_id,
        "seed_id": seed_id,
        "synthon_id": synthon_id,
        "reaction_id": reaction_id,
        "reaction_name": reaction_name,
        "synthon_smiles": synthon_smiles,
        "product_smiles": product_smiles,
        "ligand_id": ligand_id,
        "ligand_file": out_pdbqt.name,
        "shard": shard_name,
        "out_pdbqt": str(out_pdbqt),
        "duplicate": False,
        "error": "",
    }

    if out_pdbqt.exists():
        
        return {**base_result, "status": "skipped"}

    try:
        if not product_valid:
            raise ValueError(failure_reason or "product_valid is False")
        if not product_smiles:
            raise ValueError("Missing product_smiles")

        #mol, cleaned_smiles, num_removed = clean_marked_product_smiles(product_smiles)
                
        if mol is None:
            raise ValueError("RDKit failed to parse product_smiles")
        
        
        bad_elements = sorted({
            atom.GetSymbol()
            for atom in mol.GetAtoms()
            if atom.GetSymbol() in invalid_atom_types
            })
        if bad_elements:
            raise ValueError(f"Unsupported element(s) for QuickVina2.1 GPU: {bad_elements}")

        mol.SetProp("_Name", ligand_id)
        mol = ensure_3d(mol)
        
        preparator = MoleculePreparation(rigid_macrocycles=True)
        setups = preparator.prepare(mol)
        
        if not setups:
            raise ValueError("Meeko returned no setups")
        # returns a string which can be immediatly parsed or saved using below function from meeko
        pdbqt_string, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            raise ValueError(err if err else "Meeko PDBQTWriterLegacy failed")

        out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
        with open(out_pdbqt, "w") as fh:
            fh.write(f"REMARK Name = {ligand_id}\n")
            fh.write(f"REMARK SourceParquet = {source_file}\n")
            fh.write(f"REMARK RunID = {run_id}\n")
            fh.write(f"REMARK RouteID = {route_id}\n")
            if seed_id:
                fh.write(f"REMARK SeedID = {seed_id}\n")
            if synthon_id:
                fh.write(f"REMARK SynthonID = {synthon_id}\n")
            if synthon_smiles:
                fh.write(f"REMARK SynthonSmiles = {synthon_smiles}\n")
            if reaction_id:
                fh.write(f"REMARK ReactionID = {reaction_id}\n")
                fh.write(f"REMARK ReactionID = {reaction_id}\n")
                fh.write(f"REMARK ReactionID = {reaction_id}\n")
            if reaction_name:
                fh.write(f"REMARK ReactionName = {reaction_name}\n")
            if product_smiles:
                fh.write(f"REMARK ProductSMILES = {product_smiles}\n")
                fh.write(f"REMARK CanonicalProductSMILES = {canonical_prod_smiles}\n")
            fh.write(pdbqt_string)

        return {**base_result, "status": "ok"}

    except Exception as e:
        return {**base_result, "status": "error", "error": str(e)}


# Optional receptor prep

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


# Flat PDBQT finalization for QuickVina2-GPU

def finalize_flat_pdbqt_dir(
    manifest_path: Path,
    flat_dir: Path,
    mode: str = "hardlink",
    overwrite: bool = False,
):
    """
    Create a flat directory containing all successfully prepared ligand PDBQT files.

    Recommended mode:
      hardlink  - fastest and uses no extra data blocks, but requires same filesystem

    Other modes:
      symlink   - very fast, but depends on original shard files remaining in place
      copy      - portable, but duplicates storage and is slower

    This reads ligand_manifest.tsv and only includes rows with status ok/skipped.
    """

    if mode not in {"hardlink", "symlink", "copy"}:
        raise ValueError(f"Unsupported flat PDBQT mode: {mode}")

    flat_dir.mkdir(parents=True, exist_ok=True)

    entries = []
    seen_dest_names = {}

    with open(manifest_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            if row.get("status") not in {"ok", "skipped"}:
                continue

            src = Path(row["out_pdbqt"])
            dest_name = row["ligand_file"]
            dest = flat_dir / dest_name

            if dest_name in seen_dest_names and seen_dest_names[dest_name] != src:
                raise RuntimeError(
                    "Flat-directory filename collision detected:\n"
                    f"  filename: {dest_name}\n"
                    f"  first:    {seen_dest_names[dest_name]}\n"
                    f"  second:   {src}\n"
                    "This would make QuickVina2-GPU input ambiguous. "
                    "Fix ligand_id generation or include a stronger unique component."
                )

            seen_dest_names[dest_name] = src
            entries.append((src, dest))

    linked = 0
    skipped_existing = 0
    missing = 0

    for src, dest in entries:
        if not src.exists():
            missing += 1
            continue

        if dest.exists() or dest.is_symlink():
            try:
                if dest.samefile(src):
                    skipped_existing += 1
                    continue
            except OSError:
                pass

            if not overwrite:
                raise RuntimeError(
                    "Destination already exists and is not the same file:\n"
                    f"  destination: {dest}\n"
                    f"  source:      {src}\n"
                    "Use --flat-pdbqt-overwrite if you intentionally want to replace it."
                )

            dest.unlink()

        if mode == "hardlink":
            os.link(src, dest)
        elif mode == "symlink":
            dest.symlink_to(src.resolve())
        elif mode == "copy":
            shutil.copy2(src, dest)

        linked += 1

    return {
        "flat_dir": str(flat_dir),
        "mode": mode,
        "created": linked,
        "already_present": skipped_existing,
        "missing_sources": missing,
        "expected": len(entries),
    }


# Main
def main():
    parser = argparse.ArgumentParser(
        description="Scalable ligand preparation for parquet-based virtual screening libraries."
    )
    parser.add_argument(
        "input_path",
        type=str,
        help="Path to either a single .parquet file or a directory containing one or more .parquet files.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        required=True,
        help="Output directory for prepared ligands and manifests.",
    )
    parser.add_argument(
        "--flat-pdbqt-dir",
        type=str,
        default=None,
        help = "Optional final flat ligand pdbqt directory for Qvina-GPU usage"
    )
    parser.add_argument(
        "--flat-mode",
        choices=["hardlink", "symlink", "copy"],
        default="hardlink",
        help="How to populate --flat-pdbqt-dir. hardlink is recommended when source and destination are on the same filesystem."
    )

    parser.add_argument(
        "--flat-pdbqt-overwrite",
        action="store_true",
        help="Overwrite existing files in --flat-pdbqt-dir if they are not already linked to the source.",
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
        "--parquet-batch-rows",
        type=int,
        default=5000,
        help="Number of parquet rows to stream per batch when using pyarrow.",
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

    if input_path.is_file():
        if input_path.suffix.lower() not in {".parquet", ".pq"}:
            parser.error(f"If input_path is a file, it must be a .parquet/.pq file: {input_path}")
        parquet_files = [input_path]
    elif input_path.is_dir():
        parquet_files = sorted(input_path.glob("*.parquet")) + sorted(input_path.glob("*.pq"))
        # de-duplicate in case both globs overlap strangely
        parquet_files = sorted({p.resolve() for p in parquet_files})
        parquet_files = [Path(p) for p in parquet_files]
        if not parquet_files:
            parser.error(f"No .parquet or .pq files found in directory: {input_path}")
    else:
        parser.error(f"Input path does not exist: {input_path}")
    total_rows = count_parquet_rows(parquet_files)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    shards_root = outdir / "shards"
    shards_root.mkdir(parents=True, exist_ok=True)

    # Optional receptor prep
    if args.receptor:
        
        receptor_path = Path(args.receptor)
        receptor_out = outdir / f"{receptor_path.stem}_rec.pdbqt"
        if receptor_out.exists():
         pass
      #   print(f"[receptor] exists: {receptor_out}")
        else:
            #print(f"[receptor] preparing: {receptor_path.name} -> {receptor_out.name}")
            prepare_receptor_legacy(receptor_path, receptor_out, del_water=args.del_water)
            #print(f"[receptor] done: {receptor_out}")

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
            "run_id",
            "route_id",
            "seed_id",
            "synthon_id",
            "reaction_id",
            "reaction_name",
            "synthon_smiles",
            "product_smiles",
            "ligand_id",
            "ligand_file",
            "shard",
            "out_pdbqt",
            "duplicate",
            "status",
        ]
    )
    error_writer.writerow(
        [
            "global_idx",
            "source_file",
            "source_idx",
            "run_id",
            "route_id",
            "seed_id",
            "synthon_id",
            "reaction_id",
            "reaction_name",
            "synthon_smiles",
            "product_smiles",
            "shard",
            "out_pdbqt",
            "error",
        ]
    )

    total_seen = 0
    total_ok = 0
    total_skipped = 0
    total_error = 0

    mp_ctx = mp.get_context("spawn")
    # typey type type stall anbd let it go
    pbar = tqdm(
            total=total_rows,
            desc="Preparing Ligands",
            unit="Lig",
            dynamic_ncols = True,
            )

    seen_products = {}
    with ProcessPoolExecutor(max_workers=args.nproc, mp_context=mp_ctx) as ex:
        for parquet_file in parquet_files:
            #print(f"[input] streaming {parquet_file.name}")
            # here need to change out_pdbqt file path
            stream = iter_parquet_records(parquet_file, batch_rows=args.parquet_batch_rows)

            for batch in batched(stream, args.submit_batch):
                tasks = []
                for source_idx, rec in batch:
                    global_idx = total_seen
                    shard_idx = global_idx // args.shard_size
                    shard_name = f"shard_{shard_idx:05d}"
                    shard_dir = shards_root / shard_name / "ligands"

                    ligand_id = ligand_stub_from_record(
                        source_file=parquet_file.name,
                        run_id=normalize_text(rec.get("run_id")),
                        route_id=normalize_text(rec.get("route_id")),
                    )
                    out_pdbqt = shard_dir / f"{ligand_id}.pdbqt"
                        
                    
                    product_smiles = normalize_text(rec.get("product_smiles"))
                    product_valid = normalize_bool(rec.get("product_valid"))
                    

                    try:
                        if not product_valid:
                            raise ValueError(normalize_text(rec.get("failure_reason")) or "product_valid is False")

                        if not product_smiles:
                            raise ValueError("Missing product_smiles")


                        mol, canonical_smiles, num_removed, dedupe_key = clean_marked_product_smiles(product_smiles)
                        
                    except Exception as e:
                        
                        manifest_writer.writerow(
                            [global_idx,
                            parquet_file.name,
                            source_idx,
                            normalize_text(rec.get("run_id")),
                            normalize_text(rec.get("route_id")),
                            normalize_text(rec.get("seed_id")),
                            normalize_text(rec.get("synthon_id")),
                            normalize_text(rec.get("reaction_id")),
                            normalize_text(rec.get("reaction_name")),
                            normalize_text(rec.get("synthon_smiles")),
                            product_smiles,
                            ligand_id,
                            "",
                            "",
                            "",
                            "",
                            "error",
                            
                            
                        ]
                        )

                        error_writer.writerow(
                            [
                                global_idx,
                                parquet_file.name,
                                source_idx,
                                normalize_text(rec.get("run_id")),
                                normalize_text(rec.get("route_id")),
                                normalize_text(rec.get("seed_id")),
                                normalize_text(rec.get("synthon_id")),
                                normalize_text(rec.get("reaction_id")),
                                normalize_text(rec.get("reaction_name")),
                                normalize_text(rec.get("synthon_smiles")),
                                product_smiles,
                                shard_name,
                                "",
                                str(e),
                            ]
                        )

                        total_error += 1
                        total_seen += 1
                        pbar.update(1)
                        continue

                    if dedupe_key in seen_products:
                        master = seen_products[dedupe_key]

                        manifest_writer.writerow(
                            [
                                global_idx,
                                parquet_file.name,
                                source_idx,
                                normalize_text(rec.get("run_id")),
                                normalize_text(rec.get("route_id")),
                                normalize_text(rec.get("seed_id")),
                                normalize_text(rec.get("synthon_id")),
                                normalize_text(rec.get("reaction_id")),
                                normalize_text(rec.get("reaction_name")),
                                normalize_text(rec.get("synthon_smiles")),
                                product_smiles,
                                ligand_id,
                                Path(master["out_pdbqt"]).name,
                                master["shard_name"],
                                master["out_pdbqt"],
                                True,
                                "ok",
                            ]
                        )

                        total_skipped += 1
                        total_seen += 1
                        pbar.update(1)
                        continue

                    out_pdbqt = shard_dir / f"{ligand_id}.pdbqt"

                    seen_products[dedupe_key] = {
                        "ligand_id": ligand_id,
                        "out_pdbqt": str(out_pdbqt),
                        "shard_name": shard_name,
                        "canonical_smiles": canonical_smiles,
                    }

                    tasks.append(
                        {
                            "source_file": parquet_file.name,
                            "source_idx": source_idx,
                            "global_idx": global_idx,
                            "record": rec,
                            "ligand_id": ligand_id,
                            "out_pdbqt": str(out_pdbqt),
                            "shard_name": shard_name,
                            "mol": mol,
                            "canonical_prod_smiles": canonical_smiles,
                            "prepared_ligand_id": ligand_id,
                            "duplicate_of_ligand_id": "",
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
                                result["run_id"],
                                result["route_id"],
                                result["seed_id"],
                                result["synthon_id"],
                                result["reaction_id"],
                                result["reaction_name"],
                                result["synthon_smiles"],
                                result["product_smiles"],
                                result["ligand_id"],
                                result["ligand_file"],
                                result["shard"],
                                result["out_pdbqt"],
                                result["duplicate"],
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
                                result["run_id"],
                                result["route_id"],
                                result["seed_id"],
                                result["synthon_id"],
                                result["reaction_id"],
                                result["reaction_name"],
                                result["synthon_smiles"],
                                result["product_smiles"],
                                result["shard"],
                                result["out_pdbqt"],
                                result["duplicate"],
                                result["error"],
                            ]
                        )
                        

                    pbar.update(1)
                    pbar.set_postfix(
                            ok=total_ok,
                            skipped=total_skipped,
                            error=total_error,
                            )
                manifest_fh.flush()
                errors_fh.flush()
             #   print(
              #      f"[progress] seen={total_seen:,} ok={total_ok:,} "
               #     f"skipped={total_skipped:,} error={total_error:,}"
                #)
    pbar.close()
    manifest_fh.close()
    errors_fh.close()

    
    flat_summary = None
    if args.flat_pdbqt_dir:
        flat_summary = finalize_flat_pdbqt_dir(
            manifest_path=manifest_path,
            flat_dir=Path(args.flat_pdbqt_dir),
            mode=args.flat_mode,
            overwrite=args.flat_pdbqt_overwrite,
        )


    print("\nDone.")
    print(f"Prepared/kept ligands: {total_ok + total_skipped:,}")
    print(f"Newly prepared: {total_ok:,}")
    print(f"Skipped (existing): {total_skipped:,}")
    print(f"Errors: {total_error:,}")
    print(f"Manifest: {manifest_path}")
    print(f"Errors: {errors_path}")
    print(f"Shards root: {shards_root}")

    
    if flat_summary:
        print(f"Flat PDBQT dir: {flat_summary['flat_dir']}")
        print(f"Mode: {flat_summary['mode']}")
        print(f"Created: {flat_summary['created']:,}")
        print(f"Already present: {flat_summary['already_present']:,}")
        print(f"Missing sources: {flat_summary['missing_sources']:,}") 
        print(f"Expected inputs: {flat_summary['expected']:,}")


if __name__ == "__main__":
    main()

