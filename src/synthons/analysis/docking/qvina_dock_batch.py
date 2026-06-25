#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path

import numpy as np
from rdkit import Chem


DEFAULT_BASE = Path("/home2/esi22219/p53_docking_lib")
DEFAULT_APP = Path("/home2/esi22219/apps/Vina-GPU-2.1/QuickVina2-GPU-2.1/QuickVina2-GPU-2-1")


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def compute_center_and_size(sdf_file: Path, buffer: float, min_size: float):
    """
    Read first valid molecule from reference SDF and compute:
      center = mean atom coordinate
      size   = bounding-box dimension + 2 * buffer, with min_size floor
    """
    sdf_file = Path(sdf_file)

    if not sdf_file.is_file():
        raise FileNotFoundError(f"Reference SDF does not exist: {sdf_file}")

    suppl = Chem.SDMolSupplier(str(sdf_file), sanitize=False)

    if suppl is None:
        raise FileNotFoundError(f"Could not open SDF file: {sdf_file}")

    first_mol = None
    for mol in suppl:
        if mol is not None:
            first_mol = mol
            break

    if first_mol is None:
        raise ValueError(f"No valid molecules found in SDF file: {sdf_file}")

    if first_mol.GetNumConformers() == 0:
        raise ValueError(f"Reference SDF molecule has no conformer coordinates: {sdf_file}")

    coords = first_mol.GetConformer().GetPositions()

    center_x, center_y, center_z = coords.mean(axis=0)
    max_x, max_y, max_z = coords.max(axis=0)
    min_x, min_y, min_z = coords.min(axis=0)

    size_x = max((max_x - min_x) + 2 * buffer, min_size)
    size_y = max((max_y - min_y) + 2 * buffer, min_size)
    size_z = max((max_z - min_z) + 2 * buffer, min_size)

    center = np.array([center_x, center_y, center_z], dtype=float)
    size = np.array([size_x, size_y, size_z], dtype=float)

    return center, size


def write_config(
    config_path: Path,
    receptor: Path,
    ligand_dir: Path,
    opencl_binary_path: Path,
    center,
    size,
    threads: int,
    search_depth=None,
):
    """
    Write QuickVina2-GPU config for a single shard.
    """
    config_path = Path(config_path)
    config_path.parent.mkdir(parents=True, exist_ok=True)

    with open(config_path, "w") as f:
        f.write(f"receptor = {receptor}\n")
        f.write(f"ligand_directory = {ligand_dir}\n")
        f.write(f"opencl_binary_path = {opencl_binary_path}\n")
        f.write(f"center_x = {center[0]:.3f}\n")
        f.write(f"center_y = {center[1]:.3f}\n")
        f.write(f"center_z = {center[2]:.3f}\n")
        f.write(f"size_x = {size[0]:.3f}\n")
        f.write(f"size_y = {size[1]:.3f}\n")
        f.write(f"size_z = {size[2]:.3f}\n")
        f.write(f"thread = {threads}\n")

        if search_depth is not None:
            f.write(f"search_depth = {search_depth}\n")


def count_pdbqt(ligand_dir: Path) -> int:
    if not ligand_dir.is_dir():
        return 0
    return sum(1 for _ in ligand_dir.glob("*.pdbqt"))


def discover_shards(base: Path, start_shard=None, end_shard=None):
    shards_dir = base / "shards"
    shards = sorted([p for p in shards_dir.glob("shard_*") if p.is_dir()])

    filtered = []
    for shard in shards:
        shard_id = shard.name.replace("shard_", "")

        if start_shard is not None and shard_id < start_shard:
            continue

        if end_shard is not None and shard_id > end_shard:
            continue

        filtered.append(shard)

    return filtered


def make_env(gpu_id: str):
    """
    Best-effort GPU isolation.

    QuickVina2-GPU uses OpenCL and has no runtime device argument, so this may or
    may not force separate GPUs depending on how the OpenCL runtime enumerates devices.

    Test with:
      watch -n 2 nvidia-smi
    """
    env = os.environ.copy()

    env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    env["GPU_DEVICE_ORDINAL"] = str(gpu_id)
    env["ROCR_VISIBLE_DEVICES"] = str(gpu_id)
    env["HIP_VISIBLE_DEVICES"] = str(gpu_id)

    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("NUMEXPR_NUM_THREADS", "1")

    return env


def run_shard(
    shard_dir: Path,
    gpu_id: str,
    args,
    center,
    size,
):
    base = args.base
    app = args.app

    shard_name = shard_dir.name
    shard_id = shard_name.replace("shard_", "")

    ligand_dir = shard_dir / "ligands"
    config_path = base / "configs" / f"config_{shard_id}.txt"

    log_dir = base / "logs"
    status_dir = base / "run_status" / shard_name

    log_dir.mkdir(parents=True, exist_ok=True)
    status_dir.mkdir(parents=True, exist_ok=True)

    log_path = log_dir / f"{shard_name}.gpu{gpu_id}.log"
    done_file = status_dir / "DONE"
    failed_file = status_dir / "FAILED"
    status_json = status_dir / "status.json"

    #ligand_count = count_pdbqt(ligand_dir)

    status = {
        "shard": shard_name,
        "shard_id": shard_id,
        "gpu_id": gpu_id,
        "ligand_dir": str(ligand_dir),
        #"ligand_count": ligand_count, ignore
        "config": str(config_path),
        "log": str(log_path),
        "start_time": now(),
        "end_time": None,
        "state": None,
        "returncode": None,
    }

    if done_file.exists() and not args.force:
        status["state"] = "SKIPPED_DONE"
        status["returncode"] = 0
        status["end_time"] = now()
        status_json.write_text(json.dumps(status, indent=2))
        return status

    if not ligand_dir.is_dir():
        status["state"] = "FAILED_MISSING_LIGAND_DIR"
        status["returncode"] = 2
        status["end_time"] = now()
        failed_file.touch()
        status_json.write_text(json.dumps(status, indent=2))
        return status

    #if ligand_count == 0:
     #   status["state"] = "FAILED_NO_LIGANDS"
      #  status["returncode"] = 3
       # status["end_time"] = now()
        #failed_file.touch()
        #status_json.write_text(json.dumps(status, indent=2))
        #return status
        
    # potential failure reasons
    if not args.receptor.is_file(): 
        status["state"] = "FAILED_MISSING_RECEPTOR"
        status["returncode"] = 4
        status["end_time"] = now()
        failed_file.touch()
        status_json.write_text(json.dumps(status, indent=2))
        return status

    if not app.is_file():
        status["state"] = "FAILED_MISSING_APP"
        status["returncode"] = 5
        status["end_time"] = now()
        failed_file.touch()
        status_json.write_text(json.dumps(status, indent=2))
        return status

    # Generate config fresh for this shard.
    write_config(
        config_path=config_path,
        receptor=args.receptor,
        ligand_dir=ligand_dir,
        opencl_binary_path=args.opencl_binary_path,
        center=center,
        size=size,
        threads=args.threads,
        search_depth=args.search_depth,
    )

    cmd = [str(app), "--config", str(config_path)]
    env = make_env(gpu_id)

    status["command"] = cmd
    status["center"] = [float(x) for x in center]
    status["size"] = [float(x) for x in size]
    status["threads"] = args.threads
    status["search_depth"] = args.search_depth

    if args.dry_run:
        status["state"] = "DRY_RUN"
        status["returncode"] = 0
        status["end_time"] = now()
        status_json.write_text(json.dumps(status, indent=2))
        print(f"[DRY RUN] {shard_name} GPU {gpu_id}: {' '.join(cmd)}")
        return status

    failed_file.unlink(missing_ok=True)

    with open(log_path, "w") as log:
        log.write("=" * 100 + "\n")
        log.write(f"Shard: {shard_name}\n")
        log.write(f"Assigned GPU attempt: {gpu_id}\n")
        log.write(f"Start: {status['start_time']}\n")
        log.write(f"Host: {os.uname().nodename}\n")
        log.write(f"Ligand dir: {ligand_dir}\n")
        log.write(f"Config: {config_path}\n")
        log.write(f"Receptor: {args.receptor}\n")
        log.write(f"Reference SDF: {args.reference_sdf}\n")
        log.write(f"Center: {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}\n")
        log.write(f"Size: {size[0]:.3f}, {size[1]:.3f}, {size[2]:.3f}\n")
        log.write(f"Threads: {args.threads}\n")
        log.write(f"Search depth: {args.search_depth}\n")
        log.write(f"OpenCL binary path: {args.opencl_binary_path}\n")
        log.write(f"Executable: {app}\n")
        log.write(f"Command: {' '.join(cmd)}\n")
        log.write(f"CUDA_VISIBLE_DEVICES={env.get('CUDA_VISIBLE_DEVICES')}\n")
        log.write("=" * 100 + "\n\n")
        log.flush()

        log.write("Generated config:\n")
        log.write("-" * 100 + "\n")
        log.write(config_path.read_text())
        log.write("-" * 100 + "\n\n")
        log.flush()

        log.write("nvidia-smi before run:\n")
        log.write("-" * 100 + "\n")
        log.flush()

        subprocess.run(
            ["nvidia-smi"],
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
            check=False,
        )

        log.write("\n" + "=" * 100 + "\n")
        log.write("QuickVina2-GPU output:\n")
        log.write("=" * 100 + "\n")
        log.flush()

        try:
            result = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
                check=False,
            )
            rc = result.returncode
        except Exception as e:
            log.write("\n")
            log.write("=" * 100 + "\n")
            log.write(f"Python exception while running QuickVina2-GPU: {repr(e)}\n")
            log.write("=" * 100 + "\n")
            rc = 99

        log.write("\n" + "=" * 100 + "\n")
        log.write("nvidia-smi after run:\n")
        log.write("-" * 100 + "\n")
        log.flush()

        subprocess.run(
            ["nvidia-smi"],
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
            check=False,
        )

        log.write("\n" + "=" * 100 + "\n")
        log.write(f"Exit code: {rc}\n")
        log.write(f"End: {now()}\n")
        log.write("=" * 100 + "\n")
        log.flush()

    status["returncode"] = rc
    status["end_time"] = now()

    if rc == 0:
        done_file.touch()
        failed_file.unlink(missing_ok=True)
        status["state"] = "DONE"
    else:
        failed_file.touch()
        status["state"] = "FAILED"

    status_json.write_text(json.dumps(status, indent=2))
    return status


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate QuickVina2-GPU configs from reference SDF and run ligand shards safely."
    )

    p.add_argument("--base", type=Path, default=DEFAULT_BASE)
    p.add_argument("--app", type=Path, default=DEFAULT_APP)

    p.add_argument(
        "--receptor",
        type=Path,
        required=True,
        help="Path to receptor PDBQT.",
    )

    p.add_argument(
        "--reference-sdf",
        type=Path,
        required=True,
        help="Reference ligand SDF used to compute center and box size.",
    )

    p.add_argument(
        "--opencl-binary-path",
        type=Path,
        default=DEFAULT_APP.parent,
        help="Value written as opencl_binary_path in each config.",
    )

    p.add_argument("--buffer", type=float, default=8.0)
    p.add_argument("--min-size", type=float, default=20.0)

    p.add_argument("--threads", type=int, default=1000)
    p.add_argument("--search-depth", type=int, default=None)

    p.add_argument("--gpus", default="0")
    p.add_argument("--parallel", type=int, default=2)

    p.add_argument("--start-shard", default=None, help="Example: 00000")
    p.add_argument("--end-shard", default=None, help="Example: 00099")
    p.add_argument("--limit", type=int, default=None)

    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--force", action="store_true", help="Rerun shards even if DONE exists.")

    return p.parse_args()


def main():
    args = parse_args()

    gpus = [x.strip() for x in args.gpus.split(",") if x.strip()]
    if not gpus:
        raise RuntimeError("No GPU IDs provided. An example would be to use --gpus 0,1 or --gpus 0.")

    if args.parallel > len(gpus):
        print("WARNING: --parallel is greater than the number of GPU IDs.")
        print("This can oversubscribe GPUs.")

    center, size = compute_center_and_size(
        sdf_file=args.reference_sdf,
        buffer=args.buffer,
        min_size=args.min_size,
    )

    shards = discover_shards(
        base=args.base,
        start_shard=args.start_shard,
        end_shard=args.end_shard,
    )

    pending = []
    for shard in shards:
        done_file = args.base / "run_status" / shard.name / "DONE"
        if done_file.exists() and not args.force:
            continue
        pending.append(shard)

    if args.limit is not None:
        pending = pending[: args.limit]

    # prints for readability
    print("=" * 100)
    print("QuickVina2-GPU shard runner")
    print("=" * 100)
    print(f"Base: {args.base}")
    print(f"Executable: {args.app}")
    print(f"Receptor: {args.receptor}")
    print(f"Reference SDF: {args.reference_sdf}")
    print(f"OpenCL binary path in config: {args.opencl_binary_path}")
    print(f"Buffer: {args.buffer}")
    print(f"Min size: {args.min_size}")
    print(f"Computed center: {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}")
    print(f"Computed size: {size[0]:.3f}, {size[1]:.3f}, {size[2]:.3f}")
    print(f"Threads: {args.threads}")
    print(f"Search depth: {args.search_depth}")
    print(f"Discovered shards: {len(shards)}")
    print(f"Pending shards: {len(pending)}")
    print(f"GPU IDs: {gpus}")
    print(f"Parallel jobs: {args.parallel}")
    print(f"Dry run: {args.dry_run}")
    print("=" * 100)

    if not pending:
        print("Nothing to run.")
        return

    max_workers = args.parallel
    available_gpus = gpus[:max_workers]

    active = {}
    results = []
    next_idx = 0

    # start multiprocessing
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for gpu_id in available_gpus:
            if next_idx >= len(pending):
                break

            shard = pending[next_idx]
            next_idx += 1

            print(f"[{now()}] Launching {shard.name} on GPU assignment {gpu_id}")

            future = executor.submit(
                run_shard,
                shard,
                gpu_id,
                args,
                center,
                size,
            )
            active[future] = gpu_id

        while active:
            done_futures, _ = wait(active.keys(), return_when=FIRST_COMPLETED)

            for future in done_futures:
                gpu_id = active.pop(future)

                try:
                    result = future.result()
                except Exception as e:
                    result = {
                        "shard": "UNKNOWN",
                        "gpu_id": gpu_id,
                        "state": "PYTHON_RUNNER_EXCEPTION",
                        "returncode": 100,
                        "error": repr(e),
                    }

                results.append(result)

                print(
                    f"[{now()}] Finished {result.get('shard')} "
                    f"on GPU assignment {gpu_id}: "
                    f"{result.get('state')} rc={result.get('returncode')}"
                )

                if next_idx < len(pending):
                    shard = pending[next_idx]
                    next_idx += 1

                    print(f"[{now()}] Launching {shard.name} on GPU assignment {gpu_id}")

                    new_future = executor.submit(
                        run_shard,
                        shard,
                        gpu_id,
                        args,
                        center,
                        size,
                    )
                    active[new_future] = gpu_id

    summary_path = args.base / "logs" / "run_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(results, indent=2))

    n_done = sum(1 for r in results if r.get("state") == "DONE")
    n_failed = sum(1 for r in results if str(r.get("state", "")).startswith("FAILED"))
    n_skipped = sum(1 for r in results if r.get("state") == "SKIPPED_DONE")

    # print nicely for finish
    print("=" * 100)
    print("Run complete")
    print("=" * 100)
    print(f"DONE this run: {n_done}")
    print(f"FAILED this run: {n_failed}")
    print(f"SKIPPED this run: {n_skipped}")
    print(f"Summary written to: {summary_path}")
    print("=" * 100)


if __name__ == "__main__":
    main()
