import argparse
import subprocess
from pathlib import Path
from rdkit import Chem



def read_sdf_coords(sdf_file):
    suppl = Chem.SDMolSupplier(str(sdf_file), sanitize=False)
    if not suppl:
        raise FileNotFoundError(f"Could not open SDF file: {sdf_path}")

    # Get the first molecule safely
    first_mol = None
    for mol in suppl:
        if mol is not None:  # Skip invalid molecules
            first_mol = mol
            break
    cx, cy, cz = first_mol.GetConformer().GetPositions().mean(0)
    max_x, max_y, max_z = first_mol.GetConformer().GetPositions.max(0)
    min_x, min_y, min_z = first_mol.GetConformer().GetPositions.min(0)
    
def compute_center_and_size(sdf_file, buffer, min_size):
    suppl = Chem.SDMolSupplier(str(sdf_file), sanitize=False)
    if not suppl:
        raise FileNotFoundError(f"Could not open SDF file: {sdf_path}")

    # Get the first molecule safely
    first_mol = None
    for mol in suppl:
        if mol is not None:  # Skip invalid molecules
            first_mol = mol
            break
    center_x, center_y, center_z = first_mol.GetConformer().GetPositions().mean(0)
    max_x, max_y, max_z = first_mol.GetConformer().GetPositions().max(0)
    min_x, min_y, min_z = first_mol.GetConformer().GetPositions().min(0)

    size_x = max((max_x - min_x) + 2 * buffer, min_size)
    size_y = max((max_y - min_y) + 2 * buffer, min_size)
    size_z = max((max_z - min_z) + 2 * buffer, min_size)

    return center_x, center_y, center_z, size_x, size_y, size_z


def write_config(
    config_path,
    receptor,
    ligand_dir,
    qvina_bin,
    center,
    size,
    threads,
):
    with open(config_path, "w") as f:
        f.write(f"receptor = {receptor}\n")
        f.write(f"ligand_directory = {ligand_dir}\n")
        f.write(f"opencl_binary_path = {qvina_bin}\n")
        f.write(f"center_x = {center[0]:.3f}\n")
        f.write(f"center_y = {center[1]:.3f}\n")
        f.write(f"center_z = {center[2]:.3f}\n")
        f.write(f"size_x = {size[0]:.3f}\n")
        f.write(f"size_y = {size[1]:.3f}\n")
        f.write(f"size_z = {size[2]:.3f}\n")
        f.write(f"thread = {threads}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor", required=True)
    parser.add_argument("--ligand_dir", required=True)
    parser.add_argument("--reference_ligand", required=True)
    parser.add_argument("--qvina_bin", required=True)
    parser.add_argument("--buffer", type=float, default=7.0)
    parser.add_argument("--min_size", type=float, default=20.0)
    parser.add_argument("--threads", type=int, default=5000)
    parser.add_argument("--config_prefix", required=True)

    args = parser.parse_args()

    cx, cy, cz, sx, sy, sz = compute_center_and_size(
        args.reference_ligand,
        buffer=args.buffer,
        min_size=args.min_size,
    )
    
    conf_path = f"{args.config_prefix}_config.txt"
    write_config(
        conf_path,
        args.receptor,
        args.ligand_dir,
        args.qvina_bin,
        center=(cx, cy, cz),
        size=(sx, sy, sz),
        threads=args.threads,
    )
    
    binary = args.qvina_bin + "/QuickVina2-GPU-2-1"
    subprocess.run(
        [binary, "--config", conf_path],
        check=True,
    )


if __name__ == "__main__":
    main()
