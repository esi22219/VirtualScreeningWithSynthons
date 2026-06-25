import re
import rdkit
from rdkit import Chem
from pathlib import Path


def parse_vina_score_from_remark(remark):
    """Extract Vina affinity from an SDF REMARK field."""
    if not remark:
        return None

    match = re.search(r"VINA RESULT:\s*([-+]?\d*\.?\d+)", str(remark))
    return float(match.group(1)) if match else None


def sdf_prop(mol, prop_name):
    """Safely get an SDF property."""
    return mol.GetProp(prop_name).strip() if mol is not None and mol.HasProp(prop_name) else None


def load_pose_sdf_lookup(poses_sdf_path):
    """
    Load combined.sdf once and build lookup dictionaries.

    Returns:
      pose_lookup: key -> RDKit mol with 3D coordinates
      pose_scores: key -> Vina score, if found
    """
    supplier = Chem.SDMolSupplier(str(poses_sdf_path), removeHs=False)

    pose_lookup = {}
    pose_scores = {}

    for mol in supplier:
        if mol is None:
            continue

        remark = sdf_prop(mol, "REMARK")
        vina_score = parse_vina_score_from_remark(remark)

        keys = []

        # RDKit molecule name, if present
        name = sdf_prop(mol, "_Name")
        if name:
            keys.append(name)

        # Direct SDF props, if they exist as separate fields
        for prop in ("Name", "CanonicalProductSMILES", "ProductSMILES", "SMILES"):
            value = sdf_prop(mol, prop)
            if value:
                keys.append(value)

        # Your example stores key/value metadata inside REMARK
        if remark:
            for line in str(remark).splitlines():
                line = line.strip()
                if " = " not in line:
                    continue

                k, v = line.split(" = ", 1)

                if k in {"Name", "CanonicalProductSMILES", "ProductSMILES", "SMILES"}:
                    keys.append(v.strip())

        for key in keys:
            key = str(key).strip()
            pose_lookup.setdefault(key, mol)

            if vina_score is not None:
                pose_scores.setdefault(key, vina_score)

    if not pose_lookup:
        raise ValueError(f"No usable docked poses found in SDF: {poses_sdf_path}")

    return pose_lookup, pose_scores


def has_3d_conformer(mol):
    """Return True when an RDKit mol has at least one conformer with coordinates."""
    return mol is not None and mol.GetNumConformers() > 0