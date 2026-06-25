import time
import argparse
from argparse import ArgumentParser
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, QED, rdMolDescriptors, DataStructs
from pathlib import Path
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from pharmacophore import get_feature_factory, build_pharmacophore, compute_pharmacophore_score
from utilities import *
# Argument Utilities

def valid_path(path_str):
    path = Path(path_str)
    if not path.exists():
        raise argparse.ArgumentTypeError(f"Path does not exist: {path}")
    return path

# Molecule Utilities

def manifest_pose_keys(row):
    """
    Generate likely lookup keys for one manifest row.
    Primary match should be out_pdbqt filename stem.
    """
    keys = []

    if "out_pdbqt" in row and pd.notna(row["out_pdbqt"]):
        keys.append(Path(str(row["out_pdbqt"])).stem)

    for col in (
        "Name",
        "name",
        "CanonicalProductSMILES",
        "canonical_product_smiles",
        "clean_smiles",
        "ProductSMILES",
        "product_smiles",
    ):
        if col in row and pd.notna(row[col]):
            keys.append(str(row[col]).strip())

    if "product_smiles" in row and pd.notna(row["product_smiles"]):
        try:
            _, clean_smi, _ = clean_marked_product_smiles(str(row["product_smiles"]))
            keys.append(clean_smi)
        except Exception:
            pass

    return [k for k in keys if k]


def find_pose_for_manifest_row(row, pose_lookup):
    """Return matching docked pose mol and the key used."""
    for key in manifest_pose_keys(row):
        mol = pose_lookup.get(key)
        if mol is not None:
            return mol, key

    return None, None
def clean_marked_product_smiles(product_smiles: str, placeholder_symbols=("W","V","Rf", "Db", "Sg", "Rh5",), placeholder_atomic_nums=(0,)):
    """
    Convert an enumerated product SMILES containing residual synthon marks / placeholders
    into a dockable RDKit molecule.

    Current behavior:
    - removes placeholder atoms such as terminal [W]
    - removes atom-map numbers (e.g. :20)
    - repairs aromatic [n]([W]) -> [nH] when the placeholder was acting as a cap

    Returns:
      (clean_mol, clean_smiles, removed_placeholder_count)
    """
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

    return mol, clean_smiles, len(remove_idxs)

def clean_products_from_manifest(df):
    """Clean product_smiles once and return aligned clean_smiles/product_mols lists."""
    clean_smiles = []
    product_mols = []

    for smi in df["product_smiles"]:
        if pd.isna(smi):
            clean_smiles.append(None)
            product_mols.append(None)
            continue

        try:
            m, clean_smi, _ = clean_marked_product_smiles(smi)
        except Exception:
            m, clean_smi = None, None

        clean_smiles.append(clean_smi if m is not None else None)
        product_mols.append(m)

    return clean_smiles, product_mols

def smiles_to_mol(smi_path):
    mols = []
    with open(smi_path) as inp:
        for line in inp:
            sline = line.strip()
            if sline:
                for smiles in line.split()[0].split("."):
                    mols.append(Chem.MolFromSmiles(smiles))
    if not mols:
        raise ValueError("No valid mols in smi file")
    return mols

def sdf_to_mol(sdf_path):
    supplier = Chem.SDMolSupplier(sdf_path)
    mols = [m for m in supplier if m is not None]
    if not mols:
        raise ValueError("No valid mols in sdf file")
    return mols

def compute_fingerprint(mol, radius=2, fp_size=2048):
    morgan_gen = GetMorganGenerator(radius=radius, fpSize=fp_size)
    return morgan_gen.GetFingerprint(mol)

# 2D Similarity scoring function

def compute_2d_similarity(ref_fp, mol):
    if mol is None:
        return None
    fp = compute_fingerprint(mol)
    return DataStructs.TanimotoSimilarity(ref_fp, fp)

# Core Scoring Pipeline

def score_manifest(args):

    # Load manifest
    df = pd.read_csv(args.manifest, sep="\t")
    if args.dedupe:
        df = df[df.duplicate == False]
    
    #df = df.head(10000)
    # Prepare reference
    if any(item == ".sdf" for item in args.reference.suffixes):
        try:
            ref_mols = sdf_to_mol(args.reference)
        except OSError as e:
            raise OSError(f"Could not read sdf path {args.refernece}")
            
    elif any(item == ".smi" for item in args.reference.suffixes):
        ref_mols = smiles_to_mol(args.reference)
    
    if ref_mols is None:
        raise ValueError("Could Not Read reference mol(s).")

    # 2D scoring
    
    if "2d" in args.method: # simple tanimoto similarity score
        print("Starting 2D Similarity scoring\nMay take a couple minutes depending on the size of your dataset")
        start_time = time.time()
        refs = [compute_fingerprint(ref_mol) for ref_mol in ref_mols]
        
        clean_smiles = []
        fps = []
        clogp_list = []
        mw_list = []
        tpsa_list = []
        qed_list = []

        for smi in df["product_smiles"]:
            if smi is None:
                fps.append(None)
                clogp_list.append(None)
                mw_list.append(None)
                tpsa_list.append(None)
                qed_list.append(None)

            else:
                m, clean_smi, _ = clean_marked_product_smiles(smi)
                #m = Chem.MolFromSmiles(smi)
                if m:
                    # descriptors
                    clean_smiles.append(clean_smi)
                    fp = compute_fingerprint(m)
                    fps.append(fp)
                    clogp_list.append(Crippen.MolLogP(m))
                    mw_list.append(Descriptors.MolWt(m))
                    tpsa_list.append(rdMolDescriptors.CalcTPSA(m))
                    qed_list.append(QED.qed(m))
                else:
                    clean_smiles.append(None)
                    fps.append(None)
                    clogp_list.append(None)
                    mw_list.append(None)
                    tpsa_list.append(None)
                    qed_list.append(None)
        df["clean_smiles"] = clean_smiles
        valid_fps = [fp for fp in fps if fp is not None]

        all_sims = [
            DataStructs.BulkTanimotoSimilarity(ref_fp, valid_fps)
            for ref_fp in refs
        ]
        max_sims = [max(vals) for vals in zip(*all_sims)]
        #sims = max(DataStructs.BulkTanimotoSimilarity(ref_fp, valid_fps) for ref_fp in refs)

        it = iter(max_sims)
        df["tanimoto_sim"] = [
            next(it) if fp is not None else None
            for fp in fps
            ]
        
        df["clogp"] = clogp_list
        df["mw"] = mw_list
        df["tpsa"] = tpsa_list
        df["qed"] = qed_list

        end_time = time.time()
        total_time = end_time - start_time
        print(f"2D scoring Complete.\nTotalTime: {round(total_time, 3)} for {len(df)} samples")
        #print(f"avg per 100\n{(1000/total_time) * 100}")
    
    # 3D placeholder
    
    if "3d" in args.method:
        df["score_3d_similarity"] = None

    # Pharmacophore
    if "pharmacophore" in args.method:
        print("Starting pharmacophore scoring")
        start_time = time.time()

        if not all(has_3d_conformer(m) for m in ref_mols):
            raise ValueError(
                "Pharmacophore scoring requires 3D reference molecule(s). "
                "Use a 3D SDF reference; SMILES references do not carry coordinates."
            )

        if args.poses_sdf is None:
            raise ValueError(
                "Pharmacophore scoring requires a combined SDF of top docked poses. "
                "Pass it with --poses_sdf /path/to/combined.sdf."
            )

        factory = get_feature_factory()
        pharm_feats = build_pharmacophore(ref_mols, factory)

        pose_lookup, pose_scores = load_pose_sdf_lookup(args.poses_sdf)

        pharm_scores = []
        docking_scores = []
        pose_keys = []
        missing_pose_count = 0

        for _, row in df.iterrows():
            pose_mol, pose_key = find_pose_for_manifest_row(row, pose_lookup)

            if pose_mol is None or not has_3d_conformer(pose_mol):
                pharm_scores.append(None)
                docking_scores.append(None)
                pose_keys.append(None)
                missing_pose_count += 1
                continue

            pharm_scores.append(
                compute_pharmacophore_score(pharm_feats, pose_mol, factory)
            )
            docking_scores.append(pose_scores.get(pose_key))
            pose_keys.append(pose_key)

        df["score_pharmacophore"] = pharm_scores
        df["docking_score"] = docking_scores
        df["pharmacophore_pose_key"] = pose_keys

        end_time = time.time()
        total_time = end_time - start_time

        print(
            f"Pharmacophore scoring complete.\n"
            f"TotalTime: {round(total_time, 3)} for {len(df)} samples\n"
            f"Missing poses: {missing_pose_count}"
        )

    # Save updated manifest
    if args.output_path is None:
        output_path = args.manifest.with_name(f"{args.manifest.stem}_scored{args.manifest.suffix}")

    #df.to_csv(output_path, sep="\t", index=False) cahn

    return df


if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument("--manifest",
                        required = True,
                        type = valid_path,
                        help = "path to manifest.tsv to update records."
                        )
    parser.add_argument("--reference",
                        required = True,
                        type = valid_path,
                        help = "path to reference ligand for similarity and strucutre matching. SMI/SDF"
                        )
    parser.add_argument("--method",
                        nargs = "+",
                        choices = ['2d', '3d', 'pharmacophore', 'docking'], 
                        help = "Scoring Method, options include ['2d', '3d', 'pharmacophore', 'docking'].",
                        required = True
                        )
    
    parser.add_argument("--dedupe", 
                        action= "store_true",
                        help = "manuallly ignore duplicate products in manifest if not already"
                        )

    parser.add_argument("--output_path",
                        "-o",
                        default = None
                        )
    parser.add_argument("--poses_sdf",
                        type=valid_path,
                        default=None,
                        help="combined SDF containing top docked poses for pharmacophore scoring"
                    )
    args = parser.parse_args()
   



    df = score_manifest(args)
    print(df.head(15))
    print(df.clean_smiles)
    print("Scoring complete.")

