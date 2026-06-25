import os
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import os
import numpy as np


def get_feature_factory():
    # Build from RDkit default features
    fdef = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    return ChemicalFeatures.BuildFeatureFactory(fdef)

g
def extract_features(mol, factory):
    feats = factory.GetFeaturesForMol(mol)

    feature_list = []
    for f in feats:
        pos = f.GetPos()
        feature_list.append((
            f.GetFamily(),  # Donor, Acceptor, Aromatic, Hydrophobe, etc. or not even mention if npt brought up otherwise
            (pos.x, pos.y, pos.z)
        ))

    return feature_list

def build_pharmacophore(ref_mols_3d, factory):
    # extracts features from mol and builds a list of features and coordinates
    pharm_features = []

    for mol in ref_mols_3d:
        feats = extract_features(mol, factory)
        pharm_features.extend(feats)

    return pharm_features

def distance(p1, p2):
    return np.linalg.norm(np.array(p1) - np.array(p2))

def match_features(pharm_feats, ligand_feats, tolerance=1.5):
    # takes two features lists and counts overlaps of features within <tolerance> distanc from each other.
    matches = 0

    for p_type, p_pos in pharm_feats:

        for l_type, l_pos in ligand_feats:
            if p_type != l_type:
                continue

            if distance(p_pos, l_pos) <= tolerance:
                matches += 1
                break
    return matches


def compute_pharmacophore_score(pharm_feats, ligand_mol, factory):
    # makes a standardized/normalized score irrelevant of number of features
    if ligand_mol is None:
        return None

    ligand_feats = extract_features(ligand_mol, factory)

    if not ligand_feats:
        return 0.0

    matches = match_features(pharm_feats, ligand_feats)

    score = matches / len(pharm_feats) if pharm_feats else 0.0

    return score