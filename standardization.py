
"""Seed standardization utilities.

Phase 1 scope:
- Canonicalize seed SMILES using RDKit canonical SMILES while preserving stereochemistry.
- Provide a lightweight StandardizedSeed record.

Notes:
- Site IDs are defined on the canonicalized seed only.
- Stereochemistry is preserved (isomericSmiles=True).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

from rdkit import Chem


@dataclass(frozen=True)
class StandardizedSeed:
    """Canonical representation of a seed molecule."""

    seed_input_smiles: str
    seed_canonical_smiles: str
    seed_id: str


def canonicalize_seed_smiles(
    seed_smiles: str,
    seed_id: Optional[str] = None,
) -> Optional[Tuple[StandardizedSeed, Chem.Mol]]:
    """Canonicalize a seed SMILES string and return (StandardizedSeed, RDKit Mol).

    Preserves stereochemistry.

    Returns None if the SMILES cannot be parsed.
    """
    if seed_smiles is None:
        return None

    mol = Chem.MolFromSmiles(seed_smiles)
    if mol is None:
        return None

    # Canonical SMILES with stereochemistry preserved.
    can = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

    # Re-parse canonical SMILES to ensure atom ordering corresponds to canonical form.
    can_mol = Chem.MolFromSmiles(can)
    if can_mol is None:
        return None

    sid = seed_id if seed_id is not None else can

    return StandardizedSeed(seed_input_smiles=seed_smiles, seed_canonical_smiles=can, seed_id=str(sid)), can_mol
