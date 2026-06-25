from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import logging
import re

from rdkit import Chem

logger = logging.getLogger(__name__)

# class to easily handle and save info regarding reactive sites on a molecule
@dataclass(frozen=True)
class ReactiveSite:
    """A single reactive site on a canonicalized seed."""

    site_id: int
    mark_type: str              # e.g., 'C:10', 'c:10', 'N:20'
    canonical_atom_rank: int
    atom_idx: int
    mark_occurrence_index: int  # ordinal among same mark_type in this seed
    site_label: str             # e.g., 'C:10@0'





def _atom_to_mark_type(atom: Chem.Atom) -> Optional[str]:
    """Convert an RDKit atom with atom-map number to the mark_type string.

    Uses lowercase element symbol for aromatic atoms to match compatibility maps (e.g., 'c:10').
    Returns None if atom has no atom-map number.
    """
    amap = atom.GetAtomMapNum()
    if not amap:
        return None

    sym = atom.GetSymbol()
    if atom.GetIsAromatic():
        sym = sym.lower()

    return f"{sym}:{amap}"


def list_reactive_sites(can_mol: Chem.Mol) -> List[ReactiveSite]:
    """List reactive sites on a canonicalized RDKit Mol.

    Deterministic ordering:
    - Compute canonical atom ranks with breakTies=False.
    - Sort sites by (canonical_rank, atom_idx).
    - Assign site_id sequentially.
    - Compute mark_occurrence_index within each mark_type in the sorted order.
    """
    if can_mol is None:
        return []

    # Canonical ranks: stable on canonicalized molecule; breakTies=False keeps symmetric atoms tied.
    ranks = Chem.CanonicalRankAtoms(can_mol, breakTies=False)

    raw_sites: List[Tuple[int, int, str]] = []  # (rank, atom_idx, mark_type)
    for atom in can_mol.GetAtoms():
        mark = _atom_to_mark_type(atom)
        if mark is None:
            continue
        raw_sites.append((int(ranks[atom.GetIdx()]), atom.GetIdx(), mark))

    # Sort deterministically.
    raw_sites.sort(key=lambda x: (x[0], x[1]))

    # Assign occurrence index per mark_type in sorted order.
    occ: Dict[str, int] = {}
    sites: List[ReactiveSite] = []
    for site_id, (rank, atom_idx, mark) in enumerate(raw_sites):
        k = occ.get(mark, 0)
        occ[mark] = k + 1
        label = f"{mark}@{k}"
        sites.append(
            ReactiveSite(
                site_id=site_id,
                mark_type=mark,
                canonical_atom_rank=rank,
                atom_idx=atom_idx,
                mark_occurrence_index=k,
                site_label=label,
            )
        )

    return sites


def list_sites_pretty(can_mol: Chem.Mol) -> List[dict]:
    """Return a JSON-friendly listing of sites (useful for CLI/logging)."""
    return [
        {
            "site_id": s.site_id,
            "mark_type": s.mark_type,
            "site_label": s.site_label,
            "canonical_atom_rank": s.canonical_atom_rank,
            "atom_idx": s.atom_idx,
            "mark_occurrence_index": s.mark_occurrence_index,
        }
        for s in list_reactive_sites(can_mol)
    ]


def resolve_allowed_sites(
    sites: Sequence[ReactiveSite],
    allowed_sites: Optional[Sequence[int]] = None,
    allowed_site_specs: Optional[Sequence[str]] = None,
    invalid_policy: str = "error",
) -> List[int]:
    """Resolve user-provided allowed site selectors to concrete site_id values.
        If both allowed_sites and allowed_site_specs are provided, their union is used.
        If no constraints are provided, an empty list is returned.

    sites:
        Output of list_reactive_sites().
    allowed_sites:
        Explicit site_id values.
    allowed_site_specs:
        Convenience selectors of the form 'MARK@ORDINAL', e.g. 'C:10@0'.
    invalid_policy:
        One of: 'error', 'warn_skip_seed', 'warn_fallback_all'.

    Returns
    List[int]
        Resolved site_id list
    """

    _MARK_SPEC_RE = re.compile(r"^(?P<mark>[A-Za-z\*]:\d+)@(?P<ord>\d+)$") # mark regex pattern
    n = len(sites)
    if n == 0:
        return []

    valid_ids = {s.site_id for s in sites}

    resolved: List[int] = []

    if allowed_sites:
        for sid in allowed_sites:
            if sid in valid_ids:
                resolved.append(int(sid))
            else:
                _handle_invalid(
                    f"Invalid site_id {sid}. Valid site_id range: 0..{n-1}",
                    invalid_policy,
                    sites,
                )
                if invalid_policy == "warn_skip_seed":
                    return []

    if allowed_site_specs:
        # build lookup from (mark_type, ordinal) -> site_id
        lookup: Dict[Tuple[str, int], int] = {(s.mark_type, s.mark_occurrence_index): s.site_id for s in sites}

        for spec in allowed_site_specs:
            m = _MARK_SPEC_RE.match(spec.strip()) if spec is not None else None
            if not m:
                _handle_invalid(
                    f"Invalid site spec '{spec}'. Expected format 'MARK@ORDINAL', e.g. 'C:10@0'.",
                    invalid_policy,
                    sites,
                )
                if invalid_policy == "warn_skip_seed":
                    return []
                continue

            mark = m.group('mark')
            ord_ = int(m.group('ord'))
            key = (mark, ord_)
            if key in lookup:
                resolved.append(lookup[key])
            else:
                # provide helpful info on available ordinals for this mark
                ords = sorted([s.mark_occurrence_index for s in sites if s.mark_type == mark])
                msg = (
                    f"Invalid site spec '{spec}': mark '{mark}' has available ordinals {ords} "
                    f"in this seed." if ords else f"Invalid site spec '{spec}': mark '{mark}' not found in this seed."
                )
                _handle_invalid(msg, invalid_policy, sites)
                if invalid_policy == "warn_skip_seed":
                    return []

    # deduplicate and keep deterministic order (ascending site_id)
    resolved = sorted(set(resolved))

    if not resolved and (allowed_sites or allowed_site_specs):
        # User requested constraints but none resolved; apply policy.
        if invalid_policy == "warn_fallback_all":
            logger.warning("No valid sites resolved; falling back to all sites.")
            return sorted(valid_ids)

    return resolved

# determine reaction to invalid argument using passed policy
def _handle_invalid(message: str, policy: str, sites: Sequence[ReactiveSite]) -> None:
    """Handle invalid site selection per policy."""
    # raise error
    if policy == "error":
        pretty = [
            {"site_id": s.site_id, "mark_type": s.mark_type, "site_label": s.site_label, "atom_idx": s.atom_idx}
            for s in sites
        ]
        raise ValueError(message + f" Available sites: {pretty}")

    # log issue but try to continue
    if policy in ("warn_skip_seed", "warn_fallback_all"):
        logger.warning(message)
        logger.warning("Available sites: %s", [
            {"site_id": s.site_id, "mark_type": s.mark_type, "site_label": s.site_label, "atom_idx": s.atom_idx}
            for s in sites
        ])
        return

    # Unknown policy -> default to error
    raise ValueError(f"Unknown invalid_policy '{policy}'.")
