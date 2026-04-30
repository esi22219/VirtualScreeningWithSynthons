
"""Synthon records and mark-based indexing for fast lookup.

Phase 2+ scope:
- Build an inverted index: mark_type -> synthon_ids
- Keep minimal metadata in memory (SMILES, marks). RDKit Mol objects are created lazily.

Assumptions
-----------
- Marks in synthon SMILES are always in bracket atom-map form (e.g., [C:10], [n:20]).
- The mark_type is represented as '{symbol}:{mapnum}' where symbol is aromatic lowercase when appropriate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple
import re





def extract_marks_from_smiles(smiles: str) -> Tuple[str, ...]:
    """Extract mark_type tokens like 'C:10' or 'c:21' from a SMILES string."""
    
    # Capture bracket atoms with an atom-map number. We keep the first alphabetic symbol (c/n/o/s etc.).
    _MARK_RE = re.compile(r"\[\s*([A-Za-z]{1}|[cnopsb])[^\]]*?:(\d+)\s*\]")
    
    if not smiles:
        return tuple()
    marks: List[str] = []
    for sym, num in _MARK_RE.findall(smiles):
        # sym is only C, N, O, S, in reaction types
        # so we keep first char of symbol only.
        # Preserve case as found (aromatic typically lowercase in SMILES).
        marks.append(f"{sym}:{num}")
    return tuple(marks)

# synthon class containing incremental id and smiles + marks information
@dataclass(frozen=True)
class SynthonRecord:
    synthon_id: int
    synthon_smiles: str
    marks: Tuple[str, ...]

# synthons indexed by marks type for easy retrieval and lookup later
class SynthonIndex:
    """In-memory synthon store with inverted index by mark_type for efficiency look ups as opposed to testng all."""

    def __init__(self):
        self._records: Dict[int, SynthonRecord] = {}
        self._by_mark: Dict[str, List[int]] = {}

    # add a synthon to index, sort appropriatly
    def add(self, record: SynthonRecord) -> None:
        self._records[record.synthon_id] = record
        for m in set(record.marks):
            self._by_mark.setdefault(m, []).append(record.synthon_id)
    
    # retrieve by id number and not mark
    def get(self, synthon_id: int) -> SynthonRecord:
        return self._records[synthon_id]
        
    # retrieve by mark index (number mark in index)
    def ids_with_mark(self, mark_type: str) -> Sequence[int]:
        return self._by_mark.get(mark_type, [])
    
    # dunder method for length
    def __len__(self) -> int:
        return len(self._records)

    @classmethod
    def from_smiles_iter(cls, smiles_iter: Iterable[str], start_id: int = 0) -> 'SynthonIndex':
        """Build index from an iterable of SMILES strings instead of smi file
        
            Strips whitespace and assumes first entity on line is the smies string
        """
        idx = cls()
        sid = start_id
        for line in smiles_iter:
            if line is None:
                continue
            s = str(line).strip()
            if not s:
                continue
            smi = s.split()[0]

            # dont try to convert to mol for efficiency sake
            marks = extract_marks_from_smiles(smi)
            # could add a filter here for mark type if wanted
            idx.add(SynthonRecord(synthon_id=sid, synthon_smiles=smi, marks=marks))
            sid += 1
        return idx
    
    # feeding file as iterable
    @classmethod
    def from_smi_file(cls, path: str, start_id: int = 0) -> 'SynthonIndex':
        with open(path, 'r') as f:
            return cls.from_smiles_iter(f, start_id=start_id)

    @classmethod
    def from_smiles_list(cls, smiles_list, start_id = 0):
        idx = cls()
        sid = start_id
        for smiles in smiles_list:
            marks = extract_marks_from_smiles(smiles)
            idx.add(SynthonRecord(synthon_id=sid, synthon_smiles=smiles, marks=marks))
            sid+=1
        return idx
