from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple
import logging
import re
import xml.etree.ElementTree as ET

from rdkit.Chem import rdChemReactions as Reactions

logger = logging.getLogger(__name__)



# Only use for xml parsing as it is purely string ops and not chemistry realistic. Breaks rdkit mols
def _normalize_mark_token(token: str) -> Optional[str]:
    """Normalize a mark token found in Labels RHS to a canonical mark_type.

    Examples:
    - 'C:10' -> 'C:10'
    - '11C:10' -> 'C:10'
    - '11CH2:10' -> 'C:10'
    - '15N:20' -> 'N:20'
    - 'nH:20' -> 'n:20'
    - '11c:21' -> 'c:21'

    Returns None if token is empty or does not match expected pattern.
    """
    # simple clean checks
    if token is None:
        return None
    token = token.strip()
    if not token or token.lower() == 'none':
        return None
    
    # matching format using :
    _MARK_TOKEN_RE = re.compile(r"^(?P<lhs>[^:]+):(?P<num>\d+)$")

    m = _MARK_TOKEN_RE.match(token)
    if not m:
        return None

    lhs = m.group('lhs')
    num = m.group('num')

    # Find first alphabetic character as element symbol proxy.
    sym = None
    for ch in lhs:
        if ch.isalpha() or ch == '*':
            sym = ch
            break

    if sym is None or sym == '*':
        return None

    return f"{sym}:{num}"

# determine what is required on mols that are on rhs and lhs of a reaction for later checks
def _parse_labels_to_side_marks(labels: str) -> List[Set[str]]:
    """Parse Labels attribute into a list of allowed mark sets per reactant side.

    The Labels string commonly uses ';' to separate reactant label groups and ',' to separate mappings.

    Returns a list of sets, one per group in order.
    """
    if labels is None:
        return []

    labels = labels.strip()
    if not labels or labels.lower() == 'none':
        return []

    groups = [g.strip() for g in labels.split(';') if g.strip()]
    side_marks: List[Set[str]] = []

    for g in groups:
        marks: Set[str] = set()
        parts = [p.strip() for p in g.split(',') if p.strip()]
        for p in parts:
            if '->' not in p:
                continue
            rhs = p.split('->', 1)[1].strip()
            mt = _normalize_mark_token(rhs)
            if mt:
                marks.add(mt)
        side_marks.append(marks)

    return side_marks

# Holding Reactions for easy retrieval and use later. 
#contains additional info on compatible marks for a reaction
@dataclass(frozen=True)
class ReactionTemplate:
    reaction_id: str
    reaction_name: str
    reconstruction_smarts: str
    labels_raw: str
    allowed_marks_side0: frozenset[str]
    allowed_marks_side1: frozenset[str]

# main index class, similar to that of synthon index, but for reactions instead
class ReactionIndex:
    """Index to retrieve applicable reactions given a (seed_mark, synthon_mark) pair."""

    def __init__(self):
        self._templates: Dict[str, ReactionTemplate] = {}
        self._compiled: Dict[str, object] = {}
        # (markA, markB) -> list of (reaction_id, order)
        # order=0 means reactants=(seed, synthon); order=1 means reactants=(synthon, seed)
        self._pair_map: Dict[Tuple[str, str], List[Tuple[str, int]]] = {}
        
        # hard coded mapping as used in Synt-On
        self._marks_combinations = {'C:10': ['N:20', 'O:20', 'C:20', 'c:20', 'n:20', 'S:20'],
                                    'c:10': ['N:20', 'O:20', 'C:20', 'c:20', 'n:20', 'S:20'],
                                    'c:20': ['N:11', 'C:10', 'c:10'], 
                                    'C:20': ['C:10', 'c:10'],
                                    'c:21': ['N:20', 'O:20', 'n:20'], 
                                    'C:21': ['N:20', 'n:20'],
                                    'N:20': ['C:10', 'c:10', 'C:21', 'c:21', 'S:10'], 
                                    'N:11': ['c:20'],
                                    'n:20': ['C:10', 'c:10', 'C:21', 'c:21'], 
                                    'O:20': ['C:10', 'c:10', 'c:21'],
                                    'S:20': ['C:10', 'c:10'], 
                                    'S:10': ['N:20'], 
                                    'C:30': ['C:40', 'N:40'],
                                    'C:40': ['C:30'], 
                                    'C:50': ['C:50'], 
                                    'C:70': ['C:60', 'c:60'],
                                    'c:60':['C:70'], 
                                    'C:60': ['C:70'], 
                                    'N:40': ['C:30'] }
    @property
    def templates(self) -> Dict[str, ReactionTemplate]:
        return self._templates
    
    # returns reactions that can be applied based on reactant marks
    def get_applicable(self, seed_mark: str, synthon_mark: str) -> List[Tuple[ReactionTemplate, int]]:
        """Return list of (ReactionTemplate, order) for a given mark pair."""
        key = (seed_mark, synthon_mark)
        hits = self._pair_map.get(key, [])
        out: List[Tuple[ReactionTemplate, int]] = []
        for rid, order in hits:
            tpl = self._templates.get(rid)
            if tpl is not None:
                out.append((tpl, order))
        return out
    
    # returns compiled version of reaction ready for application
    def get_compiled(self, reaction_id: str):
        """Get (or cache) compiled RDKit reaction for a template."""
        if reaction_id in self._compiled:
            return self._compiled[reaction_id]
        tpl = self._templates[reaction_id]
        rxn = Reactions.ReactionFromSmarts(tpl.reconstruction_smarts)
        self._compiled[reaction_id] = rxn
        return rxn

    def filter_templates(self, allowed_set):
        """
        Idea is to filter or select specific allowed reactions and ignore the others
        """
        # ooga boooga
        # need to be able to filter by id and other features or attributes as well
        # here is an example of how it can easily be done for selective reactions, very easy
        # probably can even set up in init, i dont think a new function is even necesary but oh well here it is
        # this assumes we pass the allowed at a separate time, and not at initiation but can change that functionality
        if not allowed_set:
            self._templates = self._templates
        else:
            self._templates = {rxn_name: info for rxn_name, info in self._templates.items() if rxn_name in allowed_set}
    # way to build index from file containgin reactions and metadata
    @classmethod
    def from_setup_xml(cls, xml_path: str) -> 'ReactionIndex':
        """Class method to build index from xml file

        Args:
            xml_path (str): path to config file

        Returns:
            ReactionIndex: class instance of reaction index for easy lookup and retrieval during enumeration
        """
        # create instance
        idx = cls()

        # parse tree
        tree = ET.parse(xml_path)
        root = tree.getroot()

        # Find all reaction nodes with SMARTS/Labels/ReconstructionReaction attributes.
        # original structure has reactions nested under <AvailableReactions>.
        for node in root.iter():
            attrib = node.attrib
            if not attrib:
                continue
            if 'ReconstructionReaction' not in attrib and 'SMARTS' not in attrib:
                # look for relevant information and skip if not all there
                continue

            recon = attrib.get('ReconstructionReaction')
            labels = attrib.get('Labels', '')
            name = attrib.get('name', node.tag)

            # Skip placeholder entries if exists
            if recon is None or recon.strip().lower() == 'none':
                continue

            # Derive a stable reaction_id: prefer XML tag (e.g., R1.1) if present.
            # node.tag can be like 'R1.1'. If not, fallback to name.
            reaction_id = node.tag
                        
            # Parse labels
            side_marks = _parse_labels_to_side_marks(labels)
            if len(side_marks) < 2:
                # Some reactions may not be usable for two-reactant reconstruction
                logger.debug("Skipping reaction %s due to insufficient label groups", reaction_id)
                continue

            # check both sides of reaction for important info
            side0 = frozenset(side_marks[0])
            side1 = frozenset(side_marks[1])

            # create template and attempt compilation to catch invalid SMARTS early
            tpl = ReactionTemplate(
                reaction_id=reaction_id,
                reaction_name=name,
                reconstruction_smarts=recon,
                labels_raw=labels,
                allowed_marks_side0=side0,
                allowed_marks_side1=side1,
            )
            
            # ensure valid reaction smarts
            try:
                _ = Reactions.ReactionFromSmarts(recon)
            except Exception as e: # rdkit cant recognize reaction
                logger.warning("Failed to compile ReconstructionReaction for %s (%s): %s", reaction_id, name, e)
                continue
            
            # cache for easy retrieval later
            idx._templates[reaction_id] = tpl
            
            # populate pair mapping with order information
            for a in side0:
                for b in side1:
                    idx._pair_map.setdefault((a, b), []).append((reaction_id, 0))
            for a in side1:
                for b in side0:
                    idx._pair_map.setdefault((a, b), []).append((reaction_id, 1))

        # deduplicate pair map lists while preserving deterministic order
        for k, v in list(idx._pair_map.items()):
            seen = set()
            newv = []
            for rid, order in v:
                if (rid, order) in seen:
                    continue
                seen.add((rid, order))
                newv.append((rid, order))
            idx._pair_map[k] = newv

        logger.info("ReactionIndex loaded %d reaction templates; %d mark-pair entries", len(idx._templates), len(idx._pair_map))
        return idx
