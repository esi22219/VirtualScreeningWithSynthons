
"""Single-step synthon enumeration.

Implements:
- One synthon addition per product.
- Deterministic site listing based on canonicalized seeds (uses standardization + sites modules).
- Reaction filtering via ReactionIndex.
- Candidate selection via SynthonIndex.
- Output modes: stream | parquet | stream+parquet.

This module is intended as the production enumerator core.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Any
import hashlib
import logging
import random

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from enumeration.standardization import canonicalize_seed_smiles, StandardizedSeed
from enumeration.sites import list_reactive_sites, resolve_allowed_sites, ReactiveSite
from enumeration.reactions import ReactionIndex
from enumeration.synthons import SynthonIndex
from enumeration.output_sinks import ParquetSink, OutputSummary

logger = logging.getLogger(__name__)



def prep_mol_for_enum(molList: Chem.rdchem.Mol): # list of Mol objects
    """
    Necesary preparation of mols for Syn-On enumeration. Checks marks and replaces hydrogen with a dummy leaving group based on marks type. Ensures chemical consistency when editing mols.
    args:
        molList: list of rdkit Mol objects
    reutrns:
        newList: list of modified and prepped rdkit Mol objects
    """
    newList = []
    # labels and atomsForMarking are 1:1 matching, converted to dict for simplicity and faster lookup
    labels_to_marking = {10:23, 20:74, 30:72, 40:104, 50:105, 60:106, 70:107, 21:108, 11:109}
    

    # original lookup lists kept for reference 
    #labels = [10, 20, 30, 40, 50, 60, 70, 21, 11]
    #atomsForMarking = [23, 74, 72, 104, 105, 106, 107, 108, 109]
    atomsForMarkingForDoubleBonds = [72, 104, 105]
    for mol in molList:
        if not mol:
            continue
        mol = Chem.AddHs(mol)
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() != 0: # atomMapNum equivalent to marks
                # replacement dummy atom
                repl = labels_to_marking[atom.GetAtomMapNum()]
                replCount = 0

                # find neighboring hydrogen and replace (do not replace anythinge else)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 1:
                        mol.GetAtomWithIdx(neighbor.GetIdx()).SetAtomicNum(repl)
                        replCount += 1

                        # only replace once if not double bond, else twice
                        if repl not in atomsForMarkingForDoubleBonds and replCount == 1:
                            break
                        elif replCount == 2:
                            break
        mol = Chem.RemoveHs(mol)
        newList.append(mol)
    return newList

@dataclass(frozen=True)
class SeedSpec:
    """
    Dataclass to hoold seed information
    
    Attributes:
        seed_smiles: Canonical Smiles of seed
        seed_id: incremental id for each seed
        allowed_sites: optional user specified reactive sites to enumerate from. use list of ids from list_sites.
        allowed_site_specs; Human readbale string version of allowed sites as alternate form of argument (ex. 'CH10:@0')
    """
    seed_smiles: str
    seed_id: Optional[str] = None
    allowed_sites: Optional[Sequence[int]] = None
    allowed_site_specs: Optional[Sequence[str]] = None


def deterministic_route_id(
    seed_id: str,
    seed_site_id: int,
    synthon_id: int,
    reaction_id: str,
    reactant_order: int,
    product_index: int,
) -> str:
    """Deterministic synthesis route id based on key identifiers and properties."""
    s = f"{seed_id}|{seed_site_id}|{synthon_id}|{reaction_id}|{reactant_order}|{product_index}".encode('utf-8')
    return hashlib.blake2b(s, digest_size=16).hexdigest()


class SingleStepEnumerator:
    """
    Main Enumeration class

    Args:
        synthon_index: SynthonIndex Class object containing synthons collected by marks
        reaction_index: Reaction Index Class object containing reactions collected by marks used
        marks_compatibility: mapping of marks compatible with one another in various reaction schema
        rng_seed: random seed for reproducibility
        invalid_site_policy: how to handle user passing of an invlaid or missing site
    """
    def __init__(
        self,
        synthon_index: SynthonIndex,
        reaction_index: ReactionIndex,
        marks_compatibility: Dict[str, Sequence[str]],
        rng_seed: int = 0,
        invalid_site_policy: str = 'error',
    ):
        self.synthons = synthon_index
        self.reactions = reaction_index
        self.marks_compat = marks_compatibility
        self.rng_seed = int(rng_seed)
        self.invalid_site_policy = invalid_site_policy

    def _choose_site(self, sites: Sequence[ReactiveSite], allowed: Optional[Sequence[int]], rng: random.Random) -> Optional[ReactiveSite]:
        """
        Chooses a reaction site to enumerate from/attach synthons. randomly chosen if not passed by user in allowed
        
        Args:
            sites: list of ReactiveSite class objects for a given seed.
            allowed: optional list of integers denoting user selected sites to enumerate from
            rng: random seed for usage in random site selection
        """

        if not sites: # nothing to choose from
            return None
        if allowed:
            # choose one deterministically given RNG
            sid = rng.choice(list(allowed))
            for s in sites:
                if s.site_id == sid:
                    return s
            print(f"No sits in {allowed} were found, choosing random site\n")
            
        return rng.choice(list(sites))

    def enumerate(
        self,
        seeds: Sequence[SeedSpec],
        batch_size: int = 1000,
        output_mode: str = 'stream',
        out_dir: Optional[str] = None,
        run_id: Optional[str] = None,
    ) -> Tuple[Optional[Iterator[List[Dict[str, Any]]]], Optional[OutputSummary]]:
        """Run enumeration from seeds.
        
        Args:
            seeds: list of seeds to enumerate from
            batch_size: batch size for sharding and multiprocessing
            output_mode: how to handle results: streaming for generator (useful if using results downstream), parquet to save a parquet file (persistence),  stream + parquet for both.
            out_dir: directory for parquet saving if applicable
            run_id: optional name for run as prefix for saved files

        Returns:
            (generator, None) for output_mode='stream'
            (None, summary) for output_mode='parquet'
            (generator, summary) for output_mode='stream+parquet' (summary available after exhaustion)

        Notes:
        - In stream+parquet, the generator yields batches immediately while writing parquet as a side-effect.
        """
        # parse arguments
        output_mode = output_mode.lower().strip()
        if output_mode not in ('stream', 'parquet', 'stream+parquet'):
            raise ValueError("output_mode must be one of 'stream', 'parquet', 'stream+parquet'")

        if output_mode in ('parquet', 'stream+parquet'):
            if not out_dir:
                raise ValueError("out_dir is required for parquet output modes")
            sink = ParquetSink(out_dir=out_dir)
        else:
            sink = None

        rng = random.Random(self.rng_seed)

        summary_holder: Dict[str, OutputSummary] = {}

        def gen() -> Iterator[List[Dict[str, Any]]]:
            """
            Main generation function

            Args:
                batch: list of seeds+synthon batches for parralellization
                n_*: record keeping for logging throughout runs
            """
            batch: List[Dict[str, Any]] = []
            n_skipped_no_sites = 0
            n_seed_parse_fail = 0
            n_records = 0
            n_batches = 0

            
            for seed_idx, spec in enumerate(seeds):
               
                can_res = canonicalize_seed_smiles(spec.seed_smiles, seed_id=spec.seed_id)
                             
                if can_res is None:
                    n_seed_parse_fail += 1
                    continue

                std_seed, can_mol = can_res
                
                # find and resolve site choice for a seed
                sites = list_reactive_sites(can_mol)
                
                if not sites:
                    n_skipped_no_sites += 1
                    continue
                
                allowed = resolve_allowed_sites(
                    sites,
                    allowed_sites=spec.allowed_sites,
                    allowed_site_specs=spec.allowed_site_specs,
                    invalid_policy=self.invalid_site_policy,
                )
                
                site = self._choose_site(sites, allowed if allowed else None, rng)
                
                if site is None:
                    # If allowed_sites resolved empty due to warn_skip_seed, skip.
                    continue

                # match seed marks with compatible partner marks
                seed_mark = site.mark_type
                
                partner_marks = self.marks_compat.get(seed_mark, [])
                if not partner_marks:
                    continue
            
                
                # Seed mol used repeatedly; synthon mol created lazily.
                seed_mol = can_mol
                seed_mol = prep_mol_for_enum([seed_mol])[0]
               
                for syn_mark in partner_marks:
                    # for every compatible marks, get matching potential reactions
                    rxn_list = self.reactions.get_applicable(seed_mark, syn_mark)
                    
                    if not rxn_list:
                        continue
                    # get all synthons with compatible marks
                    syn_ids = self.synthons.ids_with_mark(syn_mark)
                     
                    if not syn_ids:
                        continue
                    
                    # try every compatible synthon
                    for syn_id in syn_ids:
                        syn_rec = self.synthons.get(syn_id)
                        syn_mol = Chem.MolFromSmiles(syn_rec.synthon_smiles)
                        try:
                            syn_mol = prep_mol_for_enum([syn_mol])[0] 
                       
                            if syn_mol is None:
                                continue
                        
                        except IndexError:
                            print(f"issue with synthon: {syn_rec}, skipping")

                            continue
                        # try every compatible reaction
                        for tpl, order in rxn_list:
                            rxn = self.reactions.get_compiled(tpl.reaction_id)
                           
                            reactants = (seed_mol, syn_mol) if order == 0 else (syn_mol, seed_mol)
                            
                            # try to react 
                            try:
                                prod_sets = rxn.RunReactants(reactants)
                                if not prod_sets:
                                    reactants_rev = tuple(reversed(reactants))
                                    prod_sets = rxn.RunReactants(reactants_rev)
                            except Exception:
                                continue
                            
                            # nothing made for this combination, move on
                            if not prod_sets:
                                continue

                            for pidx, prods in enumerate(prod_sets):
                                if not prods:
                                    continue
                                # some rdkit reactions return a tuple of products, here we treat the first as main.
                                # if multiple products are produced, each set corresponds to one route.
                                prod = prods[0]
                                product_smiles = None
                                valid = True
                                failure_reason = None
                                try:
                                    Chem.SanitizeMol(prod)
                                    product_smiles = Chem.MolToSmiles(prod, canonical=True, isomericSmiles=True)
                                    
                                except Exception as e:
                                    valid = False
                                    failure_reason = f"sanitize_or_smiles_prodcut_failed: {e}"
                                    # still attempt to produce a SMILES if possible
                                    try:
                                        product_smiles = Chem.MolToSmiles(prod, canonical=True, isomericSmiles=True)
                                    except Exception:
                                        product_smiles = ''

                                rid = deterministic_route_id(
                                    seed_id=std_seed.seed_id,
                                    seed_site_id=site.site_id,
                                    synthon_id=syn_id,
                                    reaction_id=tpl.reaction_id,
                                    reactant_order=order,
                                    product_index=pidx,
                                )
                                
                                # record
                                rec: Dict[str, Any] = {
                                    'route_id': rid,
                                    'seed_id': std_seed.seed_id,
                                    'seed_canonical_smiles': std_seed.seed_canonical_smiles,
                                    'seed_site_id': site.site_id,
                                    'seed_site_mark_type': site.mark_type,
                                    'seed_site_atom_idx': site.atom_idx,
                                    'synthon_id': syn_id,
                                    'synthon_smiles': syn_rec.synthon_smiles,
                                    'synthon_mark_type': syn_mark,
                                    'reaction_id': tpl.reaction_id,
                                    'reaction_name': tpl.reaction_name,
                                    'reactant_order': order,
                                    'product_smiles': product_smiles,
                                    'product_valid': valid,
                                    'failure_reason': failure_reason,
                                    'run_id': run_id,
                                    'rng_seed': self.rng_seed,
                                }

                                batch.append(rec)
                                n_records += 1

                                if len(batch) >= batch_size:
                                    
                                    # side-effect write
                                    if sink is not None:
                                        sink.consume(batch)
                                    yield batch
                                    n_batches += 1
                                    batch = []

            # final batch
            if batch:
                if sink is not None:
                    sink.consume(batch)
                yield batch
                n_batches += 1

            # finalize summary
            extra = {
                'run_id': run_id,
                'rng_seed': self.rng_seed,
                'n_seed_parse_fail': n_seed_parse_fail,
                'n_skipped_no_sites': n_skipped_no_sites,
            }
            if sink is not None:
                summary_holder['summary'] = sink.finalize(extra=extra, write_manifest=True)
            else:
                summary_holder['summary'] = OutputSummary(out_dir=None, shards=[], n_records=n_records, n_batches=n_batches, extra=extra)

        if output_mode == 'stream':
            return gen(), None

        if output_mode == 'parquet':
            # exhaust generator internally and return summary
            for _ in gen():
                pass
            return None, summary_holder.get('summary')

        # stream+parquet
        return gen(), summary_holder.get('summary')

    def get_last_summary(self, out_dir: str) -> Optional[OutputSummary]:
        """Convenience to read manifest from a parquet out_dir."""
        import json, os
        mp = os.path.join(out_dir, 'manifest.json')
        if not os.path.exists(mp):
            return None
        with open(mp, 'r') as f:
            data = json.load(f)
        return OutputSummary(
            out_dir=data.get('out_dir'),
            shards=data.get('shards', []),
            n_records=data.get('n_records', 0),
            n_batches=data.get('n_batches', 0),
            extra=data.get('extra', {}),
        )
