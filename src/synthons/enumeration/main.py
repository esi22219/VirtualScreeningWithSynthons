from .reactions import ReactionIndex
from .synthons import SynthonIndex, extract_marks_from_smiles
from .enumeration_single_step import SeedSpec, SingleStepEnumerator
from .sites import list_reactive_sites

from synthons.SyntOn.SyntOn_BBs import mainSynthonsGenerator

from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

import argparse
from importlib.resources import files
from argparse import ArgumentParser
from pathlib import Path
from datetime import datetime


def sdf_to_smiles(sdf_path, n=None):
    """Reads all molecules in an SDF file and returns a list of their smiles strings

    Args:
        sdf_path (path/str): path to sdf file
        n (int, optional): amount of molecules to read, None for all molecules. Defaults to None.

    Returns:
        list: list of smiles strings
    """
    smiles = []
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)

    for i, mol in enumerate(supplier):
        if mol is None:
            continue  # skip invalid molecules
        smiles.append(Chem.MolToSmiles(mol, canonical=True))

        if n is not None and len(smiles) >= n:
            break

    return smiles

def smi_to_smiles(smi_path, n=None, canonical=True):
    """Reads all molecules/lines in an smi file and returns a list of their smiles strings

    Args:
        sdf_path (path/str): path to sdf file
        n (int, optional): amount of molecules to read, None for all molecules. Defaults to None.

    Returns:
        list: list of smiles strings
    """
    smiles = []
    with open(smi_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # skip blank lines
            if line.startswith("#"):
                continue  # skip comment lines (optional)

            # first token is the SMILES, rest is ignored here
            parts = line.split(maxsplit=1)       
            smi = parts[0]
            name = parts[1].strip() if len(parts) > 1 else None

            if canonical: # canonical makes sure orfering is consistent
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue  # skip invalid smiles
                smi = Chem.MolToSmiles(mol, canonical=True)

            smiles.append((smi, name))

            if n is not None and len(smiles) >= n:
                break

    return smiles


## Pathing, naming, and run id helpers
def file_with_ext(*extensions, must_exist=True):
    """Argparsing type check that ensures an argument is of the correct format

    Args:
        must_exist (bool, optional): check whether the files exists, even if extenstion is allowed. Defaults to True.

    Raises:
        argparse.ArgumentTypeError: Error if type is not allowed or file not found
        

    Returns:
        _type_: _description_
    """
    exts = {e.lower() if e.startswith(".") else "." + e.lower() for e in extensions}

    def _checker(val):
        p = Path(val)

        # check if file exists at all (apart from just extension match)
        if must_exist:
            if not p.exists():
                raise argparse.ArgumentTypeError(f"File does not exist: {p}")
            if not p.is_file():
                raise argparse.ArgumentTypeError(f"Not a file: {p}")
        # then check extension match
        if p.suffix.lower() not in exts:
            raise argparse.ArgumentTypeError(f"Expected one of {sorted(exts)} got '{p.suffix}' for: {p}")
        return p
    return _checker


def make_run_id():
    """Makes a deterministic run id using current datetime

    Returns:
        str: run name using datetime information
    """
    # Safe for directory/file names (no ":" etc.)
    return datetime.now().strftime("%Y_%m_%d_%H%M%S")


def main(seed_path,
         synthon_path,
         config_path,
         running_mode,
         output_dir,
         batch_size,
         run_id,
         rng_seed):
    
    #config_path = "SyntOn/config/Setup.xml"
    # build reaction index
    rxn_index = ReactionIndex.from_setup_xml(config_path)
    
    # dynamically read seed input
    all_synthons = []
    all_names = []
    all_smi = []

    # .sdf
    if seed_path.suffix == ".sdf":
        pass
        smiles = sdf_to_smiles(seed_path)

        for i, (smi, seed_name) in enumerate(smiles, start=1):
            all_smi.append(smi)
            synthons = mainSynthonsGenerator(smi, returnDict=True)
            max_len = -1
            current_name = None
            for key, cls in synthons.items():
                if len(list_reactive_sites(Chem.MolFromSmiles(key))) > max_len:
                    max_len = len(list_reactive_sites(Chem.MolFromSmiles(key)))
                    current_synthons = key
                    current_name = seed_name
           
            all_synthons.append(current_synthons)
            all_names.append(current_name)
       
    # .smi and ensure only one synthon (the all encompassing) 
    else:
        smiles = smi_to_smiles(seed_path)
        for i, (smi, seed_name) in enumerate(smiles, start=1):
            all_smi.append(smi)
            synthons = mainSynthonsGenerator(smi, returnDict=True)
            max_len = -1
            current_name = None
            for key, cls in synthons.items():
                if len(list_reactive_sites(Chem.MolFromSmiles(key))) > max_len:
                    max_len = len(list_reactive_sites(Chem.MolFromSmiles(key)))
                    current_synthons = key
                    current_name = seed_name
           
            all_synthons.append(current_synthons)
            all_names.append(current_name)
            # if len(all_synthons) > 6:
            #     break

    seed_smiles = all_synthons

    # build seedspec instances
    seeds = []
    if not seed_smiles:
        raise Exception("No provided seeds were able to be synthonized and prepped for enumeration. Please check them and try again")
    
    elif len(seed_smiles) == 1:
        
        #print(list_reactive_sites(Chem.MolFromSmiles(seed_smiles[0])))
        spec = SeedSpec(seed_smiles=seed_smiles[0], seed_id="seed_1", allowed_sites=[0])
        seeds.append(spec)
        #print(seed, extract_marks_from_smiles(seed))

    else:
        for i, seed in enumerate(seed_smiles):
            #print(list_reactive_sites(Chem.MolFromSmiles(seed)))
            spec = SeedSpec(seed_smiles=seed, seed_id=f"{all_names[i]}", allowed_sites=[0])
            seeds.append(spec)
            if len(seeds) > 5:
                break
            #print(seed, list_reactive_sites((Chem.MolFromSmiles(seed))))
    
    syn_index = SynthonIndex.from_smi_file(synthon_path)
    
    # sanity check
    print(f"Number of Synthons: {len(syn_index)}")
    print("starting enumeration")

    enum = SingleStepEnumerator(
        synthon_index=syn_index,
        reaction_index=rxn_index,
        marks_compatibility=rxn_index._marks_combinations,
        rng_seed=rng_seed,
        invalid_site_policy="error",  # or "warn_skip_seed"
    )
   
    #  return both named even if one is still empty/None
    gen, summary = enum.enumerate(seeds, batch_size=batch_size, output_mode=running_mode, out_dir=output_dir, run_id=run_id)

    return gen, summary


def cli_entry_point():
    parser = ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--seeds", 
                        type = file_with_ext(".sdf", ".smi"), 
                        required=True, 
                        help = "Path to seeds/scaffolds (sdf or smi file)")
    
    parser.add_argument("--synthons", 
                        type = file_with_ext(".smi"), 
                        required=True, 
                        help = "Path to synthons (smi file)")

    parser.add_argument("--rxn_config", 
                        type = file_with_ext(".xml"),
                        default= files("synthons.SyntOn.config") / "Setup.xml",
                        required=False,
                        help = "Path to reaction setup configuration file")

    parser.add_argument("--mode", 
                        choices = ["stream", "parquet", "both"], 
                        default = "stream", 
                        help = "Choose an output format: " \
                        "streaming generator (stream), parquet files (parquet), or both")
    
    parser.add_argument("-o", "--output_dir",
                        default = None,
                        help = "Where to save parquet outputs (if applicable)")
    
    parser.add_argument("-b", "--batch_size",
                        type = int,
                        default = 500,
                        help = "Batch size for enumeration")
    
    parser.add_argument("--run_name",
                        default = make_run_id(),
                        help = "Run id/name for persistence and separation between runs")
    
    parser.add_argument("--rng_seed",
                        type = int,
                        default = 123,
                        help = "Random seed for when using random site selection during enumeration")
    
    args = parser.parse_args()
   
    if args.mode != "stream" and not args.output_dir:
        parser.error("-o/--output_dir is required for --mode parquet/both")

    if args.output_dir:
        out_path = Path(args.output_dir)

    gen, summary = main(seed_path = args.seeds,
         synthon_path = args.synthons,
         config_path = args.rxn_config,
         running_mode = args.mode,
         output_dir = out_path,
         batch_size = args.batch_size,
         run_id = args.run_name,
         rng_seed = args.rng_seed)
    if gen:
        # do additional analysis on generator here
        for batch in gen: 
            pass


if __name__ == '__main__':
    cli_entry_point()