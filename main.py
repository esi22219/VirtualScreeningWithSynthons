from reactions import ReactionIndex
from synthons import SynthonIndex, extract_marks_from_smiles
from enumeration_single_step import SeedSpec, SingleStepEnumerator
from sites import list_reactive_sites

from SyntOn.src.SyntOn_BBs import mainSynthonsGenerator

from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

import argparse
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

            # first token is the SMILES
            smi = line.split()[0]

            if canonical:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue  # skip invalid smiles
                smi = Chem.MolToSmiles(mol, canonical=True)

            smiles.append(smi)

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

        if must_exist:
            if not p.exists():
                raise argparse.ArgumentTypeError(f"File does not exist: {p}")
            if not p.is_file():
                raise argparse.ArgumentTypeError(f"Not a file: {p}")
            
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

    #sdf_file = "/home2/esi22219/pdbbind_data/9s9o/fragments_no_core.sdf"
    
    # dynamically read seed input
    all_synthons = dict()
    all_smi = []
    # .sdf
    if seed_path.suffix == ".sdf":
        smiles = sdf_to_smiles(seed_path)

        for smi in smiles:
            all_smi.append(smi)
            synthons = mainSynthonsGenerator(smi, returnDict=True)
            all_synthons.update(synthons)
       
    # .smi
    else:
        smiles = smi_to_smiles(seed_path)
        for smi in smiles:
            all_smi.append(smi)
            synthons = mainSynthonsGenerator(smi, returnDict=True)
            all_synthons.update(synthons)

    seed_smiles = list(all_synthons.keys())

    # build seedspec handles
    seeds = []
    if not seed_smiles:
        raise Exception("No provided seeds were able to be synthonized and prepped for enumeration. Please check them and try again")
    
    elif len(seed_smiles) == 1:
        
        #print(list_reactive_sites(Chem.MolFromSmiles(seed_smiles[0])))
        spec = SeedSpec(seed_smiles=seed_smiles[0], seed_id=f"seed_1")
        seeds.append(spec)
        #print(seed, extract_marks_from_smiles(seed))

    else:
        for i, seed in enumerate(seed_smiles):
            #print(list_reactive_sites(Chem.MolFromSmiles(seed)))
            spec = SeedSpec(seed_smiles=seed, seed_id=f"seed_{i}")
            seeds.append(spec)
            #print(seed, extract_marks_from_smiles(seed))
    print("processed seeds\n")

    #smi_synthon_path = "/home2/esi22219/notebooks/synton_results/full_run_Synthmode.1.smi"

    #syn_index = SynthonIndex.from_smi_file(smi_synthon_path)
    # synthon index using output from BulkSynthon
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

    gen, summary = enum.enumerate(seeds, batch_size=batch_size, output_mode=running_mode, out_dir=output_dir, run_id=run_id)

    return gen, summary

if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument("--seeds", 
                        type = file_with_ext(".sdf", ".smi"), 
                        required=True, 
                        help = "Path to seeds (sdf or smi file)")
    
    parser.add_argument("--synthons", 
                        type = file_with_ext(".smi"), 
                        required=True, 
                        help = "Path to synthons (smi file)")

    parser.add_argument("--rxn_config", 
                        type = file_with_ext(".xml"),
                        default="SyntOn/config/Setup.xml",
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

    main(seed_path = args.seeds,
         synthon_path = args.synthons,
         config_path = args.rxn_config,
         running_mode = args.mode,
         output_dir = out_path,
         batch_size = args.batch_size,
         run_id = args.run_name,
         rng_seed = args.rng_seed)