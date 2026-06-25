from pdbfixer import PDBFixer
from openmm.app import PDBFile
from argparse import ArgumentParser
# simple psbfixer calls to remove any non protein residue
def main(input_filename, output_filename):
    try:
        fixer = PDBFixer(filename=input_filename)              #"/home2/esi22219/pdbbind_data/9s9o/9S9O.pdb")
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.removeHeterogens(keepWater=True)
        fixer.addMissingHydrogens(pH=7.4)
        # need to remove line below to make it more readable and simple
       # hardcoded #with open("/home2/esi22219/pdbbind_data/9s9o/9s9o_non_covalent_no_lig.pdb", "w") as f:
       
        with open(output_filename, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        print(f"Process complete. File saved to {output_filename}")
    except Exception as e:
        print(f"Removal failed\n")
        print(f"Error: {e}")

if __name__ == "__main__":
    parser = ArgumentParser()
    
    parser.add_argument("--input_file",
                        "-i", # shortcut
                        required = True,
                        help = "Input pdb containing co-crystal structure of protein complex and ligand")
    parser.add_argument("--output",
                        "-o", # shortcut
                        required = True,
                        help = "Path to save pdb with ligand removed")
    
    args = parser.parse_args()

    main(args.input_file, args.output)