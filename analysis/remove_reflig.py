from pdbfixer import PDBFixer
from openmm.app import PDBFile
from argparse import ArgumentParser

fixer = PDBFixer(filename="/home2/esi22219/pdbbind_data/covalent_MAPC/MAPC-0296-1zslx-curated.pdb")
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.removeHeterogens(keepWater=True)
fixer.addMissingHydrogens(pH=7.4)

with open("/home2/esi22219/pdbbind_data/covalent_MAPC/MAPC_protein_no_lig.pdb", "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)
