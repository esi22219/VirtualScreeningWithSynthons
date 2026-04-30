from Bio.PDB import PDBParser, PDBIO, Select

input_pdb = "/home2/esi22219/pdbbind_data/covalent_MAPC/MAPC-0296-1zslx-curated.pdb"
ligand_out = "/home2/esi22219/pdbbind_data/covalent_MAPC/MAPC_ligand.pdb"
protein_out = "/home2/esi22219/pdbbind_data/covalent_MAPC/MAPC_protein_no_lig.pdb"

# Things commonly present in PDBs that are usually NOT the ligand of interest
COMMON_EXCLUDE = {
    "HOH", "WAT", "DOD",   # waters
    "SO4", "PO4", "CL", "BR", "IOD",
    "NA", "K", "MG", "CA", "ZN", "MN", "FE", "CU", "CO", "NI",
    "EDO", "GOL", "PEG", "PG4", "PGE", "MPD", "FMT", "ACT", "ACE",
    "EOH", "IPA", "DMS", "MES", "TRS", "BME"
}

def heavy_atom_count(residue):
    count = 0
    for atom in residue:
        element = (atom.element or "").strip().upper()
        if element != "H":
            count += 1
    return count

def find_ligand_candidates(structure, exclude=COMMON_EXCLUDE, min_heavy_atoms=6):
    candidates = []

    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag, resseq, icode = residue.id
                resname = residue.resname.strip()

                # Standard polymer residue -> skip
                if hetflag == " ":
                    continue

                # Water -> skip
                if hetflag == "W" or resname in {"HOH", "WAT", "DOD"}:
                    continue

                # Exclude common junk
                if resname in exclude:
                    continue

                n_heavy = heavy_atom_count(residue)
                if n_heavy < min_heavy_atoms:
                    continue

                candidates.append({
                    "model_id": model.id,
                    "chain_id": chain.id,
                    "resname": resname,
                    "resseq": resseq,
                    "icode": icode.strip() if isinstance(icode, str) else "",
                    "heavy_atoms": n_heavy,
                    "residue": residue,
                })

    return candidates

def choose_ligand_candidate(candidates):
    if len(candidates) == 0:
        raise ValueError(
            "No ligand candidates found. "
            "Try relaxing the filters (e.g. lower min_heavy_atoms or remove exclusions)."
        )

    print("Ligand candidates found:")
    for i, c in enumerate(candidates, 1):
        print(
            f"{i:2d}. resname={c['resname']:>4s}  chain={c['chain_id']}  "
            f"resseq={c['resseq']}  icode={c['icode'] or '-'}  "
            f"heavy_atoms={c['heavy_atoms']}"
        )

    if len(candidates) == 1:
        print("\nOnly one candidate found -> auto-selecting it.")
        return candidates[0]["residue"]

    # If multiple, choose the largest as a heuristic
    best = max(candidates, key=lambda c: c["heavy_atoms"])
    print(
        f"\nMultiple candidates found -> auto-selecting largest: "
        f"{best['resname']} chain {best['chain_id']} resseq {best['resseq']} "
        f"(heavy_atoms={best['heavy_atoms']})"
    )
    return best["residue"]

class LigandSelect(Select):
    def __init__(self, target_residue):
        self.target_residue = target_residue

    def accept_residue(self, residue):
        return residue is self.target_residue

class ProteinNoLigSelect(Select):
    def __init__(self, target_residue, keep_water=True):
        self.target_residue = target_residue
        self.keep_water = keep_water

    def accept_residue(self, residue):
        # Remove the chosen ligand
        if residue is self.target_residue:
            return False

        hetflag = residue.id[0]
        resname = residue.resname.strip()

        # Keep standard polymer residues
        if hetflag == " ":
            return True

        # Keep or drop water
        if hetflag == "W" or resname in {"HOH", "WAT", "DOD"}:
            return self.keep_water

        # Drop all other heterogens from the protein file
        return False

parser = PDBParser(QUIET=True)
structure = parser.get_structure("pdb", input_pdb)

candidates = find_ligand_candidates(structure, min_heavy_atoms=6)
target_ligand = choose_ligand_candidate(candidates)

io = PDBIO()
io.set_structure(structure)
io.save(ligand_out, LigandSelect(target_ligand))
io.save(protein_out, ProteinNoLigSelect(target_ligand, keep_water=True))

print(f"\nSaved ligand to: {ligand_out}")
print(f"Saved protein without ligand to: {protein_out}")
