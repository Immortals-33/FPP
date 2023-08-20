from Bio.PDB import PDBParser
import typing as T
from pathlib import Path
from typing import List

pdb_parser = PDBParser(QUIET=True)


def calculate_centroid(pdb: Path) -> List:
    structure = pdb_parser.get_structure("protein", pdb)
    key_residue_ids = [
        "A:262",
        "A:271",
        "A:296",
        "A:392",
        "A:396",
        "A:399",
        "A:400",
        "A:439",
        "A:440",
        "A:448",
        "A:449",
        "A:515",
        "A:518",
        "A:519",
        "A:523",
        "A:525",
    ]

    # Extract coordinates of representative atoms
    key_atom_coords = []
    for residue_id in key_residue_ids:
        chain_id, residue_number = residue_id.split(":")
        chain = structure[0][chain_id]
        residue = chain[int(residue_number)]
        # representative_atom = residue["CA"]
        for atom in residue:
            # Use CA atom as a representative
            key_atom_coords.append(atom.get_coord())

    centroid_coords = [sum(coords) / len(coords) for coords in zip(*key_atom_coords)]
    return centroid_coords


print(
    "Centroid Coordinates:",
    calculate_centroid("./af2_predictions/output_r2/af2_mutateAF_Apo_38.pdb"),
)
