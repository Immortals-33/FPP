# Complexes PDBs' format editing script.
# Written by @Immortals on August 25th, 2023.

"""
When running reverse-docking by Autodock Vina, we use `cat` function from bash to concatenate protein and ligand to a new PDB file.
However, The "END" occurs at the end of protein part, which might cause the ligand not being read by MDAnalysis.
To address this problem, I use this script to move "END" from the protein part to the end of the file.
"""

from pathlib import Path
from typing import *

input_path = Path(
    "/lustre/home/acct-clschf/clschf/zzq/ads_md/docking/reverse_docking/output_r2/"
).resolve()


def format_modified(filename):
    with open(filename, "r") as file:
        lines = file.readlines()
    try:
        lines.remove("END\n")
    except ValueError:
        print(f"Error in file: {filename}")
    lines.append("END\n")

    with open(filename, "w") as file:
        file.writelines(lines)


for subfolder in input_path.glob("*"):
    if not subfolder.is_dir():
        continue

    complex_files = subfolder.glob("fpp_*.pdb")
    complex_paths = [
        complex_path for complex_path in complex_files if complex_path.is_file()
    ]
    for i in complex_paths:
        format_modified(i)

# pdb_filename = "./output_r2/af2_mutateAF_Apo_9991/fpp_af2_mutateAF_Apo_9991_ligand_01.pdb"
# format_modified(pdb_filename)
