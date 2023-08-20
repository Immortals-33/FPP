#Extract PDBs using sample numbers from indexes.
#Note that this is different from the `extract_pdbs.py` in alphafold folder, for which is used to extract predictions from AlphaFold2 results.
import argparse
import re
import os
import shutil


def extract_pdbs(input_folder, index_file, output_folder=None):
    # Read the index numbers from the index file
    with open(index_file, "r") as index_f:
        sample_numbers = set(map(int, index_f.read().split()))

    # Extract PDB files with matching sample numbers
    if output_folder is None:
        output_folder = "./extracted_pdbs"
    os.makedirs(output_folder, exist_ok=True)

    for pdb_file in os.listdir(input_folder):
        if pdb_file.endswith(".pdb"):
            sample_number_match = re.search(r"(\d+)\.pdb$", pdb_file)
            if sample_number_match:
                sample_number = int(sample_number_match.group(1))
                if sample_number in sample_numbers:
                    pdb_path = os.path.join(input_folder, pdb_file)
                    output_pdb_path = os.path.join(output_folder, pdb_file)
                    shutil.copy(pdb_path, output_pdb_path)
                    print(
                        f"Extracted PDB {pdb_file} with sample number {sample_number}."
                    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract PDB files based on sample numbers in an index file"
    )
    parser.add_argument(
        "-i",
        "--input_folder",
        type=str,
        required=True,
        help="Input folder containing PDB files",
    )
    parser.add_argument(
        "-n", "--index", type=str, required=True, help="Path to the index txt file"
    )
    parser.add_argument(
        "-o", "--output_folder", type=str, help="Output folder for extracted PDB files"
    )
    args = parser.parse_args()

    input_folder = args.input_folder
    index_file = args.index
    output_folder = args.output_folder

    extract_pdbs(input_folder, index_file, output_folder)
