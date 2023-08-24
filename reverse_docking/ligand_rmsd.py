#This is a untested version
#Written by @Immortals on August 24th, 2023
import os
import argparse
import re
import warnings
import logging
import glob
from pathlib import Path
from typing import *

import torch

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

warnings.filterwarnings("ignore")

def create_parser():
    parser = argparse.ArgumentParser(description="Calculating ligand-RMSD based on the superimposition of pockets")
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the input folder",
    )
    parser.add_argument(
        "-o", "--output", type=Path, 
        default= 'ligand_rmsd.txt',
        help="Output file name"
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to the reference complex PDB file",
    )
    parser.add_argument(
        "-t", "--threshold", type=float, help="Ligand-RMSD threshold"
    )
    parser.add_argument(
        "-q", "--quantile", type=float, help = "Quantile for filtering using sequence recovery."
    )
    parser.add_argument(
        "root-path", type=Path, default=Path("./filter"), help="Root path for file operations" 
    )
    parser.add_argument(
        "--motif-d", type=Path, required=True, help="Pocket residues of designed complexes"
    )
    parser.add_argument(
        "--motif-r", type=Path, required=True, help="Pocket residue of reference complex."
    )
    return parser

def natural_sort_key(s: Path) -> List:
    return [
        int(text) if text.isdigit() else text.lower() for text in re.split(r"(\d+)", s.name)
    ]

def pocket_extract(motif:Path) -> str:
    with open(motif, "r") as f:
        residue_num = [line.strip().split('\t')[1] for line in f]
    pocket = "(resid " + " or resid ".join(residue_num) + ")"
    return pocket

def pocket_superimposition(ref: Path, exp:Path, pocket_r:str, pocket_d:str) -> Tuple[float, float]:
    u1, u2 = (mda.Universe(ref), mda.Universe(exp))
    _, rmsd_pocket = align.alignto(u1.select_atoms(f"name CA and {pocket_r}"), 
                                   u2.select_atoms(f"name CA and {pocket_d}"))
    rmsd_ligand = rmsd(u1.select_atoms("name C* or name P* or name O* and chainID d").positions, 
                       u2.select_atoms("name C* or name P* or name O* and chainID d").positions)
    return (rmsd_pocket, rmsd_ligand)

def run(args):
    input_path = args.input.resolve()
    output_file = args.output.resolve()
    reference_file = args.reference.resolve()
    root_path = args.root_path.resolve()
    pocket_design = args.motif_d.resolve()
    pocket_reference = args.motif_r.resolve()
    rmsd_threshold = args.threshold
    basename = args.input.stem

    results = []

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%d/%m/%y %H:%M:%S",
    )
    logger = logging.getLogger(__name__)
    logger.info(f"Calculating ligand-RMSD from {input_path} to {reference_file}......")   

    for subfolder in sorted(input_path.glob("*", key=natural_sort_key)):
        if not subfolder.is_dir():
            continue

        sample_match = re.match(r".*_(\d+)$", subfolder.name)
        if not sample_match:
            continue

        sample_number = int(sample_match.group(1))

        complex_files = sorted((subfolder.glob("fpp_*.pdb")), key=natural_sort_key)
        complex_paths = [
            complex_path for complex_path in complex_files if complex_path.is_file()
        ]

        for i, complex_path in enumerate(complex_paths, start=1):
            pocket_d = pocket_extract(pocket_design)
            pocket_r = pocket_extract(pocket_reference)
            pocket_rmsd, ligand_rmsd = pocket_superimposition(reference_file, complex_path, pocket_r, pocket_d)
            results.append(sample_number, i, ligand_rmsd, pocket_rmsd)

    with output_file.open("w") as out:
        for sample_num, ligand_num, ligand_value, pocket_value in results:
            out.write(f"{sample_num}\t{ligand_num}\t{ligand_value:.3f}\t{pocket_value:.3f}\n")

        if args.threshold and args.quantile:
            parser.error("Cannot use -t and -q flags simultaneously!")

        if args.threshold:
            filtered_nums = [
                sample_num
                for sample_num, ligand_num, ligand_value, pocket_value in results
                if ligand_value <= rmsd_threshold
            ]
            logger.info(f"Number of complexes kept after filtering: {len(filtered_nums)}")

            filter_id_file = root_path / f"filterID_{basename}.txt"
            with filter_id_file.open("w") as f:
                f.write("\n".join(map(str, sorted(filtered_nums))))

        if args.quantile:
            logger.info(f"Calculating ligand-RMSD quantile......")
            rmsd_value = [ligand_value for _, _, ligand_value, _ in results]
            quantile = args.quantile
            q_threshold = torch.quantile(torch.tensor(rmsd_value), quantile).item()
            logger.info(f"Users specified quantile: {quantile}, complexes with highest {100*quantile}% ligand-RMSD would be kept.")

            filtered_nums = [
                sample_num
                for sample_num, ligand_num, ligand_value, pocket_value in results
                if rmsd_value <= q_threshold
            ]
            logger.info(f"Numbers of complexes kept after quantile filtering: {len({filtered_nums})}.")
            filter_id_file = root_path / f"filterID_{basename}_{quantile}.txt"
            with filter_id_file.open("w") as f:
                f.write("\n".join(map(str, sorted(filtered_nums))))

def main():
    parser = create_parser()
    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    main()

                        
