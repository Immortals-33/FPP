# Ligand-RMSD Calculate script.
# Written by @Immortals on August 24th, 2023
import argparse
import re
import warnings
import logging
import glob
import json
import csv
from pathlib import Path
from typing import *

import torch

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

warnings.filterwarnings("ignore")

"""
This script is for ligand (substrate/molecule)-RMSD calculating.
Steps are included:
1. Specify the pocket(motif) residues.
2. Superimpose the pocket(motif), given the pocket(motif)-RMSD is smallest.
3. Calculate the ligand(substrate)-RMSD between reference structure and different poses. 
4. Sort and filter using either a given threshold or quantile of ligand(substrate)-RMSD distribution.

You need to prepare the following steps for running this script:
1. A motif file containing pocket(motif) information. The format should be carefully checked.
2. A reference complex structure.
3. A set of complex structure candidates.

The script will output following things:
1. A result file containing information of ligand(substrate)-RMSD and pocket(motif)-RMSD.
   Currently, the support format of outputfile includes txt, csv and json. 
2. A filterID file containing samples' ID of which successfully pass the user-specified threshold.

Usage:
[Required]
"-i", "--input": The folder containing complex PDBs candidates.
"-o", "--output": The name of result file. Currently only txt and csv file are supported.
"-r", "--reference": Path to the reference complex file.
"--motif_d": Path to the txt containing pocket(motif) information of designed complexes PDBs.
"--motif_r": Path to the txt containing pocket(motif) information of reference complex PDB.
[Optional]
"-o", "--output": The name of result file, using this flag is strongly recommended. Default = "recovery.txt"
"-t", "--threshold": Sequence recovery threshold. Sequences above this value would be kept, otherwise filtered out.
"--root-path": Root path of operation. Default = My own motif-RMSD filtering folder.
"-q", "--quantile": Quantile of candidates you want to keep. Sequences with the highest ${quantile} sequence recovery would be kept.
                    For example, if you input "-q 75", then the script will calculate the 75% quantile of overall recovery,
                    sequences with recovery higher than this value would be kept, otherwise filtered out.
                    Suggested options: [25, 50, 75]
"""


def create_parser():
    parser = argparse.ArgumentParser(
        description="Calculating ligand-RMSD based on the superimposition of pockets"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the input folder",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default="ligand_rmsd.txt", help="Output file name"
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to the reference complex PDB file",
    )
    parser.add_argument("-t", "--threshold", type=float, help="Ligand-RMSD threshold")
    parser.add_argument(
        "-q",
        "--quantile",
        type=float,
        help="Quantile for filtering using sequence recovery.",
    )
    parser.add_argument(
        "--root-path",
        type=Path,
        default=Path(
            "/lustre/home/acct-clschf/clschf/zzq/ads_md/docking/reverse_docking/filter/"
        ),
        help="Root path for file operations",
    )
    parser.add_argument(
        "--motif-d",
        type=Path,
        required=True,
        help="Pocket residues of designed complexes",
    )
    parser.add_argument(
        "--motif-r",
        type=Path,
        required=True,
        help="Pocket residue of reference complex.",
    )
    return parser


def natural_sort_key(s: Path) -> List:
    return [
        int(text) if text.isdigit() else text.lower()
        for text in re.split(r"(\d+)", s.name)
    ]


def pocket_extract(motif: Path) -> str:
    with open(motif, "r") as f:
        residue_num = [line.strip().split("\t")[1] for line in f]
    pocket = "(resid " + " or resid ".join(residue_num) + ")"
    return pocket


def pocket_superimposition(
    ref: Path, exp: Path, pocket_r: str, pocket_d: str
) -> Tuple[float, float]:
    u1 = mda.Universe(ref)
    u2 = mda.Universe(exp)
    _, rmsd_pocket = align.alignto(
        u1.select_atoms(f"name CA and {pocket_r}"),
        u2.select_atoms(f"name CA and {pocket_d}"),
    )
    rmsd_ligand = rmsd(
        u1.select_atoms("name C* or name P* or name O* and chainID d").positions,
        u2.select_atoms("name C* or name P* or name O* and chainID d").positions,
    )
    return round(rmsd_pocket, 3), round(rmsd_ligand, 3)


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

    subfolders = sorted(input_path.glob("*"), key=natural_sort_key)
    for subfolder in subfolders:
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
            pocket_rmsd, ligand_rmsd = pocket_superimposition(
                reference_file, complex_path, pocket_r, pocket_d
            )
            logger.info(
                f"Sample: {Path(complex_path).name}, ligand-RMSD: {ligand_rmsd}, pocket-RMSD: {pocket_rmsd}."
            )
            results.append((sample_number, i, ligand_rmsd, pocket_rmsd))

    with output_file.open("w", newline="") as out:
        if output_file.suffix == ".csv":
            csv_writer = csv.writer(out, delimiter=",")
            csv_writer.writerow(["Sample", "Pose", "Ligand_RMSD", "Pocket_RMSD"])
            for sample_num, ligand_num, ligand_value, pocket_value in results:
                csv_writer.writerow(
                    [sample_num, ligand_num, ligand_value, pocket_value]
                )
        elif output_file.suffix == ".txt":
            out.write("Sample\tPose\tLigand_RMSD\tPocket_RMSD\n")
            for sample_num, ligand_num, ligand_value, pocket_value in results:
                out.write(
                    f"{sample_num}\t{ligand_num}\t{ligand_value}\t{pocket_value}\n"
                )
        elif output_file.suffix == ".json":
            data_dict = {}
            for sample_num, ligand_num, ligand_value, pocket_value in results:
                pose_key = f"pose_{ligand_num}"
                if sample_num not in data_dict:
                    data_dict[sample_num] = {
                        "poses": {},
                        "mean_Ligand_RMSD": 0,
                        "Pocket_RMSD": 0,
                    }
                data_dict[sample_num]["poses"][pose_key] = ligand_value
                data_dict[sample_num]["Pocket_RMSD"] = pocket_value

            for sample_num in data_dict:
                sample_dict = data_dict[sample_num]
                poses = sample_dict["poses"]
                num_poses = len(poses)
                if num_poses > 0:
                    ligand_rmsd_sum = sum(poses[pose_key] for pose_key in poses)
                    sample_dict["mean_Ligand_RMSD"] = round(
                        ligand_rmsd_sum / num_poses, 3
                    )

            with output_file.open("w") as f:
                json.dump(data_dict, f, indent=4)
        else:
            raise ValueError("The format of output file must be CSV, TXT or JSON!")

        if args.threshold and args.quantile:
            parser.error("Cannot use -t and -q flags simultaneously!")

        if args.threshold:
            filtered_nums = [
                sample_num
                for sample_num, ligand_num, ligand_value, pocket_value in results
                if ligand_value <= rmsd_threshold
            ]
            logger.info(
                f"Number of complexes kept after filtering: {len(filtered_nums)}"
            )

            filter_id_file = root_path / f"filterID_{basename}.txt"
            with filter_id_file.open("w") as f:
                f.write("\n".join(map(str, sorted(filtered_nums))))

        if args.quantile:
            logger.info(f"Calculating ligand-RMSD quantile......")
            rmsd_value = [ligand_value for _, _, ligand_value, _ in results]
            quantile = args.quantile
            q_threshold = torch.quantile(torch.tensor(rmsd_value), quantile).item()
            logger.info(
                f"Users specified quantile: {quantile}, complexes with highest {100*quantile}% ligand-RMSD would be kept."
            )

            filtered_nums = [
                sample_num
                for sample_num, ligand_num, ligand_value, pocket_value in results
                if rmsd_value <= q_threshold
            ]
            logger.info(
                f"Numbers of complexes kept after quantile filtering: {len({filtered_nums})}."
            )
            filter_id_file = root_path / f"filterID_{basename}_{quantile}.txt"
            with filter_id_file.open("w") as f:
                f.write("\n".join(map(str, sorted(filtered_nums))))


def main():
    parser = create_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
