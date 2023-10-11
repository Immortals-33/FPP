import os
import argparse
import numpy as np
import subprocess
import argparse
import typing as T
from pathlib import Path
from typing import *
from time import time
from multiprocessing import Value
from genericpath import isfile

import MDAnalysis as mda
from MDAnalysis.analysis import rms
from tmtools import tm_align
from tmtools.io import get_residue_data, get_structure

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

'''
A class for calculating structural similarity between proteins (either a set or pairwise)
RMSD: CA-RMSD. Make sure that the number of atoms remains the same.
TM-Score:
    Mode 1: Using python package tmtools to calculate TM-score. Make sure PDB is of valid format.
    Mode 2: Using TMalign (CMD mode) to calculate TM-score. Diffusion PDBs with invalid format are compatible.
            When using Mode 2, make sure your TMalign path is set in advance.
'''

def create_parser():
    parser = argparse.ArgumentParser(
        description="Utils for calculating strucutral similarity"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the input PDB file or folder",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to the reference PDB file",
    )
    parser.add_argument(
        "--mode",
        type=int,
        choices=[1, 2],
        default=1,
        help="Mode for computing TM-Score.",
    )
    return parser


class StructuralSimilarity:
    def __init__(self) -> None:
        pass

    @staticmethod
    def parse_tmalign_output(out):
        rmsd = None
        tm_score = None

        for line in out.split("\n"):
            if "RMSD" in line:
                try:
                    rmsd = float(line.split("=")[1].strip().split(",")[0])
                except ValueError:
                    pass
            elif "TM-score" in line:
                try:
                    tm_score = float(line.split("=")[1].strip().split("(if")[0])
                except ValueError:
                    pass

        return rmsd, tm_score

    def CA_rmsd(self, esm: Path, af2: Path) -> float:
        u1 = mda.Universe(esm)
        u2 = mda.Universe(af2)
        rmsd = rms.rmsd(
            u1.select_atoms("backbone").positions,
            u2.select_atoms("backbone").positions,
            center=True,
            superposition=True,
        )
        return rmsd

    @staticmethod
    def tmtools_calculate(pdb1: Path, pdb2: Path) -> float:
        s1, s2 = (get_structure(pdb1), get_structure(pdb2))
        chain1, chain2 = (next(s1.get_chains()), next(s2.get_chains()))
        c1, seq1 = get_residue_data(chain1)
        c2, seq2 = get_residue_data(chain2)
        res = tm_align(c1, c2, seq1, seq2)
        tm_score = round((res.tm_norm_chain1 + res.tm_norm_chain2) / 2, 3)
        return tm_score

    @staticmethod
    def tmalign_wrapper(template, temp_pdbfile, force_alignment=None):
        try:
            if force_alignment == None:
                p = subprocess.Popen(
                    f'TMalign {template} {temp_pdbfile} | grep -E "RMSD|TM-score=" ',
                    stdout=subprocess.PIPE,
                    shell=True,
                )
            else:
                p = subprocess.Popen(
                    f'TMalign {template} {temp_pdbfile} -I {force_alignment} | grep -E "RMSD|TM-score=" ',
                    stdout=subprocess.PIPE,
                    shell=True,
                )
            output, __ = p.communicate()
            output = output.decode("utf-8")
            tm_rmsd = float(str(output)[:-3].split("RMSD=")[-1].split(",")[0])
            tm_score = float(str(output)[:-3].split("TM-score=")[-1].split("(if")[0])
        except ValueError as e:
            cmd = ["TMalign", template, temp_pdbfile]
            print("1")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            output = result.stdout
            print(f"output: {result.stderr}")
            tm_rmsd, tm_score = StructuralSimilarity.parse_tmalign_output(output)
        return tm_rmsd, tm_score

    def tm_score(
        self, ref: Path, pdb: Path, mode=1
    ) -> Union[float, List[Tuple[str, float]]]:
        PathLike = T.Union[str, Path]
        score = None
        if mode == 1:
            if os.path.isfile(ref) and os.path.isfile(pdb):
                score = StructuralSimilarity.tmtools_calculate(ref, pdb)
            elif os.path.isfile(ref) and os.path.isdir(pdb):
                score = []
                for file in os.listdir(pdb):
                    if file.endswith(".pdb"):
                        pdb_name = file.stem
                        tmscore = StructuralSimilarity.tmtools_calculate(ref, file)
                        score.append((pdb_name, tmscore))
        elif mode == 2:
            if os.path.isfile(ref) and os.path.isfile(pdb):
                score = StructuralSimilarity.tmalign_wrapper(ref, pdb)[1]
            elif os.path.isfile(ref) and os.path.isdir(pdb):
                score = []
                for file in os.listdir(pdb):
                    if file.endswith(".pdb"):
                        ref_path = ref.absolute()
                        pdb_path = file.absolute()
                        pdb_name = pdb_path.stem
                        _, tmscore = tmalign_wrapper(ref_path, pdb_path)
                        score.append((pdb_name, tmscore))
        return score


class Sampler:
    def __init__(self, structural_analyzer: StructuralSimilarity):
        self.structural_analyzer = structural_analyzer

    def get_metrics(self, ref: Path, pdb: Path, metrics: List[str]) -> dict:
        results = {}
        for metric in metrics:
            if metric == "rmsd":
                rmsd = self.structural_analyzer.CA_rmsd(ref, pdb)
                results["RMSD"] = rmsd
            elif metric == "TM-score_mode1":
                tmscore_1 = self.structural_analyzer.tmtools_calculate(ref, pdb)
                results["TM-Score (Mode 1)"] = tmscore_1
            elif metric == "TM-score_mode2":
                tmscore_2 = self.structural_analyzer.tmalign_wrapper(ref, pdb)
                results["TM-Score (Mode 2)"] = tmscore_2
        return results


def main():
    parser = create_parser()
    args = parser.parse_args()

    PathLike = T.Union[str, Path]

    # Create an instance of the StructuralSimilarity class
    structural_analyzer = StructuralSimilarity()

    # Check if the input is a single PDB file or a folder containing PDB files
    if os.path.isfile(args.input):
        # Compute structural similarity between a single PDB file and the reference
        results = structural_analyzer.tm_score(args.reference, args.input, args.mode)
        print("Structural Similarity Score:", results)
    elif os.path.isdir(args.input):
        # Create a Sampler instance to get metrics for multiple files
        sampler = Sampler(structural_analyzer)

        # List of metrics to compute
        metrics = ["TM-score_mode2"]

        # Loop through the PDB files in the input folder
        for file in os.listdir(args.input):
            if file.endswith(".pdb"):
                pdb_path = Path(args.input, file).resolve()
                results = sampler.get_metrics(args.reference, pdb_path, metrics)
                print(f"Metrics for {file}:")
                for metric, value in results.items():
                    print(f"{metric}: {value}")
                print("\n")
    else:
        print(
            "Invalid input. The input should be a single PDB file or a folder containing PDB files."
        )


if __name__ == "__main__":
    main()
