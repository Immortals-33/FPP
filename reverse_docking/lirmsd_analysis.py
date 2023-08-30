import json
import pandas as pd
import argparse
import numpy as np

"""
Read the ligand-RMSD results from json file and write it do a dataframe or calculate their RMSD of specified quantile.
Should be used right after `ligand_rmsd.py`.
[Required]
"-i", "--input": Input json file containing ligand-RMSD data.
[Optional]
"-o", "--output": The output csv file. Default = "test.csv".
"--order": The way you want the designs to be ordered on the csv file. 
           Choose from "best, worst, mean, pocket". Default = Sample number of designs.
"-t", "--threshold": If this flag is used, then only designs with {mean_rmsd} below this value would be kept.
"-q", "--quantile": If this flag is used, then you can calculate the 25, 50, 75 quantile of overall affinity value.
                    Choose from [25, 50, 75].
"""


parser = argparse.ArgumentParser(description="Generate ligand and pocket RMSD by CSV.")
parser.add_argument("-i", "--input",required=True, help="Input JSON file name")
parser.add_argument("-o", "--output", default="test.csv", help="Output CSV file name")
parser.add_argument("--order", choices=["best", "worst", "mean", "pocket"], default="sample_number",
                    help="Order results by highest, lowest, or mean affinity")
parser.add_argument("-t", "--threshold", type=float, help="Threshold for mean affinity filtering")
parser.add_argument("-q", "--quantile", type=int, choices=[25, 50, 75], help="Quantile value for affinity values")
args = parser.parse_args()

# Load the JSON data
with open(args.input, "r") as json_file:
    data = json.load(json_file)

# Initialize lists to store results
sample_numbers = []
best_RMSD = []
worst_RMSD = []
mean_RMSD = []
pocket_RMSD = []
pose1_RMSD = []

# Iterate through the data and calculate the metrics
for sample_number, sample_data in data.items():
    lirmsd = sample_data["poses"].values()
    best_rmsd = min(lirmsd)
    worst_rmsd = max(lirmsd)
    mean_rmsd = sample_data["mean_Ligand_RMSD"]
    pocket_rmsd = sample_data["Pocket_RMSD"]
    pose1_rmsd = sample_data["poses"]["pose_1"]

    sample_numbers.append(sample_number)
    best_RMSD.append(best_rmsd)
    worst_RMSD.append(worst_rmsd)
    mean_RMSD.append(mean_rmsd)
    pocket_RMSD.append(pocket_rmsd)
    pose1_RMSD.append(pose1_rmsd)


df = pd.DataFrame({
    "Sample Number": sample_numbers,
    "Best ligand-RMSD": best_RMSD,
    "Pose1 RMSD": pose1_RMSD,
    "Wrost ligand-RMSD": worst_RMSD,
    "Mean ligand-RMSD": mean_RMSD,
    "Pocket-RMSD": pocket_RMSD
})

def quantile_calculate(quantile):
    all_mean_RMSD = []
    all_best_RMSD = []
    all_worst_RMSD = []
    all_pocket_RMSD = []

    for sample_data in data.values():
        all_mean_RMSD.append(sample_data["mean_Ligand_RMSD"])
        all_best_RMSD.extend(sample_data["poses"].values()) #Something wrong here
        all_worst_RMSD.extend(sample_data["poses"].values()) #Something wrong here
        all_pocket_RMSD.append(sample_data["Pocket_RMSD"])

    mean_quantile = np.percentile(all_mean_RMSD, quantile)
    best_quantile = np.percentile(all_best_RMSD, quantile)
    worst_quantile = np.percentile(all_worst_RMSD, quantile)
    pocket_quantile = np.percentile(all_pocket_RMSD, quantile)

    print(f"Mean ligand-RMSD of {quantile}th percentile: {mean_quantile:.3f}.")
    print(f"Highest ligand-RMSD of {quantile}th percentile: {best_quantile:.3f}.")
    print(f"Lowest ligand-RMSD of {quantile}th percentile: {worst_quantile:.3f}.")
    print(f"Pocket-RMSD of {quantile}th percentile: {pocket_quantile:.3f}")
    
# Sort the DataFrame based on the specified order
if args.order == "best":
    df = df.sort_values(by="Best ligand-RMSD", ascending=True)
elif args.order == "mean":
    df = df.sort_values(by="Mean ligand-RMSD", ascending=True)
elif args.order == "worst":
    df = df.sort_values(by="Worst ligand-RMSD", ascending=True)
elif args.order == "pocket":
    df = df.sort_values(by="Pocket-RMSD", ascending=True)

if args.threshold is not None:
    df = df[df["Mean ligand-RMSD"] < args.threshold]

if args.quantile is not None:
    q = args.quantile
    quantile_calculate(q)

df.to_csv(args.output, index=False)

print(f"CSV file '{args.output}' created.")
