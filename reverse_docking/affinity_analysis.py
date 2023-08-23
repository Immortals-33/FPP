import json
import pandas as pd
import argparse
import numpy as np

"""
Read the affinity results from json file and write it do a dataframe or calculate their affinity value of specified quantile.
Should be used right after `get_affinity.py`.
[Required]
"-i", "--input": Input json file containing affinity data.
[Optional]
"-o", "--output": The output csv file. Default = "test.csv".
"--order": The way you want the designs to be ordered on the csv file. 
           Choose from "highest, lowest, mean". Default = Sample number of designs.
"-t", "--threshold": If this flag is used, then only designs with {mean affinity} below this value would be kept.
"-q", "--quantile": If this flag is used, then you can calculate the 25, 50, 75 quantile of overall affinity value.
                    Choose from [25, 50, 75].
"""


parser = argparse.ArgumentParser(description="Generate affinity metrics CSV.")
parser.add_argument("-i", "--input",required=True, help="Input JSON file name")
parser.add_argument("-o", "--output", default="test.csv", help="Output CSV file name")
parser.add_argument("--order", choices=["highest", "lowest", "mean"], default="sample_number",
                    help="Order results by highest, lowest, or mean affinity")
parser.add_argument("-t", "--threshold", type=float, help="Threshold for mean affinity filtering")
parser.add_argument("-q", "--quantile", type=int, choices=[25, 50, 75], help="Quantile value for affinity values")
args = parser.parse_args()

# Load the JSON data
with open(args.input, "r") as json_file:
    data = json.load(json_file)

# Initialize lists to store results
sample_numbers = []
highest_affinities = []
lowest_affinities = []
mean_affinities = []

# Iterate through the data and calculate the metrics
for sample_number, sample_data in data.items():
    affinities = sample_data["poses"].values()
    highest_affinity = max(affinities)
    lowest_affinity = min(affinities)
    mean_affinity = sample_data["mean_affinity"]

    sample_numbers.append(sample_number)
    highest_affinities.append(highest_affinity)
    lowest_affinities.append(lowest_affinity)
    mean_affinities.append(mean_affinity)

# Create a DataFrame
df = pd.DataFrame({
    "Sample Number": sample_numbers,
    "Highest Affinity": highest_affinities,
    "Lowest Affinity": lowest_affinities,
    "Mean Affinity": mean_affinities
})

def quantile_calculate(quantile):
    all_mean_affinities = []
    all_lowest_affinities = []
    all_highest_affinities = []

    for sample_data in data.values():
        all_mean_affinities.append(sample_data["mean_affinity"])
        all_lowest_affinities.extend(sample_data["poses"].values())
        all_highest_affinities.extend(sample_data["poses"].values())

    mean_quantile = np.percentile(all_mean_affinities, quantile)
    lowest_quantile = np.percentile(all_lowest_affinities, quantile)
    highest_quantile = np.percentile(all_highest_affinities, quantile)

    print(f"Mean affinity of {quantile}th percentile: {mean_quantile:.3f}.")
    print(f"Highest affinity of {quantile}th percentile: {highest_quantile:.3f}.")
    print(f"Lowest affinity of {quantile}th percentile: {lowest_quantile:.3f}.")
    
# Sort the DataFrame based on the specified order
if args.order == "highest":
    df = df.sort_values(by="Highest Affinity", ascending=True)
elif args.order == "lowest":
    df = df.sort_values(by="Lowest Affinity", ascending=True)
elif args.order == "mean":
    df = df.sort_values(by="Mean Affinity", ascending=True)

# Filter the DataFrame based on the threshold if provided
if args.threshold is not None:
    df = df[df["Mean Affinity"] < args.threshold]

if args.quantile is not None:
    q = args.quantile
    quantile_calculate(q)

# Write the DataFrame to a CSV file
df.to_csv(args.output, index=False)

print(f"CSV file '{args.output}' created.")
