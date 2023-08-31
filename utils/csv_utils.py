import pandas as pd
import typing as T
from typing import *
from pathlib import Path

PathLike = T.Union[str, Path]

# Merge different csv files containing different types of data of individual items. Default index = "Sample Number".
def merge_csv(csv1:Path, csv2:Path, output:Path, joint="Sample Number"):
    df1 = pd.read_csv(csv1)
    df2 = pd.read_csv(csv2)
    merged_df = pd.merge(df1, df2, on=f"{joint}")
    merged_df.to_csv(output, index=False)
    return merged_df

# Filter data by specified properties of items. 
# Add some items by indexes of {Sample Number} is also supported.
def csv_filter(csv:Path, output:Path, *args, **kwargs):
    df = pd.read_csv(csv)
    filtered_df = df.query(
        "`Pose1 RMSD` < 3 and \
        `Best ligand-RMSD` < 2.5 or \
        `Best ligand-RMSD` < 2"
    )

    additional_index = [6866, 3363, 3267, 6247, 5838, 6485]
    additional_df = df[df['Sample Number'].isin(additional_index)]
    filtered_df = pd.concat([filtered_df, additional_df], ignore_index=True)
    filtered_df = filtered_df.filter(items=['Sample Number', 'Best ligand-RMSD', 'Pose1 RMSD', 'Mean Affinity', "Sequence Recovery"])

    filtered_df.to_csv(output, index=True)
    return filtered_df

# Write {Sample Number} to a txt file for sequence extractions.
def save_sample_numbers(df: pd.DataFrame, output_txt: Path):
    sample_numbers = df["Sample Number"]
    with open(output_txt, "w") as txt_file:
        for sample_number in sample_numbers:
            txt_file.write(str(sample_number) + "\n")

merge_csv("./pose1_lirmsd.csv", "./affinity_75.csv", "./merged.csv")
merge_csv("./merged.csv", "./r2_true_recovery.csv", "./merged.csv")
filtered_df = csv_filter("./merged.csv", "update_filtered.csv")
save_sample_numbers(filtered_df, "./filterID.txt")

print(filtered_df)
