import os
import re
import json
import argparse

"""
Get the results of Autodock Vina and save it as json-dict style.

[Required]
"-i", "--input": The folder containing many .log files. 
[Optional]
"-o", "--output": The output name of json file. Default = "log_file.json".
"""

parser = argparse.ArgumentParser(description="Parse AutoDock Vina log files and create a JSON summary.")
parser.add_argument("-i", "--input", required=True, help="Input directory containing log files")
parser.add_argument("-o", "--output", default="log_data.json", help="Output JSON file name")
args = parser.parse_args()

log_dir = args.input
output_file = args.output

data_dict = {}

for log_file in os.listdir(log_dir):
    if log_file.endswith(".log"):
        sample_number_match = re.search(r"_(\d+)\.log", log_file)
        if sample_number_match:
            sample_number = int(sample_number_match.group(1))
        else:
            continue  # Skip files that don't match the pattern
        
        sample_dict = {"poses": {}, "mean_affinity": 0.0}
        
        with open(os.path.join(log_dir, log_file), "r") as f:
            read_data = False
            affinity_sum = 0.0
            num_poses = 0
            for line in f:
                if "-----+------------+----------+----------" in line:
                    read_data = True
                    continue
                if read_data and line.strip() != "":
                    cols = line.strip().split()
                    if len(cols) >= 3 and cols[0].isdigit():
                        model = int(cols[0])
                        affinity = float(cols[1])
                        pose_key = f"ligand_{model}"
                        sample_dict["poses"][pose_key] = affinity
                        affinity_sum += affinity
                        num_poses += 1
        
        if num_poses > 0:
            sample_dict["mean_affinity"] = round(affinity_sum / num_poses, 3)
        
        data_dict[sample_number] = sample_dict

# Sort the data_dict by sample_number
sorted_data_dict = dict(sorted(data_dict.items()))

# Save the data as JSON
with open(output_file, "w") as f:
    json.dump(sorted_data_dict, f, indent=4)

print("Data extraction and saving complete.")
