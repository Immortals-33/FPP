import argparse
import json

from Bio import SeqIO

'''
This script performs two functions on a ProteinMPNN output fasta file:
 1. Rename the fastas based on user-specified format. 
 2. Extract the score and global-score metrics output by ProteinMPNN, which is a negative log-probability value 
    estimated by ProteinMPNN. We consider a lower value corresponds to a slightly more reliable sequence.
    
[Required]
"-i", "--input": ProteinMPNN output fasta.
[Optional]
"-o", "--output": Fasta file after renamed. If not provided, the original fasta file would be covered. 
"-p", "--prefix": User-specified prefix. Each sequences in fasta file would be renamed by {prefix}_{sample_number}.
"-s", "--score-output": Extract the score and global-score from original sequence headers.
"--sort": If provided, then sort the samples based on the ascending order of global score.
'''

def extract_scores(header):
    parts = header.split(", ")
    score = float(parts[2].split("=")[1])
    global_score = float(parts[3].split("=")[1])
    return score, global_score


def parse_args():
    parser = argparse.ArgumentParser(description="Rename sequences in a fasta file.")
    parser.add_argument("-i", "--input", required=True, help="Input fasta file.")
    parser.add_argument("-o", "--output", help="Output fasta file. If not provided, input file will be used.")
    parser.add_argument("-p", "--prefix", help="Prefix for renaming. Default is the original fasta name.")
    parser.add_argument("-s", "--score-output", help="JSON output file for ProteinMPNN scores.")
    parser.add_argument("--sort", help="Sort the samples by global score ascending order.")

    return parser.parse_args()

def main():
    args = parse_args()

    sequence_file = args.input
    output_file = args.output if args.output else sequence_file

    renamed_sequences = []
    sample_number = 1
    score_dict = {}

    with open(sequence_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            header = record.description
            header_parts = record.description.split(", ")
            if not header_parts[0].startswith("T"):
                continue
            
            score, global_score = extract_scores(header)
            
            new_prefix = args.prefix if args.prefix else sequence_file.split("/")[-1].split(".")[0]
            new_header = f"{new_prefix}_{sample_number}"
            renamed_header = f">{new_header}"
            record.id = renamed_header
            record.description = ""
            renamed_sequences.append(record)
            
            score_dict[f"sample_{sample_number}"] = {"score": score, "global score": global_score}
            
            sample_number += 1

    with open(output_file, "w") as file:
        SeqIO.write(renamed_sequences, file, "fasta-2line")

    print("Sequences renamed and saved successfully.")
    
    
    if args.score_output:
        output_file = args.score_output
        if args.sort:
            score_dict = dict(sorted(score_dict.items(), key=lambda x: x[1]["global score"]))
        with open(args.score_output, 'w') as json_file: 
            json.dump(score_dict, json_file, indent = 4)
            
    print(f'Scores saved into {output_file} successfully.')

if __name__ == "__main__":
    main()
