#!/bin/bash

<< EOF
Extract desired outputs from AlphaFold2 outputs.
AF2 might have multiple outputs, of which we're just interested in some of them.
This script is used to extract these outputs to a new folder and rename them.
[Usage]
sh af2pdbs_extract.sh -s ${source_dir} -t ${target_dir} -p ${file_pattern}
For example, if we want to extract all "ranked_0.pdb" from a dir containing many r2_${i} subdir, then:
sh af2pdbs_extract.sh -s ${source_dir} -t ${target_dir} -p ranked_0.pdb 
EOF  


# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--source)
            SOURCE_DIR="$2"
            shift 2
            ;;
        -t|--target)
            TARGET_DIR="$2"
            shift 2
            ;;
        -p|--pattern)
            FILE_PATTERN="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

if [[ -z "$SOURCE_DIR" || -z "$TARGET_DIR" || -z "$FILE_PATTERN" ]]; then
    echo "Usage: $0 -s source_dir -t target_dir -p file_pattern"
    exit 1
fi

echo "Source Directory: $SOURCE_DIR"
echo "Target Directory: $TARGET_DIR"
echo "File Pattern: $FILE_PATTERN"

# Loop through r2_* folders
for r2_folder in "${SOURCE_DIR}"/r2_*; do
    echo "Processing r2 folder: $r2_folder"
    # Loop through subdirectories in each r2_* folder
    for sub_dir in "${r2_folder}"/*; do
        if [[ -d "$sub_dir" ]]; then
            #echo "Processing subdirectory: $sub_dir"
            #echo "Contents of subdirectory: $(ls "$sub_dir")"
            # Extract files matching the pattern and rename them
            for file_path in ${sub_dir}//"${FILE_PATTERN}"; do
                echo "Filepath is ${file_path}"
                if [[ -f "$file_path" ]]; then
                    # Get the parent folder name and extract the desired part
                    parent_folder=$(basename "$(dirname "$file_path")")
                    new_name=$(echo "$parent_folder" | cut -d'_' -f1-3).pdb
                    echo "Copying: $file_path to ${TARGET_DIR}/${new_name}"
                    cp "$file_path" "${TARGET_DIR}/${new_name}"

                    # Additional debugging output
                    if [ $? -eq 0 ]; then
                        echo "Copy successful!"
                    else
                        echo "Copy failed!"
                    fi
                fi
            done
        fi
    done
done






