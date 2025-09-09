#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# Turn on command echoing for debugging.
set -x

# Script for De Novo Assembly using SPAdes.
# This script processes samples from subdirectories, identifies FASTQ files
# based on sequencing technology, and runs SPAdes for assembly.

# --- Script Configuration and Argument Parsing ---

# Check if the required arguments are provided.
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <path_to_work_dir> <sequencing_tech>"
    echo "  <path_to_work_dir>: The main directory containing sample subdirectories."
    echo "  <sequencing_tech>: Can be 'illumina' or 'nanopore'."
    exit 1
fi

# Get the absolute path of the work directory.
# This makes the script portable regardless of where it's executed.
WORK_DIR=$(realpath "$1")
TECH="$2"
OUTPUT_DIR="$WORK_DIR/De_Novo_Assemblies"
LOG_FILE="$WORK_DIR/De_Novo_Assembly_log.txt"

# Create the main output directory if it doesn't exist.
mkdir -p "$OUTPUT_DIR"
echo "Creating output directory: $OUTPUT_DIR"

# Redirect all stdout and stderr to the log file.
# This ensures that all script output is captured.
exec &> "$LOG_FILE"

echo "--- De Novo Assembly Log Start: $(date) ---"

# --- Main Processing Loop ---

# Navigate to the working directory to ensure relative paths work.
cd "$WORK_DIR" || { echo "Error: Cannot navigate to $WORK_DIR. Exiting."; exit 1; }

# Find all sample subdirectories and process them in a robust loop.
# The 'find ... -print0' and 'while read -d '' ...' pattern handles
# directory names that contain spaces or special characters correctly.
find . -mindepth 1 -maxdepth 1 -type d -print0 | while IFS= read -r -d '' sample_dir; do
    SAMPLE_ID=$(basename "$sample_dir")
    echo "--------------------------------------------------------"
    echo "Processing sample: $SAMPLE_ID"

    # Define a sample-specific output directory for SPAdes.
    SAMPLE_OUTPUT_DIR="$OUTPUT_DIR/$SAMPLE_ID"
    mkdir -p "$SAMPLE_OUTPUT_DIR"

    # --- Illumina Data Assembly ---
    if [ "$TECH" == "illumina" ]; then
        echo "Detected sequencing technology: Illumina"
        # Use 'find' to robustly locate R1 and R2 files.
        # This will handle filenames with spaces correctly.
        FASTQ_R1=$(find "$sample_dir" -maxdepth 1 -name "*_R1_*.fastq.gz" -print -quit)
        FASTQ_R2=$(find "$sample_dir" -maxdepth 1 -name "*_R2_*.fastq.gz" -print -quit)

        # Check if both paired-end files were found.
        if [ -z "$FASTQ_R1" ] || [ -z "$FASTQ_R2" ]; then
            echo "Warning: Paired-end files (*_R1_*.fastq.gz or *_R2_*.fastq.gz) not found for sample '$SAMPLE_ID'. Skipping."
            continue
        fi

        echo "Found R1: $FASTQ_R1"
        echo "Found R2: $FASTQ_R2"
        echo "Running SPAdes for Illumina data..."

        # Run SPAdes with the paired-end flags.
        spades.py --careful --pe1-1 "$FASTQ_R1" --pe1-2 "$FASTQ_R2" -o "$SAMPLE_OUTPUT_DIR"

    # --- Nanopore Data Assembly ---
    elif [ "$TECH" == "nanopore" ]; then
        echo "Detected sequencing technology: Nanopore"
        # Find all Nanopore fastq files and store them in an array.
        # This is more efficient and correct than concatenating them.
        readarray -t FASTQ_FILES < <(find "$sample_dir" -maxdepth 1 \( -name "*.fastq.gz" -o -name "*.fastq" \))

        if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
            echo "Warning: No fastq or fastq.gz files found for sample '$SAMPLE_ID'. Skipping."
            continue
        fi

        echo "Found Nanopore files:"
        printf " - %s\n" "${FASTQ_FILES[@]}"
        echo "Running SPAdes for Nanopore data..."

        # Pass all found files directly to SPAdes using the --nanopore flag.
        spades.py --nanopore --only-assembler --reads "${FASTQ_FILES[@]}" -o "$SAMPLE_OUTPUT_DIR"

    # --- Handle Unrecognized Technology ---
    else
        echo "Error: Sequencing technology '$TECH' not recognized. Please use 'illumina' or 'nanopore'. Skipping."
        continue
    fi

    echo "SPAdes assembly complete for $SAMPLE_ID"
done

echo "--------------------------------------------------------"
echo "De Novo Assembly process finished."
echo "Log file saved to: $LOG_FILE"
echo "--- De Novo Assembly Log End: $(date) ---"
