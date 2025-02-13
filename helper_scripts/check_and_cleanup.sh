#!/bin/bash

# Set the run name dynamically
RUN_NAME="run_1"  # Change this for different runs
INPUT_FILE="samples_${RUN_NAME}.txt"  # The input file with sample IDs

# Directories
FASTQ_DIR="/gpfs/data/fs71707/mwielsch1/SA_EPI/NCBI_FASTQ/${RUN_NAME}/"  # Where FASTQ files are stored
RESULT_DIR="/gpfs/data/fs71707/mwielsch1/SA_EPI/RESULTS/${RUN_NAME}/"  # Where results will be saved

FAILED_OUTPUT="${RESULT_DIR}${RUN_NAME}_download_failed.txt"  # Output file for failed downloads

# Ensure the result directory exists
mkdir -p "$RESULT_DIR"

# Create or overwrite the failed output file with the header
echo -e "ID\treason_for_fail" > "$FAILED_OUTPUT"

echo "🔍 Checking for failed downloads and conversions..."
while IFS=$'\t' read -r sample_id srr_id; do
    # Skip empty lines
    if [[ -z "$sample_id" || -z "$srr_id" ]]; then
        continue  
    fi

    echo "🧐 Checking: $sample_id ($srr_id)"

    # Initialize status flags
    folder_exists=false
    fastq_exists=false
    paired_end=false
    single_end=false

    # Define expected FASTQ file paths
    FASTQ_1="$FASTQ_DIR/${sample_id}_1.fastq"
    FASTQ_2="$FASTQ_DIR/${sample_id}_2.fastq"
    FASTQ_SINGLE="$FASTQ_DIR/${sample_id}.fastq"

    # Check if a folder exists for SRR ID or Sample ID (SAMN)
    if [ -d "$FASTQ_DIR/$srr_id" ] || [ -d "$FASTQ_DIR/$sample_id" ]; then
        folder_exists=true
    fi

    # Check if FASTQ files exist
    if [ -f "$FASTQ_1" ] || [ -f "$FASTQ_2" ] || [ -f "$FASTQ_SINGLE" ]; then
        fastq_exists=true
    fi

    # Check for paired-end files
    if [ -f "$FASTQ_1" ] && [ -f "$FASTQ_2" ]; then
        paired_end=true
    fi

    # Check if a single-end FASTQ file exists
    if [ -f "$FASTQ_SINGLE" ]; then
        single_end=true
    fi

    # Determine failure reason
    if ! $folder_exists && ! $fastq_exists; then
        echo -e "$sample_id\tdownload_failed" >> "$FAILED_OUTPUT"
        echo "    ❌ Marked as: download_failed"
    elif $folder_exists && ! $fastq_exists; then
        echo -e "$sample_id\tfastq_conversion_failed" >> "$FAILED_OUTPUT"
        echo "    ❌ Marked as: fastq_conversion_failed (Folder exists but no FASTQ files)"
    elif $fastq_exists && ! $paired_end; then
        echo -e "$sample_id\tno_paired_end_data" >> "$FAILED_OUTPUT"
        echo "    ❌ Marked as: no_paired_end_data"
    else
        echo "    ✅ No issues detected"
    fi
done < "$INPUT_FILE"

echo "✅ Finished checking! Results saved to: $FAILED_OUTPUT"

# ---- CLEANUP ----
echo "🧹 Cleaning up unpaired FASTQ files..."
find "$FASTQ_DIR" -type f -name "*.fastq" ! -name "*_1.fastq" ! -name "*_2.fastq" -print -delete

echo "🗂 Removing leftover folders..."
find "$FASTQ_DIR" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +

echo "🎉 Cleanup complete!"
