#!/bin/bash

INPUT_FILE="samples_run1.txt"  # File containing list of Sample IDs and SRR IDs
OUTPUT_DIR="/gpfs/data/fs71707/mwielsch1/SA_EPI/NCBI_FASTQ/run_1/"  # Directory to save FASTQ files

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

echo "📥 Downloading FASTQ files for runs listed in $INPUT_FILE..."

# Process each line, extracting only the SRR ID (second column)
while IFS=$'\t' read -r sample_id srr_id; do
    # Skip empty lines
    [[ -z "$srr_id" ]] && continue
    
    echo "🚀 Processing: $srr_id (Sample: $sample_id)"

    # Download the SRA file with a 1GB size limit
    prefetch "$srr_id" --max-size 1G  

    # Ensure prefetch was successful
    if [ $? -eq 0 ]; then
        echo "📦 Converting $srr_id to FASTQ..."

        # Convert directly from the SRA cache
        fasterq-dump "$srr_id" --outdir "$OUTPUT_DIR" --split-files --progress

        # Verify FASTQ files were created
        if [ -f "$OUTPUT_DIR/${srr_id}_1.fastq" ] || [ -f "$OUTPUT_DIR/${srr_id}_2.fastq" ]; then
            echo "✅ FASTQ files created for: $srr_id"
        else
            echo "❌ Warning: FASTQ conversion failed for $srr_id."
        fi

        # Cleanup: Remove downloaded SRA files to save space
        echo "🧹 Removing temporary SRA files..."
        rm -f "$srr_id.sra"
    else
        echo "❌ Error: Failed to download $srr_id. Skipping..."
    fi

    echo "✅ Finished processing: $srr_id"
done < <(awk '{print $1}' "$INPUT_FILE")

echo "🎉 All downloads complete! FASTQ files are saved in $OUTPUT_DIR."

