#!/bin/bash

# Set the directory containing the FASTQ files
FASTQ_DIR="."  # Change this if your files are in a different directory

echo "🔍 Searching for FASTQ files that do NOT end in '_1.fastq' or '_2.fastq'..."

# Find and delete FASTQ files that do not match the expected pattern
find "$FASTQ_DIR" -type f -name "*.fastq" ! -name "*_1.fastq" ! -name "*_2.fastq" -print -delete

echo "✅ Cleanup complete! Only paired-end files (_1.fastq and _2.fastq) remain."
