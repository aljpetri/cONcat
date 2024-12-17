#!/bin/bash

# written by Mai Nguyen: This is a bash script to run the Python program for DNA fragment covering with multi-threads
# altered by Alexander Petri to be usable on the Rust implementation
# Check if the correct number of arguments is passed
if [ "$#" -lt 4 ] || [ "$#" -gt 6 ]; then
    exit 1
fi

# Assign command-line arguments to variables
FASTQ_DIR="$1"
RUST_BINARY="$2" # the absolute path to the rust binary
EXPECTED_CSV="$3"
OUTPUT_DIR="$4"
THREADS="${5:-10}"  # Default threads to use: 10
ERROR_FILE="${6:-error_reads.csv}"  # Default to 'error_reads.csv' if not provided

# Create the output directory if it does not exist
mkdir -p "$OUTPUT_DIR"

# Create the subdirectories
mkdir -p "$OUTPUT_DIR/per_read"
mkdir -p "$OUTPUT_DIR/per_alignment"
mkdir -p "$OUTPUT_DIR/alignments"

# Check if the error file exists and delete it
if [ -f "$ERROR_FILE" ]; then
    rm "$ERROR_FILE"
fi

# Use find to locate all .fastq files and pipe them to xargs for parallel processing
find "$FASTQ_DIR" -name "*.fastq" | xargs -I {} -P "$THREADS" bash -c '
    fastq_file="$1"
    rust_binary="$2"
    expected_csv="$3"
    output_dir="$4"
    error_file="$5"

    # Extract base name from the fastq file
    base_name=$(basename "$fastq_file" .fastq)
    
    # Define output paths explicitly using the full path of $output_dir
    #output_aln_csv="${output_dir}/per_alignment/${base_name}_aln.csv"
    output_read_csv="${output_dir}/per_read/${base_name}_read.csv"
    #output_aln_file="${output_dir}/alignments/${base_name}.aln"

    # Run the Python script for each fastq file, specifying both output CSVs and the error file
    /usr/bin/time -v "$rust_binary" --expected "$expected_csv" --fastq "$fastq_file" --outfile "$output_read_csv"  
' _ {} "$RUST_BINARY" "$EXPECTED_CSV" "$OUTPUT_DIR" "$ERROR_FILE"
