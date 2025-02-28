import argparse
import random

def generate_random_dna(length: int) -> str:
    """Generates a random DNA sequence of a given length."""
    return ''.join(random.choices("ACGT", k=length))

def generate_fastq(num_reads: int, read_length: int, output_file: str):
    """Generates a FASTQ file with random DNA sequences."""
    with open(output_file, "w") as f:
        for i in range(1, num_reads + 1):
            sequence = generate_random_dna(read_length)
            quality = "I" * read_length  # Phred quality score placeholder ('I' represents a high-quality score)
            f.write(f"@sim_rand_{i}\n")
            f.write(f"{sequence}\n")
            f.write(f"+\n")
            f.write(f"{quality}\n")

def main():
    parser = argparse.ArgumentParser(description="Generate a FASTQ file with random DNA sequences.")
    parser.add_argument("num_reads", type=int, help="Number of reads to generate.")
    parser.add_argument("read_length", type=int, help="Length of each read.")
    parser.add_argument("--output", type=str, default="random_sequences.fastq", help="Output FASTQ file name.")
    
    args = parser.parse_args()
    generate_fastq(args.num_reads, args.read_length, args.output)
    print(f"FASTQ file '{args.output}' generated with {args.num_reads} reads of length {args.read_length}.")

if __name__ == "__main__":
    main()

