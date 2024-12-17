import csv
import os, sys
import random
import argparse
import errno
import math



def read_sim_infos(filename):
        """
        Parses a file with the given >sim|sim|<number>: ['item1', 'item2', ...] format.

        Args:
            file_path (str): Path to the input file.

        Returns:
            dict: A dictionary where keys are 'sim|sim|<number>' and values are lists of strings.
        """
        result = {}
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if not line or not line.startswith('>'):  # Skip empty lines or invalid lines
                    continue

                # Split key and value part at the first ':'
                key_part, value_part = line.split(':', 1)

                # Remove the leading '>' and strip whitespace
                key = key_part[1:].strip()

                # Parse the list of items (assumes it's well-formed)
                value = eval(value_part.strip())  # Using eval to parse the list safely

                # Store in the dictionary
                result[key] = value

        return result



def main(args):
    sim_info_dict = read_sim_infos(args.read_infos)
    print(sim_info_dict)
    
    #print("Generating "+str(args.nr_reads)+" different reads")
    #isoforms=generate_isoforms(args, genome_out.name)
    #reads = generate_reads(args)
    print("Hello World")
    #simulate_reads(args, reads)
            #print("Simulating reads")
    #sys.exit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate reads for our experiments")
    parser.add_argument('--read_infos', type=str,
                        help='Path to read infos file from simulation script')
    parser.add_argument('--py_results', type=str, help='PAth to results file of python implementation')
    parser.add_argument('--rs_results', type=str, help='PAth to results file of python implementation.')
    parser.add_argument('--outfile', type = str, help = 'path to and name of output file')
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    main(args)
