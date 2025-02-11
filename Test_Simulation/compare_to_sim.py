import csv
import os, sys
import random
import argparse
import errno
import math

#run this file via: python compare_to_sim.py --read_infos /home/alexanderpetri/Project3/SimulationResults/Simulation_Results/readinfos.txt --py_results /home/alexanderpetri/Project3/SimulationResults/Simulation_Results/Kristoffer_out --rs_results /home/alexanderpetri/Project3/SimulationResults/Simulation_Results/rust_out --outfile /home/alexanderpetri/Project3/SimulationResults/Simulation_Results/out.txt


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
                key_part, frag_list = line.split(':', 1)

                # Remove the leading '>' and strip whitespace
                key = key_part[1:].strip().split("|")[-1]

                # Parse the list of items (assumes it's well-formed)
                value = eval(frag_list.strip())  # Using eval to parse the list safely

                # Store in the dictionary
                result[key] = value

        return result


def parse_python_result(file_path):
    """
    Parses a file and returns a dictionary with read accessions as keys and
    fragment IDs as a list sorted by their occurrence in the file.

    Args:
        file_path (str): Path to the input file.

    Returns:
        dict: A dictionary where keys are read accessions and values are lists of fragment IDs.
    """
    result = {}
    with open(file_path, 'r') as file:
        # Read the header line and ignore it
        header = file.readline().strip().split(',')

        # Read the rest of the file
        for line in file:
            line = line.strip()
            if not line:  # Skip empty lines
                continue

            # Split the line into its components
            parts = line.split(',')
            read_acc = parts[0].split("|")[-1]
            fragment_id = parts[1]
            startpos =parts[2]
            endpos = parts[3]
            # Append fragment_id to the corresponding read_acc list
            if read_acc not in result:
                result[read_acc] = []
            result[read_acc].append((fragment_id,startpos,endpos))

    return result
    
    
def get_missing_fragments_ordered(dict1, dict2):
    """
    Finds all fragments that were present in dict1 but not in dict2,
    preserving order and handling cases where dict2 contains (fragment_id, start, stop) tuples.

    Args:
        dict1 (dict): The reference dictionary (baseline), containing lists of fragment IDs.
        dict2 (dict): The dictionary to compare against, containing lists of (fragment_id, start, stop) tuples.

    Returns:
        dict: A dictionary with read accessions as keys and missing fragments (in order) as values.
    """
    # Convert dict2 values from tuples to lists of fragment IDs
    dict2_processed = {key: [frag[0] for frag in value] for key, value in dict2.items()}

    # Now find the shared accessions based on the processed dict2
    shared_accessions = set(dict1.keys()) & set(dict2_processed.keys())

    missing_fragments = {}

    for read_acc in shared_accessions:
        list1 = [frag[0] for frag in dict1[read_acc]]
        list2 = [frag[0] for frag in dict2[read_acc]]  # Copy lists to avoid modifying originals
        print("L1: {} \nL2:{}".format(list1,list2))
        missing = []
        for fragment in list1:
            if fragment in list2:
                list2.remove(fragment)  # Remove first occurrence to respect repetition
            else:
                missing.append(fragment)  # If not in list2, it's missing
        print("Missing {}\n".format(missing))
        if missing:
            missing_fragments[int(read_acc)] = missing  # Store ordered missing elements
    
    return missing_fragments

def parse_rs_results(file_path):
    """
    Parses a file and returns a dictionary where keys are read accessions
    and values are lists of fragment IDs in order of occurrence.

    Args:
        file_path (str): Path to the input file.

    Returns:
        dict: A dictionary with read accessions as keys and lists of fragment IDs as values.
    """
    fragments_dict = {}

    with open(file_path, 'r') as file:
        # Skip the misleading header
        file.readline()

        # Process each line
        for line in file:
            line = line.strip()
            if not line:  # Skip empty lines
                continue

            # Split the line into components
            parts = line.split(',')
            if len(parts) != 5:
                raise ValueError(f"Unexpected number of fields in line: {line}")

            # Clean up and extract fields
            read_acc = parts[0].strip().split("|")[-1]
            fragment_id = parts[1].strip()
            startpos = parts[2]
            endpos = parts[3]
            # Add the fragment_id to the corresponding read accession
            if read_acc not in fragments_dict:
                fragments_dict[read_acc] = []
            fragments_dict[read_acc].append((fragment_id,startpos,endpos))

    return fragments_dict

#TODO: add more analysis measures: how many reads were perfectly reconstructed
#total number of missed fragments globally
#% of missed fragments 5/full number of fragments
#runtime
def compare_fragment_dicts(reference_dict, tested_dict):
    """
    Compares two dictionaries containing read accessions as keys and ordered fragment ID lists as values.
    Outputs meaningful metrics such as the rate of shared information.

    Args:
        dict1 (dict): First dictionary (baseline).
        dict2 (dict): Second dictionary to compare against.

    Returns:
        dict: Comparison results with detailed metrics for each read accession.
    """
    comparison_results = {}
    #shared_accessions = set(dict1.keys()) & set(dict2.keys())
    preserved_cter: int = 0
    for read_acc in reference_dict.keys():
        #print("Tested_dict",tested_dict)
        list_of_ref_infos = reference_dict[read_acc]
        list1 = [elem[0] for elem in list_of_ref_infos]
        list_of_test_infos = tested_dict[read_acc]
        list2 = [elem[0] for elem in list_of_test_infos]
        # Calculate the set of shared fragments
        set1, set2 = set(list1), set(list2)
        shared_fragments = set1 & set2
        total_fragments = len(list1)  # Union of both sets
        #print(shared_fragments)
        #print(total_fragments)
        # Calculate shared fragments' proportion
        shared_rate = len(shared_fragments) / total_fragments if total_fragments > 0 else 0
        #if shared_rate == 1.0:
            #preserved_cter += 1
        # Check if the shared fragments appear in the same order
        shared_in_order = [f for f in list1 if f in shared_fragments]
        shared_in_order_other = [f for f in list2 if f in shared_fragments]
        list_of_undetected=[]
        if shared_in_order == shared_in_order_other:
            #preserved_cter += 1
            order_preserved = True
            shared_rate = 1.0
            preserved_cter += 1
        else:
            order_preserved = False
            list_of_undetected.append(set1-set2)
        # Store results for the current read accession
        comparison_results[read_acc] = {
            'fragments_1': list1,
            'fragments_2': list2,
            'shared_fragments': list(shared_fragments),
            'shared_rate': round(shared_rate, 2),
            'order_preserved': order_preserved,
            'undetected' : list_of_undetected,
        }

    # Calculate overall shared rate across all accessions
    overall_shared_rates=[
        result['shared_rate'] for result in comparison_results.values()
    ]
    overall_rate = sum(overall_shared_rates) / len(overall_shared_rates) if overall_shared_rates else 0

    # Summary
    summary = {
        'overall_shared_rate': round(overall_rate, 2),
        #'read_accessions_compared': len(shared_accessions),
        'fully_preserved_reads':preserved_cter,
        'comparison_details': comparison_results,
    }

    return summary


def main(args):
    sim_info_dict = read_sim_infos(args.read_infos)
    #py_result_dict = parse_python_result(args.py_results)
    rs_result_dict = parse_rs_results(args.rs_results)
    #print(sim_info_dict)
    #print(py_result_dict)
    #print(rs_result_dict)
    #comparison_python={}
    #comparison_python = compare_fragment_dicts(sim_info_dict, py_result_dict)
    comparison_rust = compare_fragment_dicts(sim_info_dict, rs_result_dict)
    missing_frags = get_missing_fragments_ordered(sim_info_dict, rs_result_dict)
    print("Missing:")
    print(missing_frags)
    outfile = open(args.outfile, "w")
    print("Comp_rust {}".format(comparison_rust))
    print("\nComparison with Rust results")
    outfile.write("Rust:\n")
    for read_acc, details in comparison_rust['comparison_details'].items():
        if int(read_acc) in missing_frags.keys():
    	    missing_frag = missing_frags[int(read_acc)] 
        else:
    	    missing_frag={}
        #print(f"\nRead Accession: {read_acc}")
        #print(f"  Shared Fragments: {details['shared_fragments']}")
        #print(f"  Shared Rate: {details['shared_rate']}")
        #print(f"  Order Preserved: {details['order_preserved']}")
        outfile.write(
            "Read Accession: {0}\nShared Fragments: {1}\nOriginal: {2}\nOther: {3}\nShared Rate: {4}\n Order Preserved: {5}\n Undetected {6}\n missing {7}\n".format(read_acc,
                                                                                                         details[
                                                                                                             'shared_fragments'],
                                                                                                         details['fragments_1'],
                                                                                                         details['fragments_2'],
                                                                                                         details[
                                                                                                             'shared_rate'],
                                                                                                         details[
                                                                                                             'order_preserved'],details['undetected'], missing_frag))
    outfile.write("# reads fully preserved: {0} \n Overall Shared Rate: {1}".format(comparison_rust['fully_preserved_reads'],comparison_rust['overall_shared_rate']))
    print("\n # reads fully preserved:", comparison_rust['fully_preserved_reads'])
    print("\nOverall Shared Rate:", comparison_rust['overall_shared_rate'])
    #outfile.write("Overall Shared Rate:{0}\n".format( comparison_python['overall_shared_rate']))
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
    #parser.add_argument('--py_results', type=str, help='PAth to results file of python implementation')
    parser.add_argument('--rs_results', type=str, help='PAth to results file of python implementation.')
    parser.add_argument('--outfile', type = str, help = 'path to and name of output file')
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    main(args)
