import random
import argparse
import sys
import edlib
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

def parse_fragments(file_path):
    fragments = {}

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # skip empty lines
            identifier, sequence = line.split(',')
            prefix = identifier.split('_')[0]

            fragments[identifier]=sequence

    return fragments
    
    
def visualize_eds(eds):
    # Convert to numpy array for plotting
    matrix = np.array(eds)

    # Create heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(matrix, annot=True, fmt="d", cmap="coolwarm", square=True, cbar_kws={"label": "Distance"})
    plt.title("Distance Matrix Heatmap")
    plt.xlabel("Fragment Index")
    plt.ylabel("Fragment Index")
    plt.tight_layout()
    plt.show()   
    
    
     
def visualize_eds2(eds, fragment_names,folder,name):
    plt.figure(figsize=(12, 10))
    sns.heatmap(eds, annot=True, fmt="d", cmap="coolwarm", square=True,
                xticklabels=fragment_names, yticklabels=fragment_names,
                cbar_kws={"label": "Edit Distance"})
    plt.title("Edit Distance Matrix Heatmap")
    plt.xlabel("Fragments")
    plt.ylabel("Fragments")
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    filename = name +".png"
    outfile = os.path.join(folder,filename)
    plt.savefig(outfile, dpi=300)
    plt.close()
    
def main(args):    
    file_path = args.infile
    parsed_fragments = parse_fragments(file_path)
    fragment_names = list(parsed_fragments.keys())
    flat_fragments = list(parsed_fragments.values()) # if parsed_fragments is a dict of single sequences
    # If each entry is a list of dicts (as in previous parser suggestion), flatten it:
    eds = []
    for frag_i in flat_fragments:
        eds_i = []
        for frag_j in flat_fragments:
             #eds_i.append(edlib.align(frag_i, frag_j)['editDistance'])
             ed1 = edlib.align(frag_i, frag_j, task="distance", mode="HW")['editDistance']
             ed2 = edlib.align(frag_j, frag_i, task="distance", mode="HW")['editDistance']
             eds_i.append(min(ed1,ed2))
        eds.append(eds_i)

    # Call visualization
    visualize_eds2(eds, fragment_names,args.folder,args.name)

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate references", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', type=str, help='Input file with fragments')
    parser.add_argument('--folder', type=str, help='output_folder.')
    parser.add_argument('--name', type=str, help='name of the output file.')

    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args) 
        
