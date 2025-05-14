
import random
import argparse
import sys
import edlib
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
# import re
# import os
# import errno
# import shutil
# import glob



# def mkdir_p(path):
#     print("creating", path)
#     try:
#         os.makedirs(path)
#     except OSError as exc:  # Python >2.5
#         if exc.errno == errno.EEXIST and os.path.isdir(path):
#             pass
#         else:
#             raise



# def reverse_complement(string):
#     rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}
#     rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
#     return(rev_comp)

def edit_distance(candidate_frag, frags):
    """
    Check if candidate_frag is identical to any of the fragments in frags
    """
    #for frag in frags:
    #    if edlib.align(candidate_frag, frag)['editDistance'] == 0:
    #        return False
    #return True
    for frag in frags:
        ed1 = edlib.align(candidate_frag, frag, task="distance", mode="HW")['editDistance']
        ed2 = edlib.align(frag, candidate_frag, task="distance", mode="HW")['editDistance']
        if min(ed1, ed2) == 0:
            return False
    return True

def allowed_muts(base):
    lst = ['A','C','G','T']
    lst.remove(base)
    return lst
    
    
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

def main(args):
    frags = [''.join([random.choice('ACGT') for i in range(args.fraglen)])]*20 # 20 copies of the same fragment
    # print(frags)
    for i, frag in enumerate(frags):
        muts = set(random.sample(range(len(frag)),args.mutations)) # sample fraglen*mut_rate mutation sites
        frag_copy = "".join([frag[i] if i not in muts else random.choice(allowed_muts(frag[i])) for i in range(len(frag))])
        frags[i] = frag_copy
    
    # check ed
    eds = []
    for i, frag_i in enumerate(frags): 
        eds_i = []
        for j, frag_j in enumerate(frags): 
            eds_i.append(edlib.align(frag_i, frag_j)['editDistance'])
        eds.append(eds_i)
    for eds_i in eds:
        print(eds_i)
    visualize_eds(eds)

    frags_out = open(args.outfile, 'w')
    for i, frag in enumerate(frags):
        frags_out.write('{0},{1}\n'.format('frag_'+ str(i), frag))
    frags_out.close()




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate references", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outfile', type=str, help='Output file with fragments')
    parser.add_argument('--mutations', type=int, default=4, help='mutations.')
    parser.add_argument('--fraglen', type=int, default=20, help='segment length.')
    parser.add_argument('--nr_frags', type=int, default=20, help='Number of fragments.')

    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)
