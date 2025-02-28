import csv
import os, sys
import random
import argparse
import errno
import math


#adds errors to the reads

def simulate_errors(sequence: str, error_profile: list) -> tuple:
    """Simulates errors in a given DNA sequence based on an error profile and generates quality values."""
    mutated_seq = []
    qual_scores = []
    seq_length = len(sequence)
    avg_error_rate = sum(error_profile[-50:]) / 50  # Average error rate for last 50 bases
    
    for i, base in enumerate(sequence):
        error_rate = error_profile[i] if i < 100 else avg_error_rate
        phred_score = max(1, round(-10 * math.log10(error_rate)))  # Ensure valid Phred score
        
        if random.uniform(0, 1) < error_rate:
            r = random.uniform(0, 1)
            if r < 0.4:  # Deletion
                continue
            elif 0.4 <= r < 0.7:  # Substitution
                mutated_seq.append(random.choice("ACGT".replace(base, '')))
                qual_scores.append(phred_score)
            else:  # Insertion
                mutated_seq.append(base)
                qual_scores.append(phred_score)
                mutated_seq.append(random.choice("ACGT"))
                qual_scores.append(phred_score)
        else:
            mutated_seq.append(base)
            qual_scores.append(phred_score)
    
    quality_string = "".join(chr(q + 33) for q in qual_scores)  # Convert to Phred ASCII
    return "".join(mutated_seq), quality_string
def generate_fastq_from_reads(reads: list, output_file: str, error_profile: list):
    """Generates a FASTQ file from given reads, adding errors based on an error profile."""
    with open(output_file, "w") as f:
        for i, sequence in enumerate(reads, start=1):
            mutated_sequence, quality = simulate_errors(sequence, error_profile)
            f.write(f"@sim_rand_{i}\n")
            f.write(f"{mutated_sequence}\n")
            f.write(f"+\n")
            f.write(f"{quality}\n")

def simulate_reads(args, isoforms):
    outfile = open(os.path.join(args.outfolder,"simulated.fastq"), "w")
    is_fastq = True

    reads = {}
    #different error rates we can add to our reads
    #error_lvls =[1.0]  # 0% error rate
    #error_lvls = [0.99] # 1% error rate
    #error_lvls = [0.95] # 5% error rate
    #error_lvls = [0.9] # 10% error rate
    error_lvls = [0.85] # 15% error rate
    #error_lvls = [0.8] # 20% error rate
    #error_lvls = [0.75] # 25% error rate
    #error_lvls= [0.7]  # 30% error rate
    #error_lvls = [0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96, 0.98, 0.99, 0.995]
    #error_lvls=[0.9,0.95,0.96,0.98,0.99,0.995]#3.94%error rate
    #error_lvls=[0.8, 0.875,0.9,0.92,0.96,0.98,0.99,0.995]#7%error rate
    
    
    
    for i_acc, isoform in isoforms.items():
        read = []
        qual = []
        was_del = False

        for l, n in enumerate(isoform):
            p_correct_reading = random.choice(error_lvls)
            p_error = 1.0 - p_correct_reading
            r = random.uniform(0, 1)
            if r > p_correct_reading:
                error = True
            else:
                error = False

            if error:
                r = random.uniform(0, 1)
                if r < 0.4:  # deletion(those values depend on current base caller )
                    was_del = p_error
                    pass
                elif 0.4 <= r < 0.7: #substitution
                    #we do not want to substitute the same base, so we drop the current base from sub_possibilities
                    sub_possibilities="ACGT".replace(n,'')
                    read.append(random.choice(sub_possibilities))
                    if p_error > 0:
                        qual.append(round(-math.log(p_error, 10) * 10))
                    else:
                        qual.append(1) 

                else: #insertion
                    read.append(n)
                    if p_error > 0:
                        qual.append(round(-math.log(p_error, 10) * 10))
                    else:
                        qual.append(1) 

                    r_ins = random.uniform(0, 1)
                    ins_len=1
                    while r_ins >= 0.7:
                        ins_len += 1
                        read.append(random.choice("ACGT"))
                        r_ins = random.uniform(0, 1)
                        qual.append(round(-math.log(0.7, 10) * 10))

            else:
                if was_del:  # adding uncertainty from previous deleted base
                    read.append(n)
                    qual.append(round(-math.log(was_del, 10) * 10))
                else:
                    read.append(n)
                    if p_error > 0:
                        qual.append(round(-math.log(p_error, 10) * 10))
                    else:
                        qual.append(1)   
                was_del = False
        if not read:
            continue
        read_seq = "".join([n for n in read])
        qual_seq = "".join([chr(q + 33) for q in qual])
        reads[i_acc] = (read_seq, qual_seq)
    if is_fastq:
            for acc, (read_seq, qual_seq) in sorted(reads.items(), key=lambda x: len(x[1]), reverse=True):
                outfile.write("@sim|err|{0}\n{1}\n{2}\n{3}\n".format(acc, read_seq, "+", qual_seq))
    else: #write in fasta format
            for acc, (read_seq, qual_seq) in sorted(reads.items(), key=lambda x: len(x[1]), reverse=True):
                outfile.write(">{0}\n{1}\n".format(acc, read_seq))
    outfile.close()
    return(reads)


def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def read_fragments_csv_file(filename):
    fragment_dict = {}
    print("Reading the input (.csv) file")
    with open(filename, 'r') as theFile:
        reader = csv.reader(theFile)
        next(reader)
        fragment_cter = 0
        for line in reader:
            print("Line[0]: "+line[0]+", Line[1]"+ line[1])
            fragment_dict[fragment_cter] = (line[0], line[1].replace("/5Phos/",""))
            fragment_cter += 1
            print(line)
    return fragment_dict


def generate_reads(args):
    print("Generating the reads")
    fragments_dict = read_fragments_csv_file(args.fragments)
    print(fragments_dict)
    max_nr = len(fragments_dict.items())-1
    reads_out = open(os.path.join(args.outfolder, 'clean_reads.fasta'), "w")
    read_infos = open(os.path.join(args.outfolder, 'readinfos.txt'), "w")
    reads = {}
    
    for i in range(0, args.nr_reads):
        startpos=0
        sequence = ''
        readinfos = []
        for j in range(0, args.nr_frags_per_read):
            fragment_id = random.randint(0, max_nr)
            sequence += fragments_dict[fragment_id][1]
            endpos = startpos + len(fragments_dict[fragment_id][1])
            readinfos.append((fragments_dict[fragment_id][0],startpos,endpos))
            startpos=endpos+1
        reads[i] = sequence
        reads_out.write(">sim|sim|{0}\n{1}\n".format(i, sequence))
        read_infos.write(">sim|sim|{0}: {1}\n".format(i, readinfos))

    return reads
    
def generate_fastq_from_reads(reads: list, outfolder: str, error_profile: list):
    print(reads)
    """Generates a FASTQ file from given reads, adding errors based on an error profile."""
    with open(os.path.join(args.outfolder,"simulated.fastq"), "w") as f:
        for i, sequence in reads.items():
            print(sequence)
            mutated_sequence, quality = simulate_errors(sequence, error_profile)
            f.write(f"@sim_rand_{i}\n")
            f.write(f"{mutated_sequence}\n")
            f.write(f"+\n")
            f.write(f"{quality}\n")

def main(args):
    mkdir_p(args.outfolder)
    error_profile = [0.51544048, 0.4633226, 0.41962698, 0.36491346, 0.37311439, 0.38849357, 0.41621851, 0.41683087, 0.41575156, 0.40372351, 0.3887138, 0.36476175,
                     0.33501046, 0.30455581, 0.27064415, 0.23457386, 0.20131726, 0.17431225, 0.15349413, 0.13793746, 0.12560718, 0.11617571, 0.10932199, 0.10399,
                     0.09931417, 0.09415376, 0.08774131, 0.07982219, 0.07201981, 0.06512025, 0.05915216, 0.05341532, 0.04848053, 0.0442655, 0.04047693, 0.03681451,
                     0.0334961, 0.03060112, 0.02820919, 0.02613725, 0.02442628, 0.02307134, 0.02214342, 0.02153789, 0.02130499, 0.0213524, 0.02154084, 0.02191075,
                     0.0223825, 0.02296152, 0.02362494, 0.02419089, 0.02464837, 0.02488377, 0.02482372, 0.02460933, 0.02430405, 0.02389571, 0.0234503, 0.02297452,
                     0.02250101, 0.02202724, 0.02161443, 0.02129141, 0.02116458, 0.02108961, 0.02120504, 0.02142005, 0.02169635, 0.02199924, 0.02225144, 0.02252378,
                     0.02270676, 0.02284495, 0.02293972, 0.022927, 0.02272802, 0.02245576, 0.02213147, 0.02186582, 0.02155404, 0.0212303, 0.02098836, 0.02080077,
                     0.02066826, 0.02060962, 0.02061027, 0.02063234, 0.02074351, 0.02086724, 0.02098664, 0.02112046, 0.02119616, 0.02121688, 0.02117194, 0.0211246,
                     0.02104943, 0.02096742, 0.02089215, 0.02085319]
    print("Generating "+str(args.nr_reads)+" different reads")
    #isoforms=generate_isoforms(args, genome_out.name)
    reads = generate_reads(args)
    generate_fastq_from_reads(reads, args.outfolder, error_profile)     
    #simulate_reads(args, reads)
            #print("Simulating reads")
    sys.exit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate reads for our experiments")
    parser.add_argument('--fragments', type=str,
                        help='Path to csv file with a nucleotide sequence containing the fragments to be used for the simulation')
    parser.add_argument('--nr_reads', type=int, default=200, help='Number of reads we want to simulate')
    parser.add_argument('--outfolder', type=str, help='Outfolder.')
    parser.add_argument('--nr_frags_per_read', type = int, default = 10, help = 'Number of fragments we want to add to our read')
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    main(args)
