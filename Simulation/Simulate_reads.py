import csv
import os, sys
import random
import argparse
import errno
import math


#adds errors to the reads
def simulate_reads(args, isoforms):
    outfile = open(os.path.join(args.outfolder,"simulated.fastq"), "w")
    is_fastq = True

    reads = {}
    #different error rates we can add to our reads
    #error_lvls = [0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96, 0.98, 0.99, 0.995]
    #error_lvls=[0.9,0.95,0.96,0.98,0.99,0.995]#3.94%error rate
    error_lvls=[0.8, 0.875,0.9,0.92,0.96,0.98,0.99,0.995]#7%error rate
    #error_lvls=[0.97]
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
                    qual.append(round(-math.log(p_error, 10) * 10))

                else: #insertion
                    read.append(n)
                    qual.append(round(-math.log(p_error, 10) * 10))

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
                    qual.append(round(-math.log(p_error, 10) * 10))
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
        fragment_cter = 0
        for line in reader:
            fragment_dict[fragment_cter] = (line[0], line[1])
            fragment_cter += 1
            print(line)
    return fragment_dict


def generate_reads(args):
    print("Generating the reads")
    fragments_dict = read_fragments_csv_file(args.fragments)
    print(fragments_dict)
    max_nr = len(fragments_dict.items())-1
    reads_out = open(os.path.join(args.outfolder, 'clean_reads.fasta'), "w")
    read_infos=open(os.path.join(args.outfolder, 'readinfos.txt'), "w")
    reads = {}
    for i in range(0, args.nr_reads):
        sequence = ''
        readinfos = []
        for j in range(0, args.nr_frags_per_read):
            fragment_id = random.randint(0, max_nr)
            sequence += fragments_dict[fragment_id][1]
            readinfos.append(fragments_dict[fragment_id][0])
        reads[i] = sequence
        reads_out.write(">sim|sim|{0}\n{1}\n".format(i, sequence))
        read_infos.write(">sim|sim|{0}\n{1}\n".format(i, readinfos))

    return reads


def main(args):
    mkdir_p(args.outfolder)
    print("Generating "+str(args.nr_reads)+" different reads")
    #isoforms=generate_isoforms(args, genome_out.name)
    reads = generate_reads(args)
    print("Hello World")
    simulate_reads(args, reads)
            #print("Simulating reads")
    sys.exit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate reads for our experiments")
    parser.add_argument('--fragments', type=str,
                        help='Path to csv file with a nucleotide sequence containing the fragments to be used for the simulation')
    parser.add_argument('--nr_reads', type=int, default=200, help='Number of reads we want to simulate')
    parser.add_argument('--outfolder', type=str, help='Outfolder.')
    parser.add_argument('--nr_frags_per_read', type=int, default=0, help='Number of fragments we want to add to our read')
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    main(args)