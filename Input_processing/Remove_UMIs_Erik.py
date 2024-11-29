from Bio import SeqIO
from Bio import Align

import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd


#### THIS code is simplifed, it only looks for the UMI by finding the FWD primer and then looking for the next 10 nt (should be the UMI) and likewise for the backward UMI)
def align_data(name_data,name_out,Primer_fwd_seq,Primer_rev_seq,Primer_after_second_UMI_fwd):

    UMI_for_saving = []

    initial_align_window = 90# the amount of the begining sequence used to find the first position of the primers

    read_sequences = []

    with open(name_data) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_sequences.append(record.seq)


    sequence_lenghts = [] #keep track of the reads inside the primers
    First_UMI = []
    Second_UMI = []
    Full_sequence = []
    Adapter_free_sequence = []
    Fused_UMI = []
    FWD_primers = []
    REV_primers = []

    for n, i in enumerate(read_sequences):

        sequence_used_in_round = i # I do this to be able to swap it with its reverse compement if its a backward read

        if len(i) < 80:
            print('strand was to short')
            continue

        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = 2
        aligner.mismatch_score = -1.
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -2
        alignments_fwd = aligner.align(sequence_used_in_round[0:initial_align_window], Primer_fwd_seq)
        alignments_bwd = aligner.align(sequence_used_in_round[0:initial_align_window], Primer_rev_seq)

        try:  # sometimes the sequence is just GGGG or something like that and the aligner finds no score at all causing issues with downstream assesemtn, thus try statement throws out the sequence if no score an be found
            if alignments_fwd[0].score < len(Primer_fwd_seq)*1.5 and alignments_bwd[0].score < len(Primer_fwd_seq)*1.5:  # i throw out initial alignemnts that are where i cannot get a good alignment with either primer, they may be other junk from sequencing
                print('found no alignment')
                continue
            else:
                pass
        except:
            continue

        if alignments_fwd[0].score < alignments_bwd[0].score:
            sequence_used_in_round = sequence_used_in_round.reverse_complement() # if the read starts with the backward primer, we are reading the wrong strand. make a reverse complemement and proceed
            alignments_fwd = aligner.align(sequence_used_in_round[0:initial_align_window], Primer_fwd_seq)
            print('trying to align to the end of the reverse complement')
            print(alignments_fwd[0])
            if alignments_fwd[0].score < len(Primer_fwd_seq)*1.5:
                print('was unable to align the end of the reverse complement')
                continue
            else:
                pass

            print('found fwd aligment')
            print(alignments_fwd[0])
        start_of_read = int(alignments_fwd[0].aligned[0][-1][1])

        print('first UMI is: ',sequence_used_in_round[start_of_read:start_of_read+10])

        alignments_fwd_end = aligner.align(sequence_used_in_round[-initial_align_window:], Primer_after_second_UMI_fwd )

        end_of_second_UMI = int(alignments_fwd_end[0].aligned[0][-1][0]) + len(sequence_used_in_round) - initial_align_window


        try:  # sometimes the sequence is just GGGG or something like that and the aligner finds no score at all causing issues with downstream assesemtn, thus try statement throws out the sequence if no score an be found
            if alignments_fwd_end[0].score < len(Primer_after_second_UMI_fwd) * 1.5:  # i throw out initial alignemnts that are where i cannot get a good alignment with either primer, they may be other junk from sequencing
                print('found no  end alignment')
                continue
            else:
                pass
        except:
            continue

        print('second UMI is: ', str(sequence_used_in_round[end_of_second_UMI-10:end_of_second_UMI]))

        #### Store data only if all steps are gone through successfully
        First_UMI.append(str(sequence_used_in_round[start_of_read:start_of_read+10]))
        Second_UMI.append(str(sequence_used_in_round[end_of_second_UMI-10:end_of_second_UMI]))
        Fused_UMI.append(str(sequence_used_in_round[start_of_read:start_of_read+10])+str(sequence_used_in_round[end_of_second_UMI-10:end_of_second_UMI]))

        #### HERE I figure out the primers that would be needed to amplify the species of interest
        FWD_primers.append(str(sequence_used_in_round[start_of_read-10:start_of_read+10]))

        rev_primer = sequence_used_in_round[end_of_second_UMI-10:end_of_second_UMI+10].reverse_complement()

        REV_primers.append(str(rev_primer))

        Full_sequence.append(str(sequence_used_in_round))
        Adapter_free_sequence.append(str(sequence_used_in_round[start_of_read:end_of_second_UMI]))

        sequence_lenghts.append(len(sequence_used_in_round[start_of_read:end_of_second_UMI]))

    # dictionary of lists
    dict = {'First UMI': First_UMI, 'Second UMI': Second_UMI,'Fused_UMI' : Fused_UMI, 'Full Sequence': Full_sequence, 'Adapter free sequence': Adapter_free_sequence, 'Forward primer needed': FWD_primers, 'Reverse primers needed': REV_primers, 'adapter free sequence length': sequence_lenghts}

    df = pd.DataFrame(dict)

    df['Total Occ of fused UMI'] = df.groupby('Fused_UMI').Fused_UMI.transform('count')

    df.to_csv(name_out)

    print('number of sequences processed: ',str(len(read_sequences)))
    UMI_len = []
    for i in First_UMI: UMI_len.append(len(i))

    sns.histplot(data=UMI_len)
    plt.show()

    UMI_len = []
    for i in Second_UMI: UMI_len.append(len(i))

    sns.histplot(data=UMI_len)
    plt.show()

    sns.histplot(data=sequence_lenghts)
    plt.xlim(0,1000)
    plt.show()


if __name__ == '__main__':

     ###### THESE are the primers when we have UMIs on Dengue  #####
    #Primer_fwd_seq = 'ACGTGTGTGTGCTATGATAC'
    #Primer_rev_seq = 'TGTACACACGCATACACACA'
    #hairpin_after_UMI_fwd = 'CCGCGATGATGCGATACGTGTAAG'
    #Primer_fwd_seq_end = 'CTTATATCTGTGTATGTGTGCAGG'
    #Primer_after_second_UMI_fwd = 'TGTGTGTATGCGTGTGTACA'

    ###### THESE are the primers when we have UMIs on Cell uptae #####
    # The full sequence for the UMI Loop is /5Phos/ACTCCATACGTATCTANNNNNNNNNNACATATACATATCACGTGCTGGCCTACACGCGTACATGTATAGANNNNNNNNNNTAGATACGTATGG

    # Primer_fwd_seq = 'CCAGCACGTGATATGTATATGT'
    # Primer_rev_seq = 'CCTACACGCGTACATGTATAGA'
    #
    # Primer_fwd_seq = 'ATATGTATATGT'
    # Primer_rev_seq = 'TACATGTATAGA'
    #
    #
    # hairpin_after_UMI_fwd = 'TAGATACGTATGGAGT'
    #
    # #hairpin_after_UMI_fwd = 'CCATACGTATCTA'
    # Primer_fwd_seq_end = 'TAGATACGTATGG'
    # Primer_fwd_seq_end = 'CCATACGTATCTA'
    # Primer_after_second_UMI_fwd = 'TCTATACATGTACGCGTGTAGG'


    #align_data('DengueUMI_r4.fastq')
    align_data('HEK293_R8.fastq','HEK293_R8_simplified.csv','CCAGCACGTGATATGTATATGT','CCTACACGCGTACATGTATAGA','TCTATACATGTACGCGTGTAGG')