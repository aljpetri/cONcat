mod side_functions;
mod structs;

use csv::ReaderBuilder;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use bio::io::fastq;
use edlib_rs::edlibrs::*;
use itertools::Itertools;
use rustc_hash::FxHashMap;
use crate::side_functions::reverse_complement;
use crate::structs::Interval;
use clap::Parser;

fn read_csv_to_map(filename: String) -> Vec<(String, String)> {
    // Open the file
    let file = File::open(filename).unwrap();
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_reader(BufReader::new(file));
    // Initialize the FxHashMap
    let mut expected_fragments_map: Vec<(String, String)> = vec![];
    // Iterate over each record in the CSV file
    for result in rdr.records() {
        let record = result.unwrap();
        if record.len() != 2 {
            eprintln!("Skipping invalid record: {:?}", record);
            continue;
        }
        let seqname= record.get(0).unwrap();
        //let seqname_o=seqname.as_bytes();
        let mut sequence = record.get(1).unwrap();
        if sequence.as_bytes().starts_with(b"/5Phos/"){
            sequence = sequence.strip_prefix("/5Phos/").unwrap();
        }
        let this_tuple=(seqname.to_string(),sequence.to_string());
        expected_fragments_map.push(this_tuple);
    }
    expected_fragments_map
}
fn write_outfile(filename: String, outfile_hashmap: FxHashMap<String,Vec<(String, usize, usize, i32)>>){
    let f = File::create(filename).expect("unable to create file");
    let mut buf_write = BufWriter::new(&f);
    write!(buf_write ,"read_acc,fragment_id,start,stop,aln_length,matches,fragment_length\n").expect("We should be able to write the entries");
    for entry in outfile_hashmap{
        let read_header=entry.0;
        for fragment in entry.1{
            write!(buf_write ,"{},{},{},{}\n", read_header, fragment.0,fragment.1,fragment.2).expect("We should be able to write the entries");
        }
    }

    buf_write.flush().expect("Failed to flush the buffer");
}

fn add_rev_comp_frags(expected_fragments_vec: &mut Vec<(String,String)>){
    let mut rev_comp_frags:Vec<(String,String)>=vec![];
    for fragment in expected_fragments_vec.iter(){
        let mut frag_name= fragment.0.clone();
        let frag_seq= &fragment.1;
        let frag_seq_rc = reverse_complement(frag_seq);
        frag_name.push_str("_rc");
        let rc_frag=(frag_name,frag_seq_rc.clone());
        rev_comp_frags.push(rc_frag);
    }
    expected_fragments_vec.append(&mut rev_comp_frags);
}


fn add_bw_frags (expected_fragments_vec: &mut Vec<(String,String)>){
    let mut rev_comp_frags:Vec<(String,String)>=vec![];
    for fragment in expected_fragments_vec.iter() {
        let mut frag_name = fragment.0.clone();
        let mut frag_seq = &fragment.1;
        let rev_seq=frag_seq.chars().rev().collect::<String>();
        frag_name.push_str("_bw");
        let rc_frag=(frag_name,rev_seq.clone());
        rev_comp_frags.push(rc_frag);
    }
    expected_fragments_vec.append(&mut rev_comp_frags);
}


//recursive function that is used to detect the positions of the seperator in the read
fn collect_sep_positions(sequence:&[u8], k:i32, expected_fragments_vec:&Vec<(String,String)>, previous_frags: usize,min_identity: f64) ->Vec<(i32,i32,i32,String,usize,usize)>{
    println!("prev_frags {}",previous_frags);
    let mut sep_positions= vec![];
    let mut config = EdlibAlignConfigRs::default();
    config.mode = EdlibAlignModeRs::EDLIB_MODE_HW;
    config.k = k;
    // best_alignments stores the best hits for each fragment, so we do not have to recompute every alignment
    let mut best_alignments = vec![];
    let mut best_dist= 1000;
    let mut best_frag: &str = "";
    let mut longest_frag= 0;
    let mut best_start= 0;
    let mut best_end= 0;
    let mut best_frag_str= best_frag.to_string();
    for fragment in expected_fragments_vec {
        let fragment_name = fragment.0.clone();
        //println!("Fragment: {}",fragment_name);
        let fragment_seq = fragment.1.as_bytes();
        let fragment_len = fragment_seq.len();
        let align_res = edlibAlignRs(fragment_seq, sequence, &config);
        if align_res.endLocations.is_some() {
            let end_locs = align_res.endLocations.unwrap();
            //println!("endlocs {:?}", end_locs);
            for end_loc in end_locs {
                //TODO::From here on we need to change from best_dist to best_identity, identity = (alignment length - edit distance)/fragment_length
                if align_res.editDistance < best_dist || align_res.editDistance == best_dist && fragment_len > longest_frag {
                    best_frag = &fragment_name;
                    best_frag_str = best_frag.to_string();
                    if end_loc as usize > fragment_len{
                        best_start = end_loc as usize - (fragment_len - 1); //TODO: this might be a bug: we need to better grasp how long the alignment is
                    }
                    else{
                        best_start = 0;
                    }
                    best_end = end_loc as usize + 1;
                    best_dist = align_res.editDistance;
                    longest_frag = fragment_len;
                    let best_alignment_frag = (best_frag_str.clone(), best_dist, best_start, best_end);
                    best_alignments.push(best_alignment_frag);
                }
            }
        }
    }
    if best_dist < 1000{
        println!("best_fragment {}, start {}  end {} ED {}", best_frag_str, best_start, best_end, best_dist);
        println!("BFS: {}",best_frag_str);
        let best_total_pos= best_end + previous_frags;
        let sep_pos= ( (best_start + previous_frags) as i32, (best_end + previous_frags) as i32,  best_dist, best_frag_str, longest_frag,best_total_pos);
        sep_positions.push(sep_pos);
    }
    sep_positions
}

#[derive(Parser,Debug)]
#[command(name = "DNA fragment covering")]
#[command(author = "Alexander J. Petri <alexjpetri@gmail.com>")]
#[command(version = "0.0.1")]
#[command(about = "Clustering of long-read sequencing data into gene families", long_about = "isONclust is a tool for clustering either PacBio Iso-Seq reads, or Oxford Nanopore reads into clusters, where each cluster represents all reads that came from a gene." )]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[arg(long, short,help="Path to expected fragments file in csv format")]
    expected: String,
    #[arg(long, short, help="Path to input fastq file")]
    fastq: String,
    #[arg(long, short, help="Path to input fastq file")]
    outfile: String
    /*#[arg(short,  help="Kmer length")]
    k: Option<usize>,
    #[arg(short, help=" window size")]
    w: Option<usize>,
    #[arg(short, help=" syncmer length")]
    s: Option<usize>,
    #[arg(short, help=" minimum syncmer position")]
    t: Option<usize>,
    #[arg(long, short, help="Path to outfolder")]
    outfile: String,
    #[arg(long,short,default_value_t= 1, help="Minimum number of reads for cluster")]
    n: usize,
    #[arg(long, short,help="Path to gff3 file (optional parameter), requires a reference added by calling --init-cl <REFERENCE.fasta>")]
    gff: Option<String>,
    #[arg(long,help="we do not want to use canonical seeds")]
    noncanonical: bool,
    #[arg(long,help="Run mode of isONclust (pacbio or ont")]
    mode: String,
    //TODO:add argument telling us whether we want to use syncmers instead of kmers, maybe also add argument determining whether we want to use canonical_minimizers
    #[arg(long,help="seeding approach we choose")]
    seeding: Option<String>,
    #[arg(long,help="quality threshold used for the data (standard: 0.9) ")]
    quality_threshold: Option<f64>,
    #[arg(long,help="print additional information")]
    verbose: bool,
    #[arg(long,help="Run the post clustering step during the analysis (small improvement for results but much higher runtime)")]
    post_cluster: bool,
    #[arg(long,help="Do not write the fastq_files (no write_fastq in isONclust1)")]
    no_fastq: bool,
    #[arg(long,help="Minimum overlap threshold for reads to be clustered together (Experimental parameter)")]
    min_shared_minis: Option<f64>*/
}



fn main() {
    let cli = Cli::parse();
    let min_identity= 0.75;
    let outfilename=cli.outfile;// ="/home/alexanderpetri/Project3/SimulationResults/first4_fastqs/Alex/outfile";
    let expected_fragments_filename= cli.expected;//"/home/alexanderpetri/Project3/Expected_fragments.csv";
    let mut expected_fragments_vec: Vec<(String,String)> = vec![];
    expected_fragments_vec = read_csv_to_map(expected_fragments_filename);
    add_rev_comp_frags(&mut expected_fragments_vec);
    add_bw_frags(&mut expected_fragments_vec);
    println!("{:?}",expected_fragments_vec);
    expected_fragments_vec.sort_by(|a, b|b.1.len().partial_cmp(&a.1.len()).expect("REASON"));
    println!("{:?}",expected_fragments_vec);
    let long_fraglen= expected_fragments_vec.first().unwrap().1.len();
    let short_fraglen= expected_fragments_vec.last().unwrap().1.len();
    let filename=cli.fastq;// "/home/alexanderpetri/Project3/Fastqs/1_first_4.fastq";
    let k_inter:f64 = long_fraglen as f64 * (1f64 / 1f64 - min_identity);
    let k= k_inter.ceil() as i32;
    let mut frag_positions: Vec<(i32,i32,i32,String,usize,usize)>=vec![];
    let previous_frags: usize = 0;
    println!("fragments range from {} to {} bp length",short_fraglen,long_fraglen);
    let mut reader = fastq::Reader::from_file(Path::new(&filename)).expect("We expect the file to exist");
    let mut config = EdlibAlignConfigRs::default();
    let mut range_end = 0;
    let mut aligned_seg;
    let mut this_pos=vec![];
    let mut outfile_hashmap = FxHashMap::default();
    config.mode = EdlibAlignModeRs::EDLIB_MODE_HW;
    config.k = 8;
    for record in reader.records() {
        let mut it_ct=0;
        let mut covering_vec = vec![];
        let mut prev_end = 0;
        let seq_rec = record.expect("invalid record");
        let header_new = seq_rec.id();
        println!("{}", header_new);
        let mut sequence = seq_rec.seq();
        let mut seq_cover_vec: Vec<Interval> = vec![];
        println!("Readlength {}", sequence.len());
        //part_frags_map maps the fragment position to the part it was found in as wee need to keep track of the parts.
        let mut part_frags_map: FxHashMap<u64,(String, Vec<(i32,i32,i32,String,usize,usize)>)> = FxHashMap::default();
        //get the initial fragment positions
        frag_positions = collect_sep_positions(sequence, k, &expected_fragments_vec, previous_frags, min_identity).clone();
        part_frags_map.insert(0,(std::str::from_utf8(sequence).unwrap().to_string(), frag_positions.clone()));
        //TODO: How to I only split the part sequences and not the overall read sequence anymore?
        while !frag_positions.is_empty()  && it_ct < 3 {
            it_ct += 1;
            let mut frag_pos_to_delete= vec![];
            //iterate over the separator positions to find new separators First iteration: All coordinates are the true coordinates on the read, TODO: get correct coords for any further iterations
            for frag_pos in frag_positions.clone().into_iter() {
                let endpos = frag_pos.1 as usize;
                let frag_name = frag_pos.3.to_owned();
                let frag_len = frag_pos.4.to_owned();
                let range_start = frag_pos.0 as usize;
                let edit_dist = frag_pos.2;
                println!("frag_pos_0: {}, frag_pos_3: {}", endpos, frag_len);
                range_end = endpos + 1;
                println!("Frag start: {}, frag end {} , edit distance {}, frag {}", range_start, range_end, edit_dist, frag_name);
                let sub_range_seq = &sequence[prev_end..range_start];
                let sub_range_len = sub_range_seq.len();
                println!("Range starts at {}, length {}", prev_end, sub_range_len);
                println!("before: {:?}", covering_vec);
                covering_vec.push((frag_name, range_start, range_end, edit_dist));
                println!("after: {:?}", covering_vec);
                frag_pos_to_delete.push(frag_pos);
                aligned_seg = std::str::from_utf8(&sequence[range_start..range_end]).unwrap();
                let tmp_inter = Interval { start: range_start, end: range_end };
                seq_cover_vec.push(tmp_inter);

            }
            println!("fragments to delete: {:?}",frag_pos_to_delete);
            for frag_to_delete in frag_pos_to_delete.into_iter(){
                frag_positions.retain(|x| *x != frag_to_delete);
            }
            println!("fragments after delete: {:?}",frag_positions);
            let mut int_start = 0 as usize;
            seq_cover_vec.sort_by_key(|interval| interval.start);
            println!("{:?}",seq_cover_vec);
            let mut interlen;
            let mut add_ranges=vec![];
            for part in &seq_cover_vec {
                if part.start > int_start{
                    interlen = part.start - int_start;
                }
                else {
                    interlen = 0;
                }
                println!("Range from {} to {}", int_start, part.start);
                if interlen > 5 {
                    this_pos = collect_sep_positions(&sequence[int_start..part.start], k, &expected_fragments_vec, int_start, min_identity);
                    println!("Added this pos: {:?}", this_pos);
                    if this_pos.is_empty(){
                        add_ranges.push( Interval {start: int_start,end: part.start})
                    }
                    else {
                        frag_positions.append(&mut this_pos);
                    }

                }
                int_start = part.end;
            }
            seq_cover_vec.append(&mut add_ranges);
            println!("int_start after for loop {}", int_start);
            println!("Range from {} to {}", int_start, sequence.len());
            if sequence.len() - int_start > 5{
                this_pos = collect_sep_positions(&sequence[int_start..sequence.len()], k, &expected_fragments_vec, int_start, min_identity);
                if this_pos.is_empty(){
                    seq_cover_vec.push( Interval {start: int_start,end: sequence.len()})
                }
                else {
                    frag_positions.append(&mut this_pos);
                }
                println!("Added this pos: {:?}", this_pos);
                frag_positions.append(&mut this_pos);
            }
        }
        covering_vec.sort_by_key(|interval| interval.1);
        println!("covering: {:?}", covering_vec);
        outfile_hashmap.insert(header_new.to_string(),covering_vec);
    }
    write_outfile(outfilename,outfile_hashmap);

}
