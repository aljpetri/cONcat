mod side_functions;
mod structs;

use std::fs;
use rayon::prelude::*;
use csv::ReaderBuilder;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
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


fn write_outfile_from_vec(outfolder: &str, outfile_vector: Vec<(String, Vec<(String, usize, usize, i32)>)>, id_len_map: FxHashMap<String,usize>){
    fs::create_dir_all(outfolder).expect("The outfolder does not exist");
    let mut main_name = PathBuf::from(outfolder);
    let mut other_name = PathBuf::from(outfolder);
    main_name.push("fragments.txt");
    other_name.push("covering.txt");
    let f = File::create(main_name).expect("unable to create file");
    let f2 =File::create(other_name).expect("unable to create file");
    let mut buf_write = BufWriter::new(&f);
    let mut buf_write2 = BufWriter::new(&f2);


    write!(buf_write ,"read_acc,fragment_id,start,stop,edit_distance\n").expect("We should be able to write the entries");
    write!(buf_write2,"read_acc, length, #covered, %coverage\n").expect("We should be able to write the entries");
    for entry in outfile_vector{
        let mut bases_covered= 0;
        let read_header = entry.0;
        for fragment in entry.1{
            bases_covered = bases_covered + (fragment.2-fragment.1);
            write!(buf_write ,"{}, {}, {}, {}, {}\n", read_header, fragment.0, fragment.1, fragment.2, fragment.3).expect("We should be able to write the entries");


        }
        let readlen= id_len_map.get(read_header.as_str()).unwrap();
        write!(buf_write2,"{}, {}, {}, {}\n",read_header,readlen,bases_covered, (bases_covered as f64/ *readlen as f64) * 100_f64 ).expect("REASON")
    }

    buf_write.flush().expect("Failed to flush the buffer");
}
fn write_outfile(filename: String, outfile_hashmap: FxHashMap<String,Vec<(String, usize, usize, i32)>>){

    let f = File::create(filename).expect("unable to create file");

    let mut buf_write = BufWriter::new(&f);

    write!(buf_write ,"read_acc,fragment_id,start,stop,edit_distance,matches,fragment_length\n").expect("We should be able to write the entries");

    for entry in outfile_hashmap{
        let read_header= entry.0;
        for fragment in entry.1{

            write!(buf_write ,"{}, {}, {}, {}, {}\n", read_header, fragment.0, fragment.1, fragment.2, fragment.3).expect("We should be able to write the entries");

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
fn collect_sep_positions(sequence:&[u8], k:i32, expected_fragments_vec:&Vec<(String,String)>, previous_frags: usize,min_identity: f64, verbose: bool) ->Vec<(i32,i32,i32,String,usize,usize)>{
    if verbose{
        println!("prev_frags {}",previous_frags);
    }

    let mut sep_positions= vec![];
    let mut config = EdlibAlignConfigRs::default();
    config.mode = EdlibAlignModeRs::EDLIB_MODE_HW;
    config.task = EdlibAlignTaskRs::EDLIB_TASK_PATH;
    config.k = k;
    // best_alignments stores the best hits for each fragment, so we do not have to recompute every alignment
    let mut best_results = vec![];
    let mut best_identity= 0f64;
    let mut best_len= sequence.len();
    let mut best_frag: &str = "";
    let mut longest_frag= 0;
    let mut best_start= 0;
    let mut best_end= 0;
    let mut best_frag_str= best_frag.to_string();
    let mut best_alignment: Vec<u8> = vec![];
    let mut best_dist= 1000;
    //for all fragments use edlib and get the best
    for fragment in expected_fragments_vec {
        let fragment_name = fragment.0.clone();
        println!("Fragment: {}",fragment_name);
        let fragment_seq = fragment.1.as_bytes();
        let align_res = edlibAlignRs(fragment_seq, sequence, &config);
        if align_res.endLocations.is_some() {
            let end_locs = align_res.endLocations.clone().unwrap();
            //println!("Ends: {:?}",end_locs);
            //let mut start_locs=vec![];
            let start_locs = align_res.startLocations.clone().unwrap();
            for (loc_idx,end_loc) in end_locs.iter().enumerate() {
                let alignment_start= *start_locs.get(loc_idx).unwrap();
                let alignment_len:usize = (end_loc - alignment_start) as usize;
                let this_identity= (alignment_len as i32 - align_res.editDistance) as f64/ alignment_len as f64;
                println!("identity: {}, startpos {}, endpos {}",this_identity,alignment_start,end_loc);
                if this_identity > min_identity && this_identity > best_identity{
                    best_frag = &fragment_name;
                    best_frag_str = best_frag.to_string();
                    best_dist = align_res.editDistance;
                    best_identity = this_identity;
                    best_len = alignment_len;
                    best_start = alignment_start;
                    best_alignment = align_res.getAlignment().unwrap().clone();
                    //println!("Best identity {}",best_identity);
                    best_end = *end_loc as usize;
                    let best_alignment_frag = (best_frag_str.clone(), best_dist, alignment_start, best_end);
                    best_results.push(best_alignment_frag);
                }
            }
        }
    }
    if best_identity > 0f64{
        println!("best_fragment {}, start {}  end {} ED {}", best_frag_str, best_start, best_end, best_dist);
        println!("BFS: {}",best_frag_str);
        let best_total_pos= best_end + previous_frags;
        let sep_pos= ( (best_start + previous_frags as i32) as i32, (best_end + previous_frags) as i32,  best_dist, best_frag_str, longest_frag,best_total_pos);
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
    outfolder: String,
    #[arg(long,help="print additional information")]
    verbose: bool,
    #[arg(long,help="identity threshold used for the data (standard: 0.9) ")]
    identity_threshold: Option<f64>,
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


    #[arg(long,help="Run the post clustering step during the analysis (small improvement for results but much higher runtime)")]
    post_cluster: bool,
    #[arg(long,help="Do not write the fastq_files (no write_fastq in isONclust1)")]
    no_fastq: bool,
    #[arg(long,help="Minimum overlap threshold for reads to be clustered together (Experimental parameter)")]
    min_shared_minis: Option<f64>*/
}



fn main() {
    let cli = Cli::parse();
    let mut min_identity= 0.75;
    if cli.identity_threshold.is_some(){
        min_identity = cli.identity_threshold.unwrap();
    }
    println!("Min identity: {}",min_identity);
    let outfolder= cli.outfolder;// ="/home/alexanderpetri/Project3/SimulationResults/first4_fastqs/Alex/outfile";
    let expected_fragments_filename= cli.expected;//"/home/alexanderpetri/Project3/Expected_fragments.csv";
    let verbose=cli.verbose;
    let mut expected_fragments_vec: Vec<(String,String)> = vec![];
    expected_fragments_vec = read_csv_to_map(expected_fragments_filename);
    add_rev_comp_frags(&mut expected_fragments_vec); //comment this line out for current expected_fragments.csv
    //add_bw_frags(&mut expected_fragments_vec);
    if verbose{
        println!("{:?}",expected_fragments_vec);
    }

    expected_fragments_vec.sort_by(|a, b|b.1.len().partial_cmp(&a.1.len()).expect("REASON"));
    if verbose{
        println!("{:?}",expected_fragments_vec);
    }
    let long_fraglen= expected_fragments_vec.first().unwrap().1.len();
    let short_fraglen= expected_fragments_vec.last().unwrap().1.len();
    let filename= cli.fastq;// "/home/alexanderpetri/Project3/Fastqs/1_first_4.fastq";
    let k_inter:f64 = long_fraglen as f64 * (1f64 / 1f64 - min_identity);
    let k= k_inter.ceil() as i32;
    let mut frag_positions: Vec<(i32,i32,i32,String,usize,usize)>=vec![];
    let previous_frags: usize = 0;
    if verbose{
        println!("fragments range from {} to {} bp length", short_fraglen, long_fraglen);

    }
    let mut reader = fastq::Reader::from_file(Path::new(&filename)).expect("We expect the file to exist");
    let mut config = EdlibAlignConfigRs::default();
    let mut range_end = 0;
    let mut aligned_seg;
    let mut this_pos= vec![];
    //let mut outfile_hashmap = FxHashMap::default();
    let mut outfile_vector =vec![];
    config.mode = EdlibAlignModeRs::EDLIB_MODE_HW;
    config.k = k;
    let mut lendict=FxHashMap::default();
    //iterate over the reads in the fastq file
    for record in reader.records() {
        let mut covering_vec = vec![];
        let mut prev_end = 0;
        let seq_rec = record.expect("invalid record");
        let header_new = seq_rec.id();

        let mut sequence = seq_rec.seq();

        let mut seq_cover_vec: Vec<Interval> = vec![];
        lendict.insert(header_new.to_string(), sequence.len());
        if verbose{
            println!("{}", header_new);
            println!("Readlength {}", sequence.len());
            println!("Sequence {:?}",std::str::from_utf8(sequence));
        }

        //part_frags_map maps the fragment position to the part it was found in as wee need to keep track of the parts.
        let mut part_frags_map: FxHashMap<u64,(String, Vec<(i32,i32,i32,String,usize,usize)>)> = FxHashMap::default();
        //get the initial fragment positions
        frag_positions = collect_sep_positions(sequence, k, &expected_fragments_vec, previous_frags, min_identity, verbose).clone();
        part_frags_map.insert(0,(std::str::from_utf8(sequence).unwrap().to_string(), frag_positions.clone()));
        //we iterate until we do not have any more frag_positions to investigate
        while !frag_positions.is_empty()  {
            let mut frag_pos_to_delete= vec![];
            //iterate over the separator positions to find new separators First iteration: All coordinates are the true coordinates on the read, TODO: get correct coords for any further iterations
            for frag_pos in frag_positions.clone().into_iter() {
                let endpos = frag_pos.1 as usize;
                let frag_name = frag_pos.3.to_owned();
                let frag_len = frag_pos.4.to_owned();
                let range_start = frag_pos.0 as usize;
                let edit_dist = frag_pos.2;


                range_end = endpos + 1;
                //get the sequence and length of the bases located between two detected fragments
                let sub_range_seq = &sequence[prev_end..range_start];
                let sub_range_len = sub_range_seq.len();
                if verbose{
                    println!("Frag start: {}, frag end {} , edit distance {}, frag {}", range_start, range_end, edit_dist, frag_name);
                }
                //add the fragment information to our final data structure covering_vec
                covering_vec.push((frag_name, range_start, range_end, edit_dist));
                if verbose{

                    println!("frag_pos_0: {}, frag_pos_3: {}", endpos, frag_len);
                    println!("Range starts at {}, length {}", prev_end, sub_range_len);
                    println!("before: {:?}", covering_vec);
                    println!("after: {:?}", covering_vec);
                }
                frag_pos_to_delete.push(frag_pos);
                aligned_seg = std::str::from_utf8(&sequence[range_start..range_end]).unwrap();
                let tmp_inter = Interval { start: range_start, end: range_end };
                seq_cover_vec.push(tmp_inter);

            }
            if verbose{
                println!("fragments to delete: {:?}",frag_pos_to_delete);
            }
            for frag_to_delete in frag_pos_to_delete.into_iter(){
                frag_positions.retain(|x| *x != frag_to_delete);
            }

            let mut int_start = 0 as usize;
            seq_cover_vec.sort_by_key(|interval| interval.start);

            if verbose{

                println!("{:?}",seq_cover_vec);
                println!("fragments after delete: {:?}",frag_positions);
            }
            let mut interlen;
            let mut add_ranges= vec![];
            for part in &seq_cover_vec {
                if part.start > int_start{
                    interlen = part.start - int_start;
                }
                else {
                    interlen = 0;
                }
                if verbose{
                    println!("Range from {} to {}", int_start, part.start);
                }

                if interlen > 5 {
                    this_pos = collect_sep_positions(&sequence[int_start..part.start], k, &expected_fragments_vec, int_start, min_identity, verbose);
                    if verbose{
                        println!("Added this pos: {:?}", this_pos);
                    }

                    if this_pos.is_empty(){
                        if verbose{
                            println!("Empty Interval: {}, {}",int_start,part.start)
                        }
                        add_ranges.push( Interval {start: int_start,end: part.start})
                    }
                    else {
                        frag_positions.append(&mut this_pos);
                    }

                }
                int_start = part.end;
            }
            seq_cover_vec.append(&mut add_ranges);
            if verbose{
                println!("int_start after for loop {}", int_start);
                println!("Range from {} to {}", int_start, sequence.len());
            }


            if sequence.len() - int_start > 5{
                this_pos = collect_sep_positions(&sequence[int_start..sequence.len()], k, &expected_fragments_vec, int_start, min_identity, verbose);
                if this_pos.is_empty(){
                    seq_cover_vec.push( Interval {start: int_start,end: sequence.len()})
                }
                else {
                    frag_positions.append(&mut this_pos);
                }
                if verbose{
                    println!("Added this pos: {:?}", this_pos);
                }

                frag_positions.append(&mut this_pos);
            }
        }
        covering_vec.sort_by_key(|interval| interval.1);
        if verbose{
            println!("ReadID {}",header_new);

        }
        println!("covering: {:?}", covering_vec);
        let outfile_object = (header_new.to_string(),covering_vec);
        outfile_vector.push(outfile_object);
    }
    write_outfile_from_vec(&outfolder, outfile_vector,lendict)
}
