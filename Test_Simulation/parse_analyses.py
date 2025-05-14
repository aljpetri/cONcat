import os
import re
import ast
import collections
import csv

from pathlib import Path
import os
import re
import ast
import collections
import csv
from pathlib import Path
from collections import Counter, defaultdict

def parse_analysis_files(folder_path):
    folder = Path(folder_path)
    # Match floats for x, y, z in filenames
    pattern = re.compile(r"out_analysis_([+-]?[0-9]*\.?[0-9]+)_([+-]?[0-9]*\.?[0-9]+)_([+-]?[0-9]*\.?[0-9]+)\.txt")
    run_summaries = []
    stats_by_z = defaultdict(lambda: defaultdict(lambda: {
        'original_count': 0,
        'other_count': 0,
        'missed': 0,
        'unexpected': 0
    }))

    for file in folder.iterdir():
        match = pattern.match(file.name)
        if not match:
            continue

        z_val = match.group(3)  # z is the third capture group
        with open(file, "r") as f:
            lines = f.readlines()

        original = other = None

        for line in lines:
            if line.startswith("Original:"):
                original = ast.literal_eval(line.split("Original:")[1].strip())
            elif line.startswith("Other:"):
                other = ast.literal_eval(line.split("Other:")[1].strip())

                # Process after both lists are obtained
                original_counts = Counter(original)
                other_counts = Counter(other)

                for frag in original:
                    stats_by_z[z_val][frag]['original_count'] += 1
                for frag in other:
                    stats_by_z[z_val][frag]['other_count'] += 1

                all_frags = set(original_counts.keys()).union(other_counts.keys())
                for frag in all_frags:
                    missed = original_counts[frag] - other_counts[frag]
                    unexpected = other_counts[frag] - original_counts[frag]
                    if missed > 0:
                        stats_by_z[z_val][frag]['missed'] += missed
                    if unexpected > 0:
                        stats_by_z[z_val][frag]['unexpected'] += unexpected

                original = other = None
        match = pattern.match(file.name)
        if not match:
            continue

        x_val, y_val, z_val = match.groups()

        with open(file, "r") as f:
            lines = f.readlines()

        # Extract from end of file
        reads_fully_preserved = None
        fragments_missed = None
        fragments_detected = None

        for line in reversed(lines):
            if line.startswith(" # fragments detected:") and fragments_detected is None:
                fragments_detected = int(line.strip().split(":")[1])
            elif line.startswith(" # fragments missed:") and fragments_missed is None:
                fragments_missed = int(line.strip().split(":")[1])
            elif line.startswith("# reads fully preserved:") and reads_fully_preserved is None:
                 reads_fully_preserved = int(line.strip().split(":")[1])

        run_summaries.append({
            "x": x_val,
            "y": y_val,
            "z": z_val,
            "reads_fully_preserved": reads_fully_preserved,
            "fragments_missed": fragments_missed,
            "fragments_detected": fragments_detected
})

    return stats_by_z, run_summaries
    
def write_run_summary_csv(run_summaries,outfolder):
    output_file="run_summaries.csv"
    with open(os.path.join(outfolder, output_file), "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=[
            "x", "y", "z",
            "reads_fully_preserved",
            "fragments_missed",
            "fragments_detected"
        ])
        writer.writeheader()
        for entry in run_summaries:
            writer.writerow(entry)
    print(f"Run summary written to {output_file}")
    
def write_summary_csv_per_z(stats_by_z, outfolder):
    output_prefix="fragment_summary_z="
    for z_val, frag_stats in stats_by_z.items():
        filename = f"{output_prefix}{z_val}.csv"
        file_full=os.path.join(outfolder, filename)
        with open(file_full, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Fragment", "Original Count", "Other Count", "Missed Count", "Unexpected Count"])
            for frag, counts in sorted(frag_stats.items()):
                writer.writerow([
                    frag,
                    counts['original_count'],
                    counts['other_count'],
                    counts['missed'],
                    counts['unexpected']
                ])
        print(f"Written summary for z={z_val} to {filename}")
        
def main(args):
    print("Folder",args.folder)
    stats_by_z, run_summaries = parse_analysis_files(args.folder)
    write_summary_csv_per_z(stats_by_z,args.folder)
    write_run_summary_csv(run_summaries,args.folder)
    
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Parse analysis files and extract fragment statistics grouped by z.")
    parser.add_argument("folder", help="Folder containing out_analysis_x_y_z.txt files")
    args = parser.parse_args()
    
    print(args.folder)
    main(args)
