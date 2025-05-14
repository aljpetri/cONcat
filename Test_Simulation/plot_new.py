import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
def parse_runtime(runtime):
    if ':' in str(runtime):
        minutes, seconds = map(float, runtime.split(':'))
        return minutes * 60 + seconds
    return float(runtime)

def plot_fragments(df,folder):
    
    fig, ax = plt.subplots(figsize=(12, 6))
    error_rates = sorted(df['Error Rate'].unique())  # low to high
    identities = sorted(df['y'].unique())
    bar_width = 0.15
    group_width = bar_width * len(identities) + 0.5
    positions = np.arange(len(error_rates)) * group_width

    for i, identity in enumerate(identities):
        for j, er in enumerate(error_rates):
            subset = df[(df['y'] == identity) & (df['Error Rate'] == er)]
            if subset.empty:
                continue
            detected_mean = subset['Fragments Detected'].mean()
            missed_mean = subset['Fragments Missed'].mean()
            total = detected_mean + missed_mean if (detected_mean + missed_mean) > 0 else 1  # avoid divide-by-zero
            percent_detected = 100 * detected_mean / total
            percent_missed = 100 * missed_mean / total

            ax.bar(positions[j] + i * (bar_width + 0.05), percent_detected, bar_width, color='green')
            ax.bar(positions[j] + i * (bar_width + 0.05), percent_missed, bar_width, bottom=percent_detected, color='red')

    ax.set_ylabel('Fragments Coverage (%)')
    ax.set_xticks(positions + bar_width * len(identities) / 2)
    ax.set_xticklabels([f'{er:.0f}' for er in error_rates], rotation=0, ha='center')

    for pos, er in zip(positions, error_rates):
        for i, iden in enumerate(identities):
            ax.text(pos + i * (bar_width + 0.05), -5, str(iden), ha='center', va='top', rotation=45)

    mid_point = (positions[0] + positions[-1]) / 2
    ax.text(mid_point, -6, 'Min_Identity', ha='center', va='top', rotation=45, fontsize=10)
    ax.text(mid_point + 0.3, -2, 'Error Rate (%)', ha='center', va='top', rotation=0, fontsize=10)

    ax.legend(['Detected (%)', 'Missed (%)'], loc='upper right')
    ax.set_ylim(0, 105)
    plt.title('Percentage of Fragments Detected and Missed Grouped by Error Rate')
    plt.tight_layout(rect=[0, 0.0, 1, 1])
    #plt.show()
    outfile = os.path.join(args.folder,'fragments.png')
    plt.savefig(outfile, dpi=300)

def plot_with_errorbars(df, y, ylabel, palette, markers, folder):
    plt.figure(figsize=(10, 6))
    grouped = df.groupby(['Error Rate', 'y'])[y]
    
    # Calculate mean and std for error bars
    stats = grouped.agg(['mean', 'std']).reset_index()
    
    for identity in df['y'].unique():
        subset = stats[stats['y'] == identity]
        plt.errorbar(
            subset['Error Rate'],
            subset['mean'],
            yerr=subset['std'],
            label=f'Identity {identity}',
            fmt=markers.get(identity, 'o'),
            capsize=5,
            color=palette.get(identity, 'black'),
            linestyle='-'
        )
    
    plt.title(f'{ylabel} vs. Error Rate')
    plt.xlabel('Error Rate (%)')
    plt.ylabel(ylabel)
    plt.legend(title='Identity')
    plt.grid(True)
    #plt.show()
    outfile = os.path.join(folder,'ReadsMapping.png')
    plt.savefig(outfile, dpi=300)

def plot_data(args):
    infile=os.path.join(args.folder,'run_summaries.csv')
    palette = {0.25: 'darkorange', 0.50: 'blue', 0.75: 'green'}
    markers = {0.25: 's', 0.50: 'o', 0.75: 'D', 1.00: 'X'}

    df = pd.read_csv(infile)
    df.rename(columns={
        'reads_fully_preserved': '#Reads Fully Preserved',
        'fragments_missed': 'Fragments Missed',
        'fragments_detected': 'Fragments Detected'
    }, inplace=True)

    # Convert Correct Rate to Error Rate
    df['Error Rate'] = (1 - df['x']) *100
    # Ensure numeric types
    for col in ['#Reads Fully Preserved', 'Fragments Missed', 'Fragments Detected']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    plot_fragments(df,args.folder)

    # Plot individual points with error bars
    plot_with_errorbars(df, '#Reads Fully Preserved', '#Reads Fully Preserved', palette, markers,args.folder)

    # Optional: Plot Runtime and Memory if available
    for metric in ['Runtime(s)', 'Memory(kbytes)']:
        if metric not in df.columns:
            continue
        df[metric] = df[metric].apply(parse_runtime)
        plot_with_errorbars(df, metric, metric, palette, markers, args.folder)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Parse analysis files and extract fragment statistics grouped by z.")
    parser.add_argument("folder", help="Folder containing out_analysis_x_y_z.txt files")
    args = parser.parse_args()
    plot_data(args)
    print(args.folder)
    
