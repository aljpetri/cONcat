import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
def parse_runtime(runtime):
    """Convert runtime to seconds if in minutes:seconds format."""
    if ':' in str(runtime):
        minutes, seconds = map(float, runtime.split(':'))
        return minutes * 60 + seconds
    return float(runtime)
    

def plot_fragments(data):
  # Adjust colors as needed

    df = pd.DataFrame(data)
    fig, ax = plt.subplots(figsize=(12, 6))
    error_rates = df['Error Rate'].unique()
    identities = df['Identity'].unique()
    bar_width = 0.15
    group_width = bar_width * len(identities) + 0.5

    positions = np.arange(len(error_rates)) * group_width

    for i, identity in enumerate(identities):
        subset = df[df['Identity'] == identity]
        ax.bar(positions + i * (bar_width + 0.05), subset['# frags detected'], bar_width, color='green')
        ax.bar(positions + i * (bar_width + 0.05), subset['# frags missed'], bar_width, bottom=subset['# frags detected'], color='red')

    ax.set_ylabel('Number of Fragments')
    ax.set_xticks(positions + bar_width * len(identities) / 2)
    ax.set_xticklabels([f'{er}' for er in error_rates], rotation=0, ha='center')

    for pos, er in zip(positions, error_rates):
        for i, iden in enumerate(identities):
            ax.text(pos + i * (bar_width + 0.05), -500, str(iden), ha='center', va='top', rotation=45)

    mid_point = (positions[0] + positions[-1]) / 2
    ax.text(mid_point, -600, 'Min_Identity', ha='center', va='top', rotation=45, fontsize=10)
    ax.text(mid_point + 0.3, -200, 'Error Rate', ha='center', va='top', rotation=0, fontsize=10)  # Slightly shifted right

    ax.legend(['Detected', 'Missed'], loc='lower left')
    plt.title('Fragments Detected and Missed Grouped by Error Rate')
    plt.tight_layout(rect=[0, 0.0, 1, 1])  # Adjust layout to avoid overlap
    plt.show()

def plot_data(file_path):
    palette = {0.25: 'darkorange', 0.50: 'blue', 0.75: 'green'}  # Adjust colors as needed
    markers = {0.25: 's', 0.50: 'o', 0.75: 'D', 1.00: 'X'}  # Different markers for each Identity


    # Load the data
    df = pd.read_csv(file_path, sep='\t')
    plot_fragments(df)
    # Convert runtime column
    df['Runtime(s)'] = df['Runtime(s)'].apply(parse_runtime)
    
    # Plot #Reads fully mapped vs. Error Rate
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df, x='Error Rate', y='#Reads fully mapped', hue='Identity',style='Identity', markers=markers,markersize =10, palette= palette)
    plt.title('#Reads Fully Covered vs. Error Rate')
    plt.xlabel('Error Rate')
    plt.ylabel('#Reads Fully Coverd')
    plt.legend(title='Identity')
    plt.grid(True)
    plt.show()
    
    # Plot Runtime and Memory usage separately for each Identity
    for metric in ['Runtime(s)', 'Memory(kbytes)']:
        plt.figure(figsize=(10, 6))
        sns.lineplot(data=df, x='Error Rate', y=metric, hue='Identity',style='Identity', markers=markers,markersize =10, palette= palette)
        plt.title(f'{metric} vs. Error Rate')
        plt.xlabel('Error Rate')
        plt.ylabel(metric)
        plt.legend(title='Identity')
        plt.grid(True)
        plt.show()

if __name__ == "__main__":
    plot_data('Results.tsv')
