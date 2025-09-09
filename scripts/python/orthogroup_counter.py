from Bio import SeqIO
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

# Count sequences per FASTA files
def count_sequences(folder_path):
    """
    Count Sequences per .faa FASTA file of the folder path.
    Scans in the folder files that end with .faa, parses each file 
    and counts the number of sequences. Returns a disctionary with 
    orthogroup IDs and sequence counts.

    Args:
    folder_path (str): Folder path with the .faa FASTA files

    Returns:
    List (dictionary, str): A list with orthogroup base filename 
    and sequence count per FASTA file.
    """
    counts = []
    for file in os.listdir(folder_path):
        if file.endswith(".faa"):   # Process FASTA files ending with .faa, change if needed.
            file_path = os.path.join(folder_path, file)
            seq_count = sum(1 for _ in SeqIO.parse(file_path, "fasta"))
            counts.append({
                "Orthogroup": file.replace(".faa", ""),
                "SequenceCount": seq_count
            })
    return counts

# DataFrame with renaming and ordering
def custom(counts, rename, order):
    """
    DataFrame, with renaiming and plotting order.
    The list is turned into a DataFrame, replaces the names of the Orthogroup
    IDs and adds a custom plotting order. 

    Args:
    counts(dictionary, str): Output list of count_sequences.
    rename (dictionary, str): New labels for Orthogroup IDs x-axis.
    custom (list, str): Order of orthogroup labels. 

    Returns:
    DataFrame with columns: Orthogroup (str) renamed, Sequence Counts (int) per orthogroup,
    orthogroup_order (category),orthogroup_label (str).
    """
    df = pd.DataFrame(counts)
    # Rename orthogroups
    df["Orthogroup"] = df["Orthogroup"].replace(rename)
    # Re-order
    df["Orthogroup_order"] = pd.Categorical(
        df["Orthogroup"],
        categories=order,
        ordered=True
    )
    # Sort DataFrame
    df.sort_values(["Orthogroup_order"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    df["Orthogroup_label"] = df["Orthogroup"].str.replace("(","\n(", regex= False)
    return df

# Bar Chart
def plot_counter(df, output_file= "orthogroups_counts.png"):
    """
    Cretes and saves the bar chart as a png of the sequence counts
    per orthogroup. 

    Args:
    df(pandas DataFrame): DataFrame output by function def custom.
    output_file (str): Creates a Figure with Matplotlib.
    
    Returns:
    Figure (png): Creates a Sequence Count Figure.
    """
    # Define color of the bar plot
    cmap = plt.get_cmap("Blues") 
    colors = cmap(np.linspace(0.4, 0.9, len(df)))
    plt.figure(figsize=(20,9))
    bars = plt.bar(df["Orthogroup_label"], df["SequenceCount"], color=colors)

    # Labels above bars
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            height + 0.2,
            str(int(height)),
            ha= 'center', fontsize=15
        )

    # Highlight first 3 bars
    for i, bar in enumerate(bars):
        if i < 3:
            bar.set_hatch("//")
            bar.set_edgecolor("black")
        else:
            bar.set_edgecolor("black")

    # Add legends
    p1 = mpatches.Patch(
        facecolor="lightgray", edgecolor="black", hatch="//",
        label="Pathway 1: Enzymes of the Phenylpropanoid Pathway"
    )
    p2 = mpatches.Patch(
        facecolor="lightgray", edgecolor="black",
        label="Pathway 2: (-)-Epicatechin Biosynthetic Enzymes"
    )
    lgnd = plt.legend(
        handles= [p1,p2], title="Legend", loc="upper right", fontsize=13
    )
    lgnd.get_title().set_fontsize(14)

    plt.xticks(rotation=0, fontsize=15)
    plt.xlabel("Orthogroups", fontsize=15)
    plt.ylabel("Number of Sequence", fontsize=16)
    plt.tight_layout()

    #Save png and show
    plt.savefig(output_file, dpi=300)
    plt.show()

# Main
def main():
    """
    Main Workflow: count sequences, prepares DataFrame  and plots
    """
    folder_path = "Orthogroups_FASTAS" # File name of your input
    rename = {
        "OG0000380": "OG0000380 (PAL)",
        "OG0001216": "OG0001216 (4CL)",
        "OG0002553": "OG0002553 (C4H)",
        "OG0000281": "OG0000281 (CHS)",
        "OG0005496": "OG0005496 (CHI)",
        "OG0007848": "OG0007848 (F3H)",
        "OG0001314": "OG0001314 (F3'H)",
        "OG0005192": "OG0005192 (DFR)",
        "OG0010355": "OG0010355 (ANS)",
        "OG0004494": "OG0004494 (ANR)"
    }
    order = [
        "OG0000380 (PAL)",
        "OG0001216 (4CL)",
        "OG0002553 (C4H)",
        "OG0000281 (CHS)",
        "OG0005496 (CHI)",
        "OG0007848 (F3H)",
        "OG0001314 (F3'H)",
        "OG0005192 (DFR)",
        "OG0010355 (ANS)",
        "OG0004494 (ANR)",
    ]

    # Run main
    counts = count_sequences(folder_path)
    df = custom(counts, rename, order)
    plot_counter(df, output_file="orthogroups_counts.png")

# Run script
if __name__ == "__main__":
    main()
