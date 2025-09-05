import os
from collections import defaultdict
import subprocess

# Input and output folders
# cdhit_results should contain BLAST.tsv files and orthogroups.faa after CD-HIT.
# output_best stores the best-hit FASTA sequence per species.
input = "cdhit_results"
output = "output_best"
os.makedirs(output, exist_ok= True) # Creates output directory

# Finds all BLAST results if ending with .tsv from the input folder.
blast_files = [f for f in os.listdir(input) if f.endswith(".tsv")]

# Loops through the BLAST results file
for blast_file in blast_files:
    # Extracts the base name and removes the _blast.tsv name
    base = blast_file.replace("_blast.tsv", "")
    # File paths for BLAST, CD-HIT FASTA, ID list and final output
    blast_path = os.path.join(input, blast_file)
    fasta_path = os.path.join(input, f"{base}_cdhit.faa") # extracts sequence from file cdhit.faa
    best_ids = os.path.join(output, f"{base}_best.txt")   # stores sequence IDs.
    output_fasta = os.path.join(output, f"{base}.faa")   # final output FASTA file per-orthogroup
   
    # Troubleshooting ensures that FASTA files cdhit.faa exist
    if not os.path.exists(fasta_path):
        print(f"Missing FASTA file")
        continue
    # Parse BALST files for best hit per species
    # Creates a dictionary to store best hits per species per orthogroup.
    best_hit = {}
    with open(blast_path) as f:
        for line in f:
            cols = line.strip().split("\t")
            subject_id = cols[1] # Sequence ID
            bitscore = float(cols[-1]) # Bitscore
            species = subject_id.split("_")[0]
            
            # Keep the highest score per specie
            if species not in best_hit or bitscore > best_hit[species][1]:
                best_hit[species] = (subject_id, bitscore)

    # Parse BLAST output file
    # Write subject IDs to temp file, tab separted columns
    # col[1] = sequence ID
    # cols[-1] = bitscore
    with open(best_ids, "w") as out:
        for subject_id, _ in best_hit.values():
            out.write(subject_id + "\n" )
            
    # Seqtk to change txt to faa files
    subprocess.run([
        "seqtk", "subseq", fasta_path, best_ids
        ], stdout = open(output_fasta, "w"))
  
# Notify when processing is finalised
print("All fasta files done", output_fasta)
