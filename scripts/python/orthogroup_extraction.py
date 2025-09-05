import pandas as pd
from pathlib import Path
from Bio import SeqIO

proteome_dir = "proteomes"
orthogroups_tsv = "Orthogroups.tsv"
species_code_excel = "species_code1.xlsx"
output_dir = "Orthogroups_FASTAS"

# Protein IDs of Arabidopsis thaliana (EC biosynthetic Pathway)
protein_ids = {
    'NP_181241.1', 'NP_176686.1', 'NP_180607.1',
    'NP_196897.1', 'NP_191072.1', 'NP_190692.1',
    'NP_196416.1', 'NP_199094.1', 'NP_194019.1', 'NP_176365.1'
}

# Load species dictionary
def load_species(path):
    """
    Load the species code from an Excel file. 

    Args: 
    path (str): Path to the excel file.

    Returns: 
    A dictionary mapping lowercase species name to species codes.

    Raises:
    FileNotFoundError: If the excel file does not exist.
    ValueError: If the excel files has fewer than the columns needed.
    """
    df = pd.read_excel(path)
    return {row[0].strip().lower(): row[1].strip() for _, row in df.iterrows()}

# Load sequences and assign species
def load_sequences(proteome_dir, species_codes):
   """
    Load the protein sequences from the proteome FASTA files. Extracts the sequence records using Biopython. 
    For each record infers the species scientific name and species code. Then the records is store in a dictionary.
    It creates a tuple for each ID: species code, description and orthogroup.

    Args: 
    proteomes_dir(str): Directory containing proteme FASTA files.
    species_code (dict): Mapping of lowercase species name to the species code.

    Returns: 
    A dictionary mapping gene_id (sequence, species_code, description)

    Raises:
    FileNotFoundError: If the proteomes file does not exist
    ValueError: If there a no matches for the species.
    """
   data = {}
   for file in Path(proteome_dir).glob("*.faa"):
        for record in SeqIO.parse(file, "fasta"):
            gene_id = record.id
            header = record.description.lower()
            species = next((species_codes[name] for name in species_codes if name in header), None)
            if not species:
                continue
            desc = record.description.replace(gene_id, "").split("[")[0].replace("hypothetical protein", "").strip()
            desc = desc.replace(" ", "_") if desc else "NA"
            data[gene_id] = (record.seq, species, desc)
        return data

# Find orthogroups that contain target proteins
def get_orthogroups(tsv, targets):
   """
    Identify orthogroups that contain at least one of the target IDs from Arabidopsis thaliana.
    Reads the TSV file and matchs the orthogroup ID with the target protein IDs.

    Args: 
    tsv(str): Path to Orthogroups.tsv file.
    targets (set): Set of protein IDs of interest.

    Returns: 
    set: set of Orhogroup IDs (str) that contain one or more target proteins.

    Raises:
    FileNotFoundError: If TSV file does not exist.
    ValueError: If  the file is empty.
    """
   df = pd.read_csv(tsv, sep="\t")
   return {
        row[0]
        for _, row in df.iterrows()
        if any(pid in str(cell) for cell in row[1:] if not pd.isna(cell) for pid in targets)
    }

# Write matching orthogroup sequences
def extract_orthogroups(tsv, sequences, matching_set, outdir):
  """
    Extract and saves the sequences for the matching orthogroups into FASTA files. 
    Iterates over Orthogroups.tsv fle and for each orthogroup matches the protein IDs.
    It retrieves the sequence information and from the sequences dictionary creates
    a new FASTA file with standarised headers.

    Args: 
    tsv(str): Path to Orthogroups.tsv file.
    sequence(dict):Mapping of the gene_ids
    outdir (set): Output directory of FASTA files.

    Returns: 
    None

    Raises:
    FileNotFoundError: If TSV file does not exist.
    ValueError: If the output directory can not be created. 
    """
  Path(outdir).mkdir(exist_ok=True)
  df = pd.read_csv(tsv, sep="\t")
  for _, row in df.iterrows():
        og = row[0]
        if og not in matching_set:
            continue
        records = []
        for cell in row[1:]:
            if pd.isna(cell): continue
            for gid in cell.split(", "):
                if gid in sequences:
                    seq, sp, desc = sequences[gid]
                    new_id = f"{sp}_{gid}_{desc}_{og}"
                    record = SeqIO.SeqRecord(seq, id=new_id, description="")
                    records.append(record)
        if records:
            SeqIO.write(records, Path(outdir) / f"{og}.faa", "fasta")

if __name__ == "__main__":
    # Load species code dictionary from Excel 
    species_codes = load_species(species_code_excel)

    # Load all proteome sequences and save them in memory
    sequences = load_sequences(proteome_dir, species_codes)

    # Identify Orthogroups with protein_ids
    matching_set = get_orthogroups(orthogroups_tsv, protein_ids)

    # Print the matching Orthogroups with the matching target IDs
    df = pd.read_csv(orthogroups_tsv, sep="\t")  
    print("Orthogroups containing targets:")
    for og in sorted(matching_set):
        row = df[df.iloc[:, 0] == og].iloc[0]

        hits = []
        # Scan all the genes IDs of the row
        for cell in row[1:]:
            if pd.isna(cell):  
                continue
            for gid in cell.split(", "):
                if gid in protein_ids:
                    hits.append(gid)

        # Print the orthogroup names and the target ids found
        print(f" - {og}: {', '.join(hits)}")

    # Extract all sequences that matches and save them into FASTA files
    extract_orthogroups(orthogroups_tsv, sequences, matching_set, output_dir)


