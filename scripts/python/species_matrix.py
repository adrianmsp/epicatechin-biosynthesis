# ------------------------
# Libraries
# ------------------------
import os
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

# ------------------------
# User settings 
# ------------------------
fasta_directory = "Orthogroups_BLASTP"  # FASTA files folder, change if needed
species_code = "species_code2.xlsx"     # Species Code excel
enzyme_file = "enzyme_names.xlsx"       # Enzyme file with Pathway order


# Specific clade order 
clade_order = ["Bryophytes", "Lycophytes", "Ferns", "Gymnosperms", 
               "Basal Angiosperms", "Monocots", "Eudicots"]

# ------------------------
# Load species metdata (excel with 4-letter code)
# ------------------------
data_df = pd.read_excel(species_code)
data_df["Code_lower"] = data_df["Code"].str.lower()

# Dictionaries: code to matrix maps 4-letter code to species scientific names and
# matrix clade maps scientific names to clades orders
code_matrix = dict(zip(data_df["Code_lower"], data_df["Matrix Name"]))
matrix__clade =  dict(zip(data_df["Matrix Name"], data_df["Clade"]))

# ------------------------
# Load enzyme excel for EC pathway order
# ------------------------
enzyme_df = pd.read_excel(enzyme_file)
assert "Orthogroup Name" in enzyme_df.columns and "Protein Name" in enzyme_df.columns
    

# Orthogroup ID (folder/file name)
orthogroup_rename = dict(zip(enzyme_df["Orthogroup Name"], enzyme_df["Protein Name"]))

# Sort enzymes by EC pathway 
if "Pathway Step" in enzyme_df.columns: # Checks for Column in excel Pathway Step
    enzyme_df = enzyme_df.sort_values("Pathway Step")
enzyme_order = enzyme_df["Protein Name"].tolist()  # Desire enzyme order

# ------------------------
# Counts proteins per species of every orthogroup FASTA files
# ------------------------
counts = defaultdict(lambda: defaultdict(int)) 
for file in os.listdir(fasta_directory):
    orthogroup = file.replace(".faa","")
    path = os.path.join(fasta_directory, file)

    # Record.id starts with 4-letter species code
    for record in SeqIO.parse(path, "fasta"):
        code = record.id[:4].lower()
        species_name = code_matrix.get(code)
        if species_name:
            counts[orthogroup][species_name] += 1

# ------------------------
# Species-orthogroup matrix
# ------------------------
species_df = pd.DataFrame(counts). fillna(0).astype(int)
# Renames orthogroup columns with enzyme names
species_df.rename(columns=orthogroup_rename, inplace=True)
# Columns appear in the desired order of the pathway
species_df = species_df.reindex(columns=enzyme_order, fill_value=0)
# Attach clade to species name (row index is species-matrix name)
species_df.index.name = species_df.index.name or "Species"
species_df["Clade"] = species_df.index.map(matrix__clade).fillna("Uknown")

# Make Clade an ordered categorical 
species_df["Clade"] = pd.Categorical(species_df["Clade"], categories=clade_order, ordered=True)

# Sort by Clade
species_df = (
    species_df
    .assign(Species= species_df.index.astype(str))
    .sort_values(by=["Clade"], kind="stable")
)

# Keep a list with the clade strings aligned with row order 
row_clades = species_df["Clade"].astype(str).tolist()


# ------------------------
# Heatmap: presence/absence per species 
# ------------------------
absence = "#f2f2f2"
presence = "#1f77b4"
cmap = ListedColormap([absence, presence])
norm = BoundaryNorm([-0.55, 0.5, 1.5], ncolors=2)

# Convert counts to 0/1
heatmap_df = (species_df[enzyme_order] > 0).astype(int)

# Color palette for clade labels in the y-axis
clade_order = ["Bryophytes","Lycophytes","Ferns","Gymnosperms","Basal Angiosperms","Monocots","Eudicots"]
clade_palette = {
    "Bryophytes": "#FFECB3",       
    "Ferns": "#C8E6C9",             
    "Lycophytes": "#A5D6A7",        
    "Gymnosperms": "#BBDEFB",       
    "Basal Angiosperms": "#F8BBD0", 
    "Monocots": "#D1C4E9",          
    "Eudicots": "#B2FBFF"           
}

# Build the figure
fig = plt.figure(figsize=(12, 11), constrained_layout=True)
gs = fig.add_gridspec(ncols=1)
ax_heat = fig.add_subplot(gs[0])

# Heatmap with seaborn
sns.heatmap(
    heatmap_df,
    cmap=cmap,
    norm=norm,
    linewidths=0.5,
    linecolor="lightgray",
    cbar=False,
    ax=ax_heat
)

# Color the background of each tick label with its clade color
for label, clade in zip(ax_heat.get_yticklabels(), row_clades):
    label.set_bbox(dict(facecolor = clade_palette.get(clade,  "#999999"),
                        edgecolor="none", boxstyle= "round,pad=0.15")) 
    
# Axis labels names and style
ax_heat.set_xlabel("Enzymes", labelpad=10)
ax_heat.set_ylabel("Species")
plt.setp(ax_heat.get_xticklabels(), rotation=0)
plt.setp(ax_heat.get_yticklabels(), rotation=0)

# Legends: Clade colors and presence/absence 
clades_handles = [
    Patch(facecolor=clade_palette[c], label=c)
    for c in clade_order if c in row_clades
]
presence_handles = [
    Patch(facecolor=presence, label="Present (1)"),
    Patch(facecolor=absence, label="Absent (0)")
]
# Position of the lengends
legend1= ax_heat.legend(handles= clades_handles, title="Clades", 
                        loc="upper left", prop={'size': 8}, title_fontsize=9, bbox_to_anchor= (1.02,1))

ax_heat.add_artist(legend1)
ax_heat.legend(handles= presence_handles, title="Presence", 
                        loc="upper left", prop={'size': 8}, title_fontsize=9, bbox_to_anchor= (1,0.82))

# save heatmap figues as a png
plt.tight_layout()
plt.savefig("orthogroup_species.png", dpi=300)

# ------------------------
# Matrix per Clade
# ------------------------
enzyme_columns = [e for e in enzyme_order if e in species_df.columns]

# Collapse species to clades 1 if any species in the clade has the enzyme
clade_matrix = (
    species_df.groupby("Clade", observed= True)[enzyme_columns]
            .max()
            .astype(int)
 )
# Save clade matrix as an excel
clade_matrix.reset_index().to_excel("clade_matrix.xlsx", index=False)
