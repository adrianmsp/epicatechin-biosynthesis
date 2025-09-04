```{r}
#------------------------------
# Clean Environment
#------------------------------
rm(list = ls(all.names=TRUE))  # clean environment

#------------------------------
# Load libraries
#------------------------------
library(tidyverse)
library(pheatmap)
library(grid)

#------------------------------
# Load Gene Expression Analysis from Expression Atlas: E-GEOD-38612
#------------------------------
atlas <- read.delim("~/E-GEOD-38612-query-results.tpms.tsv", comment.char = "#", header = TRUE, check.names = FALSE) 
# Change folder path
# Make data frame with my Expression Atlas: E-GEOD-38612 data, lines that start with # are comments and ignored, first line is header and check.names keeps column names as they are.
head(atlas) # show head of atlas

#------------------------------
# (-)-Epicatechin Biosynthesis genes of interest (TAIR IDs)
#------------------------------
epic_genes <- c("AT2G37040", "AT1G65060", "AT2G30490", "AT5G13930", "AT3G55120", "AT3G51240", "AT5G07990", "AT5G42800", "AT4G22880", "AT1G61720") # Make Dataframe with my TAIR IDs of interest

#------------------------------
# (-)-Epicatechin Biosynthesis genes Expression Atlas 
#------------------------------
epic_atlas <- atlas %>%                    # Create data frame with TAIR IDs and leaf, flower, fruit and root data
  filter(`Gene ID` %in% epic_genes) %>%    # Filter the expression atlas with my TAIR IDs of interest
  column_to_rownames(var = "Gene ID") %>%  # Set Gene ID column to row names
  select(leaf, flower, fruit, root)        # Select leaf, flower, fruit and root data 

#------------------------------
# Rename columns
#------------------------------
colnames(epic_atlas) <- c("Leaf", "Flower", "Fruit", "Root") # Change names "capitalise words"

#------------------------------
# Data Normalization
#------------------------------
log_atlas <- log2(epic_atlas +1 ) # log2 normalization to reduce skew

#------------------------------
# Order Genes EC Pathway
#------------------------------
epic_order <- intersect(epic_genes, rownames(log_atlas)) # Order of EC Pathway, as the genes Dataframe epic_genes
log_atlas <- log_atlas[epic_order, , drop = FALSE]

#------------------------------
# Plot Pheatmap
#------------------------------
heat_plot <- pheatmap(log_atlas, scale = "row", # expression scaled by row 
         display_numbers = TRUE, # Display numbers of expression
         number_color = "black", # Color numbers black
         fontsize_number = 10,   # Fontsize numbers
         fontsize_row = 10,      # Row label font size
         fontsize_col = 10,      # Column label font size
         cluster_rows = FALSE,   # Not show dendogram for rows
         cluster_cols = FALSE,   # Not show dendogram for columns
         angle_col = 0,          # angle for column labels
         )
#------------------------------
# Save Plot
#------------------------------
ggsave(filename = "GeneExpression.svg", plot = heat_plot,    # Save plot as svg file
       width = 10, height = 7, dpi = 300)
show(heat_plot)                                              # Show plot
```
