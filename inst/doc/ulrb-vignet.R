## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ulrb)
library(dplyr)
library(tidyr)
library(ggplot2)

## -----------------------------------------------------------------------------
# Load raw OTU table from N-ICE
data("nice_raw", package = "ulrb")

# Change name of first column
nice_clean <- rename(nice_raw, Taxonomy = "X.SampleID")

# Select 16S rRNA amplicon sequencing samples
selected_samples <- c("ERR2044662", "ERR2044663", "ERR2044664","ERR2044665", "ERR2044666", "ERR2044667","ERR2044668", "ERR2044669", "ERR2044670")

# Add a column with phylogenetic units ID (OTU in this case)
nice_clean <- mutate(nice_clean, OTU = paste0("OTU_", row_number()))

# Select relevant columns
nice_clean <- select(nice_clean, selected_samples, OTU, Taxonomy)

# Separate Taxonomy column into each taxonomic level
nice_clean <- separate(nice_clean, Taxonomy, c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep=";")

# Remove Kingdom column, because it is not used for prokaryotes
nice_clean <- select(nice_clean, -Kingdom)

# Remove eukaryotes
nice_clean <- filter(nice_clean, Domain != "sk__Eukaryota")

# Remove unclassified OTUs at phylum level
nice_clean <- filter(nice_clean, !is.na(Phylum))

# Simplify name
nice <- nice_clean

# Quick look at the table
head(nice)

## -----------------------------------------------------------------------------
nice_tidy <- prepare_tidy_data(nice, sample_names = selected_samples, samples_in = "cols")

## -----------------------------------------------------------------------------
classified_table <- define_rb(nice_tidy)

# Quick output check
colnames(classified_table)

classified_table %>% 
  select(OTU, Sample, Abundance, 
         Classification, Silhouette_scores, Cluster_median_abundance, 
         pam_object) %>% 
  head()

## -----------------------------------------------------------------------------
# Simple automation example
define_rb(nice_tidy, automatic =  TRUE)

## ----fig.height = 6, fig.width = 6--------------------------------------------
# One sample as example
plot_ulrb_clustering(classified_table, 
                       sample_id = selected_samples[1],
                       taxa_col = "OTU") +
  labs(title = paste("Clustering for sample", selected_samples[1]))

# All samples, with centrality metric
plot_ulrb_clustering(classified_table,
                     taxa_col = "OTU", 
                     plot_all = TRUE, 
                     log_scaled = TRUE) +
  labs(title = "Clustering for all samples")

## ----fig.width = 10-----------------------------------------------------------
# One sample as example
plot_ulrb_silhouette(classified_table,
                     sample_id = selected_samples[1],
                     taxa_col = "OTU") +
  labs(title = paste("Silhouette plot of sample", selected_samples[1]))

#
plot_ulrb_silhouette(classified_table,
                     sample_id = selected_samples[1],
                     taxa_col = "OTU",
                     plot_all = TRUE) +
  labs(title = "Silhouette plot of all samples")


## ----fig.width = 10-----------------------------------------------------------
# For a single sample
plot_ulrb(classified_table,
          sample_id =  selected_samples[1],
          taxa_col = "OTU")

# For all samples
plot_ulrb(classified_table,
          taxa_col = "OTU",
          plot_all = TRUE) 

