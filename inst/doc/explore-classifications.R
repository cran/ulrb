## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ulrb)
library(cluster)
library(dplyr)
library(ggplot2)
library(tidyr)
# a vector with some colors
qualitative_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## -----------------------------------------------------------------------------
# Load raw OTU table from N-ICE
load("../data/nice_raw.rda")

# Change name of first column
nice_clean <- rename(nice_raw, Taxonomy = "X.SampleID")

# Select 16S rRNA amplicon sequencing samples
selected_samples <- c("ERR2044662", "ERR2044663", "ERR2044664",
                      "ERR2044665", "ERR2044666", "ERR2044667",
                      "ERR2044668", "ERR2044669", "ERR2044670")

# Add a column with phylogenetic units ID (OTU in this case)
nice_clean <- mutate(nice_clean, OTU = paste0("OTU_", row_number()))

# Select relevant collumns
nice_clean <- select(nice_clean, selected_samples, OTU, Taxonomy)

# Separate Taxonomy column into each taxonomic level
nice_clean <- separate(nice_clean, 
                       Taxonomy,
                       c("Domain","Kingdom","Phylum",
                         "Class","Order","Family",
                         "Genus","Species"), 
                       sep=";")

# Remove Kingdom column, because it is not used for prokaryotes
nice_clean <- select(nice_clean, -Kingdom)

# Remove eukaryotes
nice_clean <- filter(nice_clean, Domain != "sk__Eukaryota")

# Remove unclassified OTUs at phylum level
nice_clean <- filter(nice_clean, !is.na(Phylum))

# Simplify name
nice <- nice_clean

# Tidy data
nice_tidy <- prepare_tidy_data(nice, 
                               sample_names = selected_samples, 
                               samples_in = "cols")

## -----------------------------------------------------------------------------
rb_default <- define_rb(nice_tidy)

## ---- fig.width=5, fig.height=4-----------------------------------------------
plot_ulrb_clustering(rb_default, 
                     sample_id = "ERR2044662",
                     taxa_col = "OTU",
                     log_scaled = TRUE)

## ---- fig.width=8, fig.height=8-----------------------------------------------
plot_ulrb_clustering(rb_default, 
                     taxa_col = "OTU",
                     log_scaled = TRUE,
                     plot_all = TRUE)

## -----------------------------------------------------------------------------
rb_k2 <- define_rb(nice_tidy, classification_vector = c("Rare", "Abundant"))

## ---- fig.width = 5, fig.height= 4--------------------------------------------
plot_ulrb_clustering(rb_k2,
                     taxa_col = "OTU",
                     plot_all = TRUE, 
                     log_scaled = TRUE, 
                     colors = c("#009E73", "#F0E442"))

## ---- fig.width=15------------------------------------------------------------
gridExtra::grid.arrange(
  plot_ulrb_clustering(rb_default, 
                     taxa_col = "OTU",
                     log_scaled = TRUE,
                     plot_all = TRUE),
  plot_ulrb_clustering(rb_k2,
                     taxa_col = "OTU",
                     plot_all = TRUE, 
                     log_scaled = TRUE, 
                     colors = c("#009E73", "#F0E442")),
  nrow = 1
)

## -----------------------------------------------------------------------------
#
rb_k4 <- define_rb(nice_tidy, 
                   classification_vector = c("very rare", "rare", "abundant", "very abundant"))
#

## ---- fig.width=15, fig.height=4----------------------------------------------
# One sample as example
plot_ulrb(rb_k4,
          sample_id = "ERR2044662", 
          taxa_col = "OTU",
          colors = c("#009E73", "#F0E442", "grey","#CC79A7"),
          log_scaled = TRUE)

## ---- fig.width=15, fig.height=8----------------------------------------------
# all samples
plot_ulrb(rb_k4,
          taxa_col = "OTU",
          colors = c("#009E73", "#F0E442", "grey","#CC79A7"),
          log_scaled = TRUE,
          plot_all = TRUE)

## -----------------------------------------------------------------------------
#
rb_k5 <- define_rb(nice_tidy, 
                   classification_vector = c("very rare", "rare", "undetermined", "abundant", "very abundant"))

## ---- fig.width=15, fig.height=4----------------------------------------------
# One sample as example
plot_ulrb(rb_k5,
          sample_id = "ERR2044662", 
          taxa_col = "OTU",
          colors = qualitative_colors[1:5],
          log_scaled = TRUE)

## ---- fig.width=15, fig.height=8----------------------------------------------
# All samples
plot_ulrb(rb_k5,
          taxa_col = "OTU",
          colors = qualitative_colors[1:5],
          log_scaled = TRUE,
          plot_all = TRUE)

## -----------------------------------------------------------------------------
#
rb_k1 <- define_rb(nice_tidy, classification_vector = c("rare"))

## ---- fig.width=5, fig.height=4-----------------------------------------------
plot_ulrb_clustering(rb_k1, 
                     taxa_col = "OTU", 
                     colors = "green4", 
                     plot_all = TRUE, 
                     log_scaled = TRUE)

## -----------------------------------------------------------------------------
rb_sample1 <- nice_tidy %>% filter(Sample == "ERR2044662")

# Calculate maximum k
max_k_sample1 <- rb_sample1 %>% pull(Abundance) %>% unique() %>% length()
#
max_k_sample1
# Improvise a classification vector for maximum k
# that is just any vector with the same length
classification_vector_max_k_sample1 <- seq_along(1:max_k_sample1)
#
rb_sample1_max_k <- 
  define_rb(rb_sample1,
            classification_vector = classification_vector_max_k_sample1)
#
rb_sample1_max_k %>% select(OTU, Classification, Abundance) %>% head(10)

## -----------------------------------------------------------------------------
suggest_k(nice_tidy, detailed = TRUE)  

## -----------------------------------------------------------------------------
suggest_k(nice_tidy, detailed = TRUE, range = 10:20)  

## -----------------------------------------------------------------------------
## One sample
# To get values
check_avgSil(nice_tidy, sample_id = selected_samples[1])

# To plot results
check_avgSil(nice_tidy, sample_id = selected_samples[1], with_plot = TRUE)

## -----------------------------------------------------------------------------
## Davie-Boulding index
# To get values
check_DB(nice_tidy, sample_id = selected_samples[1])

# To plot results
check_DB(nice_tidy, sample_id = selected_samples[1], with_plot = TRUE)

## Calinsky-Harabasz index
# To get values
check_CH(nice_tidy, sample_id = selected_samples[1])

# To plot results
check_CH(nice_tidy, sample_id = selected_samples[1], with_plot = TRUE)

## -----------------------------------------------------------------------------
evaluate_sample_k(nice_tidy, sample_id = selected_samples[1], with_plot = TRUE)

## -----------------------------------------------------------------------------
## To get values
evaluate_k(nice_tidy)

## To plot
evaluate_k(nice_tidy, with_plot = TRUE)

## -----------------------------------------------------------------------------
# default option with average Silhouette score
suggest_k(nice_tidy)

# best k for Davies-Bouldin
suggest_k(nice_tidy, index = "Davies-Bouldin")

# best k for Calinsky-Harabasz
suggest_k(nice_tidy, index = "Calinsky-Harabasz")

## -----------------------------------------------------------------------------
automatic_classification <- define_rb(nice_tidy, automatic = TRUE)

# Plot automatic result

plot_ulrb(automatic_classification, taxa_col = "OTU", plot_all = TRUE)

## -----------------------------------------------------------------------------
more_complex_automatic_classification <- define_rb(nice_tidy, 
                                                   automatic = TRUE,
                                                   index = "Calinsky-Harabasz",
                                                   range = 4:6)

## -----------------------------------------------------------------------------
# Plot automatic result
plot_ulrb(more_complex_automatic_classification, 
          plot_all = TRUE, 
          taxa_col = "OTU", 
          colors = qualitative_colors[1:5],
          log_scaled = TRUE)

## -----------------------------------------------------------------------------
# Start by deciding the maximum range across the entire dataset
max_k <- nice_tidy %>%
    filter(Abundance > 0, !is.na(Abundance)) %>%
    group_by(Sample) %>%
    summarise(topK = length(unique(Abundance))) %>%
    ungroup() %>%
    pull(topK) %>%
    min()
# print maximum number of clusters allowed for all samples in the N-ICE dataset
max_k

## -----------------------------------------------------------------------------
evaluate_k(nice_tidy, with_plot = TRUE, range = 2:max_k)

