## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ulrb)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(purrr)
#
set.seed(123)

## -----------------------------------------------------------------------------
# Load raw OTU table from N-ICE
load("../data/nice_raw.rda")
load("../data/nice_env.rda")

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
nice_clean <- separate(nice_clean, Taxonomy,
                       c("Domain","Kingdom","Phylum",
                         "Class","Order","Family","Genus","Species"),
                       sep=";")

# Remove Kingdom column, because it is not used for prokaryotes
nice_clean <- select(nice_clean, -Kingdom)

# Remove eukaryotes
nice_clean <- filter(nice_clean, Domain != "sk__Eukaryota")

# Remove unclassified OTUs at phylum level
nice_clean <- filter(nice_clean, !is.na(Phylum))

# Simplify name
nice <- nice_clean

# Check if everything looks normal
head(nice)

# Change table to tidy format
# You can automatically do this with an ulrb function
nice_tidy <- prepare_tidy_data(nice, ## data to tidy 
                  sample_names = contains("ERR"), ## vector with ID samples
                  samples_in = "cols") ## samples can be in columns (cols) or rows. 


## ----fig.width = 6, fig.height = 4.5------------------------------------------
## First, check how many reads each sample got
nice_tidy %>% 
  group_by(Sample) %>% ## because data is in tidy format
  summarise(TotalReads = sum(Abundance)) %>% 
  ggplot(aes(Sample, TotalReads)) + 
  geom_hline(yintercept = 40000) + 
  geom_col() + 
  theme_bw() + 
  coord_flip()

## -----------------------------------------------------------------------------
nice_rarefid <- 
  nice_tidy %>% 
  group_by(Sample) %>%
  nest() %>% 
  mutate(Rarefied_reads = map(.x = data, 
                              ~as.data.frame(
                                t(
                                  rrarefy(
                                    .x$Abundance, 
                                    sample = 40000))))) %>% 
  unnest(c(data,Rarefied_reads))%>% 
  rename(Rarefied_abundance = "V1")  

## -----------------------------------------------------------------------------
nice_classified <- define_rb(nice_tidy, simplified = TRUE)
head(nice_classified)

## -----------------------------------------------------------------------------
nice_eco <- nice_classified %>% left_join(nice_env, by = c("Sample" = "ENA_ID"))

## -----------------------------------------------------------------------------
# Calculate diversity
nice_diversity <- nice_eco %>% 
  group_by(Sample, Classification, Month, Depth, Water.mass) %>% 
  summarise(SpeciesRichness = specnumber(Abundance),
            ShannonIndex = diversity(Abundance, index = "shannon"),
            SimpsonIndex = diversity(Abundance, index = "simpson")) %>% 
  ungroup()

## -----------------------------------------------------------------------------
## Re-organize table to have all diversity indices in a single col.
nice_diversity_tidy <- nice_diversity %>% 
  pivot_longer(cols = c("SpeciesRichness", "ShannonIndex", "SimpsonIndex"),
               names_to = "Index",
               values_to = "Diversity") %>% 
# correct order
  mutate(Month = factor(Month, levels = c("March", "April", "June"))) %>% 
# edit diversity index names
  mutate(Index = case_when(Index == "ShannonIndex" ~ "Shannon Index",
                           Index == "SimpsonIndex" ~ "Simpson Index",
                           TRUE ~ "Number of OTUs")) %>%
  # order diversity index names
  mutate(Index = factor(Index, c("Number of OTUs", "Shannon Index","Simpson Index")))

## ----fig.width = 6, fig.height = 4.5------------------------------------------
qualitative_colors <- 
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Plots
nice_diversity_tidy %>%
  ggplot(aes(Month, Diversity, col = Classification)) + 
  geom_point() + 
  facet_wrap(~Index, scales = "free_y")+
  scale_color_manual(values = qualitative_colors[c(5,3,6)]) + 
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "top") + 
  labs(col = "Classification: ",
       title = "Alpha diversity across month")

## ----fig.width = 6, fig.height = 4.5------------------------------------------
nice_diversity_tidy %>% 
  ggplot(aes(Water.mass, Diversity, col = Classification)) + 
  geom_point() + 
  facet_wrap(~Index, scales = "free_y")+
  scale_color_manual(values = qualitative_colors[c(5,3,6)]) + 
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "top") %>%
  labs(col = "Classification: ",
       title = "Alpha diversity across water mass")

## -----------------------------------------------------------------------------
rare_biosphere <- 
  nice_eco %>%  
  filter(Classification == "Rare") %>% 
  select(Sample, Abundance, OTU, Depth, Month, Water.mass) %>% 
  mutate(Sample_unique = paste(Sample, Month))

rb_env <- rare_biosphere %>% 
  ungroup() %>% 
  select(Sample_unique, Depth, Water.mass, Month) %>% 
  distinct()

rb_sp_prep <- rare_biosphere %>% 
  ungroup() %>% 
  select(Sample_unique, Abundance, OTU)

rb_sp_prep %>% head() # sanity check

rb_sp <- 
  rb_sp_prep %>% 
  tidyr::pivot_wider(names_from = "OTU",
                     values_from = "Abundance",
                     values_fn = list(count= list)) %>% 
  print() %>% 
  unchop(everything())

## -----------------------------------------------------------------------------
rb_sp[is.na(rb_sp)] <- 0

# Prepare aesthetics
rb_env <- rb_env %>% 
  mutate(col_month = case_when(Month == "March" ~ qualitative_colors[1],
                               Month == "April" ~ qualitative_colors[2],
                               TRUE ~ qualitative_colors[3])) %>% 
  mutate()

## ----fig.height = 8, fig.width= 8---------------------------------------------
cca_plot_rare_biosphere <- 
  cca(rb_sp[,-1], display = "sites", scale = TRUE) %>% 
  plot(display = "sites", type = "p", main = "Rare biosphere")
#
cca_plot_rare_biosphere
points(cca_plot_rare_biosphere,
       bg = rb_env$col_month,
       pch = 21, 
       col = "grey", 
       cex = 2)
#
with(rb_env,
     ordispider(cca_plot_rare_biosphere,
                Month, lty = "dashed", label = TRUE))

