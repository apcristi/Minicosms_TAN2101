# Pipeline for Microviz taxonomy visualization. ps_MC corresponds to the phyloseq object created in phyloseq analysis
library(stringr)
   library(ggplot2)
   library(dplyr)
   library(tidyr)
   library(tibble)
   library(readr)
   library(maps)
   library(glue)
library("Cairo")
   library(phyloseq)
library("readxl")
library("openxlsx")
library("gplots")
library(forcats)
library(mgcv)
library(microbiome)
library(microViz)
library(ggpubr)
 #if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager", type = "binary") # (not binary if you're on linux)
#BiocManager::install(c("stringr", "ggplot2", "dplyr", "tidyr", "tibble", "readr", "maps", "glue", "phyloseq", "readxl", "openxlsx", "gplots","forcats", "mgcv", "microbiome" ))



#devtools::install_github("david-barnett/microViz")
#devtools::install_github("david-barnett/microViz@0.7.2")


## visualizing composition - barplots
# Merge replicates
variable1 = as.character(get_variable(ps_MC, "Experiment"))
variable2 = as.character(get_variable(ps_MC, "Treatment"))
sample_data(ps_MC)$NewPastedVar <- mapply(paste0, variable1, variable2, 
    collapse = "_")
merge_samples(ps_MC, "NewPastedVar")

MC_merged <- merge_samples(ps_MC, "NewPastedVar", fun = mean)

# Plot merged replicated at Genus level
MC_merged %>% 
   comp_barplot(
    tax_level = "Genus",  n_taxa = 15, sample_order = "default",
    merge_other = FALSE +
   labs(x = NULL, y = NULL) +
  theme(
    axis.ticks.x = element_blank()
  ) + theme(text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1, )) 
