---
title: "Minicosms_metabarcoding"
author: "Antonia"
date: "02 Feb 2025"
output:
  pdf_document: default
  html_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r}
#load packages
library("phyloseq")
library("ggplot2")
library("dplyr")
library("tidyr")
library("tibble")
library("readxl")
library("readr")
library("stringr")
library("rmarkdown")
library("yaml")
library(tidyverse)
library(viridis)
library(forcats)
library(ggpubr)
library(vegan)
library(pairwiseAdonis)
library("Matrix")
library("reshape2")
```


```{r}
#specify color pallete
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 16
mycolors <- colorRampPalette(brewer.pal(16, "Set1"))(nb.cols)

# Create a ggplot with 18 colors 
# Use scale_fill_manual

#Set theme
theme_set(theme_bw())
```


```{r}
#Create phyloseq object
otu_mat<- read_xlsx("seqtab_nonchim_18S.xlsx")
tax_mat<- read_xlsx("taxa.xlsx")
samples_df <- read_excel("Metadata.xlsx")
```

```{r}
otu_mat <- otu_mat %>%
    tibble::column_to_rownames("OTUNumber") 

tax_mat <- tax_mat %>% 
    tibble::column_to_rownames("OTUNumber")

  samples_df <- samples_df %>% 
    tibble::column_to_rownames("Sample_ID") 
```

```{r}
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
```

```{r}
OTU = otu_table(otu_mat , taxa_are_rows = TRUE)
 TAX = tax_table(tax_mat)
 samples = sample_data(samples_df)
ps_MC<- phyloseq(OTU, TAX, samples)
ps_MC
```


```{r}
#leave only MC samples
ps_MC<- subset_samples(ps_MC, Minicosms=="yes")
```

```{r}

#filter data to keep only photosnthetic phylas
ps_MC<- subset_taxa(ps_MC, Division %in% c("Chlorophyta", "Dinoflagellata", "Cryptophyta", 
                                                "Haptophyta", "Ochrophyta", "Cercozoa", "Radiolaria"))
ps_MC <- subset_taxa(ps_MC, !(Class %in% c("Syndiniales", "Sarcomonadea")))
ps_MC
```

```{r}
#Normalize number of reads in each sample using median sequencing depth.
 total = median(sample_sums(ps_MC))
 standf = function(x, t=total) round(t * (x / sum(x)))
ps_MC = transform_sample_counts(ps_MC, standf)
```

## Subset by experiment
```{r}
MC_1<- subset_samples(ps_MC, Experiment=="E1")
MC_2<- subset_samples(ps_MC, Experiment=="E2")
MC_3<- subset_samples(ps_MC, Experiment=="E3")
```

## T0 taxonomy

```{r}
T0 <- subset_samples(ps_MC, Replicate=="T0")

ps_MC_abund <- filter_taxa(T0, function(x) sum(x > total*0.02) > 0, TRUE)

RelAb = transform_sample_counts(ps_MC_abund , function(x) x / sum(x) )

p <- plot_bar(RelAb, x="Experiment")

p + geom_bar(aes(fill=Genus), stat="identity", position="stack") + theme(text = element_text(size=20),axis.text.x = element_text(angle=45, hjust=1)) + scale_fill_manual(values = mycolors)  + ylab("18S Reads abundance")  #+coord_flip()

#ggsave(filename = "Exp_Class_reab.pdf", width = 9, height = 3, dpi= 300)
```

##Taxonomy per pexperiment

#Merge replicates
```{r}
variable1 = as.character(get_variable(ps_MC, "Experiment"))
variable2 = as.character(get_variable(ps_MC, "Treatment"))
sample_data(ps_MC)$NewPastedVar <- mapply(paste0, variable1, variable2, 
    collapse = "_")
merge_samples(ps_MC, "NewPastedVar")

MC_merged <- merge_samples(ps_MC, "NewPastedVar", fun = mean)
```


```{r}

ps_MC_abund <- filter_taxa(MC_merged, function(x) sum(x > total*0.12) > 0, TRUE)

RelAb = transform_sample_counts(ps_MC_abund, function(x) x / sum(x) )

plot_bar(RelAb, fill="Genus") + geom_bar(aes(fill=Class), stat= "identity", position="stack") + theme(text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1)) +  scale_fill_manual(values = mycolors) + coord_flip() #
 ggsave(filename = "MC_relab.png", width = 40, height = 7, dpi= 300, type = "cairo")

```


## Diversity
```{r}
##Alpha diversity

Richness <- estimate_richness(ps_MC, split = TRUE, measures=c("Shannon", "Chao1", "observed"))
write.csv(Richness, file = "Richness.csv")

# Diveristy was added to the metadata
Richeness_edited <- read_xlsx("Metadata.xlsx")
```

#Chao
# M1
```{r}
M1_MT <- filter(Richeness_edited, Experiment=="E1")

M1_diversity_Chao <- M1_MT %>% ggplot( aes(x=Number_2, y=Chao1, fill=Treatment)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.5, alpha=0.9) +
   ggtitle("MC1") + ylim(40, 100) + theme(text = element_text(size = 15))

M1_diversity_Chao <- M1_diversity_Chao + scale_fill_manual(values =c("#ffc100", "#00688b", "#e0301e","#602320","#999999"))
```

```{r}
M1_diversity_Chao
```


# M2
```{r}
M2_MT <- filter(Richeness_edited, Experiment=="E2") 


M2_diversity_Chao <- M2_MT %>% ggplot( aes(x=Number_2, y=Chao1, fill=Treatment)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    ggtitle("MC2")  + ylim(40, 100) + theme(text = element_text(size = 15))

M2_diversity_Chao <- M2_diversity_Chao  + scale_fill_manual(values =c("#ffc100", "#00688b", "#e0301e","#602320","#999999"))
```


```{r}
M2_diversity_Chao
```


# M3
```{r}
M3_MT <- filter(Richeness_edited, Experiment=="E3") 


M3_diversity_Chao <- M3_MT %>% ggplot( aes(x=Number_2, y=Chao1, fill=Treatment)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    ggtitle("MC3")  + ylim(40, 100)  + theme(text = element_text(size = 15))

  M3_diversity_Chao <- M3_diversity_Chao  + scale_fill_manual(values =c("#ffc100", "#00688b", "#e0301e","#602320","#999999"))
```


```{r}
 M3_diversity_Chao
```


#Shannon
#MC1
```{r}
M1_diversity_shannon <- M1_MT %>% ggplot( aes(x=Number_2, y=Shannon, fill=Treatment)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.5, alpha=0.9) +
   ggtitle("MC1") + ylim(1.5, 3) + theme(text = element_text(size = 15))

M1_diversity_shannon <- M1_diversity_shannon + scale_fill_manual(values =c("#ffc100", "#00688b", "#e0301e","#602320","#999999"))
```


```{r}
M1_diversity_shannon
```


```{r}
M2_MT <- filter(Richeness_edited, Experiment=="E2") 


M2_diversity_shannon <- M2_MT %>% ggplot( aes(x=Number_2, y=Shannon, fill=Treatment)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    ggtitle("MC2") + ylim(1.5, 3) + theme(text = element_text(size = 15))

M2_diversity_shannon <- M2_diversity_shannon  + scale_fill_manual(values =c("#ffc100", "#00688b", "#e0301e","#602320","#999999"))
```


```{r}
M2_diversity_shannon
```


# M3
```{r}
M3_MT <- filter(Richeness_edited, Experiment=="E3") 


M3_diversity_shannon <- M3_MT %>% ggplot( aes(x=Number_2, y=Shannon, fill=Treatment)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    ggtitle("MC3") + ylim(1.5, 3)   + theme(text = element_text(size = 15))

M3_diversity_shannon <- M3_diversity_shannon  + scale_fill_manual(values =c("#ffc100", "#00688b", "#e0301e","#602320","#999999"))
```


```{r}
M3_diversity_shannon
```


```{r}
ggarrange(M1_diversity_Chao, M2_diversity_Chao, M3_diversity_Chao, M1_diversity_shannon, M2_diversity_shannon, M3_diversity_shannon,
          labels = c("A", "B", "C", "D", "E", "F"), 
          ncol = 3, nrow = 2) 

ggsave("MC_Diveristy_Shannon.pdf",  dpi = 300,  width =10, height =5, scale=1.2)
```


## B-diveristy (NMDS)

## M1 
```{r}
MC_1_nmds <- subset_samples(MC_1, Replicate != "T0" )

ps_dada2.ord <- ordinate(MC_1_nmds, "NMDS", "bray")
```


```{r}
MC1 <- plot_ordination(MC_1_nmds, ps_dada2.ord, type="samples", color="Treatment", title="MC1",  shape="Replicate") + geom_point(size=3) + scale_color_manual(values =c("#ffc100", "#00688b","#602320", "#e0301e", "#999999"))  +xlim(-0.4,0.4) + ylim(-0.3,0.3)

ggsave(filename = "MC1_NMDs.png", width = 4, height = 3 , type="cairo")
```

## MC2
```{r}
MC_2_nmds <- subset_samples(MC_2, Replicate != "T0" )

ps_dada2.ord <- ordinate(MC_2_nmds, "NMDS", "bray")
```


```{r}
MC2 <- plot_ordination(MC_2_nmds, ps_dada2.ord, type="samples", color="Treatment", title="MC2",  shape="Replicate") + geom_point(size=3) + scale_color_manual(values =c("#ffc100", "#00688b","#602320", "#e0301e", "#999999")) +xlim(-0.4,0.4) + ylim(-0.3,0.3)

ggsave(filename = "MC2_NMDs.png", width = 4, height = 3 , type="cairo")
```

## MC3
```{r}
MC_3_nmds <- subset_samples(MC_3, Replicate != "T0" )

ps_dada2.ord <- ordinate(MC_3_nmds, "NMDS", "bray")
```


```{r}
MC3 <- plot_ordination(MC_3, ps_dada2.ord, type="samples", color="Treatment", title="MC3",  shape="Replicate") + geom_point(size=3) + scale_color_manual(values =c("#ffc100", "#00688b","#602320", "#e0301e", "#999999")) +xlim(-0.4,0.4) + ylim(-0.3,0.3)

ggsave(filename = "MC3_NMDs.png", width = 4, height = 3 , type="cairo")
```


```{r}
ggarrange(MC1, MC2, MC3,labels = c("A", "B", "C"), 
          ncol = 3, nrow = 1) 

ggsave("MC_Bdiveristy.pdf",  dpi = 300,  width =10, height =3, scale=1.2)
```

#Heatmap species abundance
```{r}
#Heatmap
ps_MC_abund <- filter_taxa(ps_MC, function(x) sum(x > total*0.03) > 0, TRUE)
ps_MC_abund

plot_heatmap(ps_MC_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", taxa.label = "Species",  taxa.order = "Class", sample.order = "Number", trans=NULL, low="white", high="blue", na.value="white") + theme(text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1.0, vjust = 1.0, size = 12)) 

#ggsave(filename = "Heatmap_MC.pdf", width = 15, height = 7 )
```

##DESeq2 
```{r}
library("DESeq2")
```

## MC1 as an example.  
```{r}
## Filter all experiments for A and D

MC1_AD <- MC_1%>%  subset_samples(Treatment %in% c("A", "B"))

# design formula
diagdds = phyloseq_to_deseq2(MC1_AD, ~ Treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res <- results(diagdds)
summary(res)
```


```{r}
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(MC1_AD)[rownames(sigtab), ], "matrix"))
head(sigtab)
```


```{r}
write.csv(sigtab, file = "Deseq_MC1_AD.csv")
dim(sigtab)
```


```{r}
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Class), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))


ggplot(sigtab, aes(x = fct_reorder(Genus, log2FoldChange, .desc = TRUE), color=Class, y = log2FoldChange)) +
  geom_point(size=3) + theme(text = element_text(size=10), axis.text.x = element_text(angle = 45, hjust = 1.0, vjust=1.0)) +
  coord_flip() + scale_colour_manual (values = mycolors)
```


```{r}
## Filter all experiments for A and C

MC1_AC <- MC_1 %>%  subset_samples(Treatment %in% c("A", "C"))

# design formula
diagdds = phyloseq_to_deseq2(MC1_AC, ~ Treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res <- results(diagdds)
summary(res)
```


```{r}
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(MC1_AC)[rownames(sigtab), ], "matrix"))
head(sigtab)
```


```{r}
write.csv(sigtab, file = "Deseq_MC1_AC.csv")
dim(sigtab)
```


```{r}
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Class), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))


ggplot(sigtab, aes(x = fct_reorder(Genus, log2FoldChange, .desc = TRUE), color=Class, y = log2FoldChange)) +
  geom_point(size=3) + theme(text = element_text(size=10), axis.text.x = element_text(angle = 45, hjust = 1.0, vjust=1.0)) +
  coord_flip() + scale_colour_manual (values = mycolors)
```


```{r}
## Filter all experiments for A and B
MC1_AB <- MC_1 %>%  subset_samples(Treatment %in% c("A", "B"))

# design formula
diagdds = phyloseq_to_deseq2(MC1_AB, ~ Treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res <- results(diagdds)
summary(res)
```


```{r}
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(MC1_AB)[rownames(sigtab), ], "matrix"))
head(sigtab)
```


```{r}
write.csv(sigtab, file = "Deseq_MC1_AB.csv")
dim(sigtab)
```


```{r}
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Class), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))


ggplot(sigtab, aes(x = fct_reorder(Genus, log2FoldChange, .desc = TRUE), color=Class, y = log2FoldChange)) +
  geom_point(size=3) + theme(text = element_text(size=10), axis.text.x = element_text(angle = 45, hjust = 1.0, vjust=1.0)) +
  coord_flip() + scale_colour_manual (values = mycolors)
```