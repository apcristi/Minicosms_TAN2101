---
title: "Minicosms_metabarcoding"
author: "Antonia"
date: "15th Jan 2026"
output: html_document
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
library(ggpubr)
library(tidyverse)
library(viridis)
library(forcats)
library(vegan)
library(tidyverse)
library("ggrepel")
library("Matrix")
library("reshape2")
library("ALDEx2")

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
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
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

#Filter by MC 

```{r}
MC_1<- subset_samples(ps_MC, Experiment=="E1")
MC_2<- subset_samples(ps_MC, Experiment=="E2")
MC_3<- subset_samples(ps_MC, Experiment=="E3")
```


#Merge replicates

```{r}
variable1 = as.character(get_variable(ps_MC, "Experiment"))
variable2 = as.character(get_variable(ps_MC, "Treatment"))
sample_data(ps_MC)$NewPastedVar <- mapply(paste0, variable1, variable2, 
    collapse = "_")
merge_samples(ps_MC, "NewPastedVar")

MC_merged <- merge_samples(ps_MC, "NewPastedVar", fun = mean)
```

## ALDEx2 for MC1

```{r}
sample_variables(MC_1)
rank_names(MC_1)
```
#Filter low abundance ASVs

```{r}
# Filter low-abundance taxa (≥10 total reads) and remove T0
MC_1_filt <- subset_samples(MC_1, Replicate != "T0" )

ps_filt_1 <- MC_1_filt  %>%
  prune_taxa(taxa_sums(.) >= 10, .)


```


```{r}
# Extract count matrix (taxa x samples)
counts <- as(otu_table(ps_filt_1), "matrix")

# Ensure taxa are rows
if (!taxa_are_rows(ps_filt_1)) {
  counts <- t(counts)
}

# Extract condition vector
meta <- data.frame(sample_data(ps_filt_1))
conds <- meta$Treatment
conds <- factor(conds)
conds
```

```{r}
conds <- as.character(meta$Treatment)

set.seed(123)
x <- aldex.clr(counts,
               conds,
               mc.samples = 512,
               denom = "all",
               verbose = TRUE)
```


```{r}
class(conds)
unique(conds)
length(conds) == ncol(counts)
```

#Overall test to check for differential abundance (Kruskal–Wallis, >2 groups)
```{r}
#Overall test (Kruskal–Wallis, >2 groups)
kw <- aldex.kw(x)


```

```{r}
#Significant asv
sig <- kw %>%
  as_tibble(rownames = "ASV") %>%
  filter(kw.eBH < 0.05)


```


```{r}
# Effect size
res_pw <- list()

for (trt in setdiff(unique(conds), "A")) {

  sel <- conds %in% c("A", trt)

  x_pw <- aldex.clr(counts[, sel],
                    conds[sel],
                    mc.samples = 512)

  res_pw[[trt]] <- data.frame(
    aldex.ttest(x_pw),
    aldex.effect(x_pw)
  )
}



write.csv(res_pw$C, file = "res_pw_MC1C.csv")

```

```{r}
sig <- res_pw$D %>%
  as_tibble(rownames = "ASV") %>%
  filter(wi.eBH < 0.05,
         abs(effect) > 1)



#use wilkinson due to small n and compositional data
#effect size > 1 according to developers to be considered meaningful
```


```{r}
ggplot(res_pw$D, aes(x = effect, y = -log10(we.eBH))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(title = "ALDEx2 Differential Abundance (18S)",
       x = "Effect size (CLR)",
       y = "-log10(FDR)")

```
#CLR + Aitchison distance + PERMANOVA FOR MC1 using counts (ps_filt)

```{r}
#install.packages("compositions")
library(compositions)
clr_mat <- t(apply(counts, 2, function(x) clr(x + 1)))

dist_aitch <- dist(clr_mat, method = "euclidean")

adonis2(dist_aitch ~ Treatment,
        data = meta,
        permutations = 999)

bd <- betadisper(dist_aitch, meta$Treatment)
permutest(bd)
```



```{r}

#Pairwise permanova
# Your filtered phyloseq object
ps <- ps_filt_1   # rename if needed

# Extract count matrix (taxa x samples)
counts <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) counts <- t(counts)

# Metadata
meta <- as.data.frame(sample_data(ps))
meta$Treatment <- as.character(meta$Treatment)

```


```{r}
meta <- sample_data(ps)
meta <- data.frame(meta, check.names = FALSE, stringsAsFactors = FALSE)

meta$Treatment <- as.factor(meta$Treatment)
```


```{r}
class(meta)
str(meta$Treatment)
```


```{r}
meta <- meta[rownames(clr_mat), , drop = FALSE]

all(rownames(meta) == rownames(clr_mat))
```


```{r}
sel <- meta$Treatment %in% c("A", "B")

meta_sub <- meta[sel, , drop = FALSE]
meta_sub <- as.data.frame(meta_sub)

dist_sub <- as.dist(as.matrix(dist_aitch)[sel, sel])

adonis2(dist_sub ~ Treatment, data = meta_sub, permutations = 999)

```


```{r}
control_name <- "A"
pairwise_results <- list()

for (trt in setdiff(levels(meta$Treatment), control_name)) {

  sel <- meta$Treatment %in% c(control_name, trt)

  meta_sub <- meta[sel, , drop = FALSE]
  meta_sub <- as.data.frame(meta_sub)

  dist_sub <- as.dist(as.matrix(dist_aitch)[sel, sel])

  pairwise_results[[paste(trt, "vs", control_name)]] <-
    adonis2(dist_sub ~ Treatment, data = meta_sub, permutations = 999)
}

pairwise_results

write.csv(pairwise_results$`D vs A`, file = "Pairwise_MC1_AD.csv")

```


## MC2
```{r}
sample_variables(MC_2)
rank_names(MC_2)
```
#Filter low abundnance ASVs

```{r}
# Filter low-abundance taxa (≥10 total reads) and remove T0
MC_2_filt <- subset_samples(MC_2, Replicate != "T0" )

ps_filt_2 <- MC_2_filt  %>%
  prune_taxa(taxa_sums(.) >= 10, .)


```


```{r}
# Extract count matrix (taxa x samples)
counts <- as(otu_table(ps_filt_2), "matrix")

# Ensure taxa are rows
if (!taxa_are_rows(ps_filt_2)) {
  counts <- t(counts)
}

# Extract condition vector
meta <- data.frame(sample_data(ps_filt_2))
conds <- meta$Treatment
conds <- factor(conds)
conds
```

```{r}
conds <- as.character(meta$Treatment)

set.seed(123)
x <- aldex.clr(counts,
               conds,
               mc.samples = 512,
               denom = "all",
               verbose = TRUE)

```


```{r}
class(conds)
unique(conds)
length(conds) == ncol(counts)
```

#Overall test to check for differential abundance (Kruskal–Wallis, >2 groups)
```{r}
#Overall test (Kruskal–Wallis, >2 groups)
kw <- aldex.kw(x)

```



```{r}
sig <- kw %>%
  as_tibble(rownames = "otu") %>%
  filter(kw.eBH < 0.05)


```

```{r}
# Effect size
res_pw <- list()

for (trt in setdiff(unique(conds), "A")) {

  sel <- conds %in% c("A", trt)

  x_pw <- aldex.clr(counts[, sel],
                    conds[sel],
                    mc.samples = 512)

  res_pw[[trt]] <- data.frame(
    aldex.ttest(x_pw),
    aldex.effect(x_pw)
  )
}



```

```{r}
write.csv(res_pw$D, file = "res_pw_MC2D.csv")
```

```{r}
sig <- res_pw$B %>%
  as_tibble(rownames = "ASV") %>%
  filter(wi.eBH < 0.05,
         abs(effect) > 1)



#use wilkinson due to small n and compositional data
```


```{r}
ggplot(res_pw$D, aes(x = effect, y = -log10(we.eBH))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(title = "ALDEx2 Differential Abundance (18S)",
       x = "Effect size (CLR)",
       y = "-log10(FDR)")

```
#CLR + Aitchison distance + PERMANOVA FOR MC2 using counts (ps_filt)

```{r}
#install.packages("compositions")
library(compositions)
clr_mat <- t(apply(counts, 2, function(x) clr(x + 1)))

dist_aitch <- dist(clr_mat, method = "euclidean")

adonis2(dist_aitch ~ Treatment,
        data = meta,
        permutations = 999)

bd <- betadisper(dist_aitch, meta$Treatment)
permutest(bd)
```



```{r}

#Pairwise permanova
# Your filtered phyloseq object
ps <- ps_filt_2   # rename if needed

# Extract count matrix (taxa x samples)
counts <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) counts <- t(counts)

# Metadata
meta <- as.data.frame(sample_data(ps))
meta$Treatment <- as.character(meta$Treatment)

```


```{r}
meta <- sample_data(ps)
meta <- data.frame(meta, check.names = FALSE, stringsAsFactors = FALSE)

meta$Treatment <- as.factor(meta$Treatment)
```


```{r}
class(meta)
str(meta$Treatment)
```


```{r}
meta <- meta[rownames(clr_mat), , drop = FALSE]

all(rownames(meta) == rownames(clr_mat))
```


```{r}
sel <- meta$Treatment %in% c("A", "D")

meta_sub <- meta[sel, , drop = FALSE]
meta_sub <- as.data.frame(meta_sub)

dist_sub <- as.dist(as.matrix(dist_aitch)[sel, sel])

adonis2(dist_sub ~ Treatment, data = meta_sub, permutations = 999)

```


```{r}
control_name <- "A"
pairwise_results <- list()

for (trt in setdiff(levels(meta$Treatment), control_name)) {

  sel <- meta$Treatment %in% c(control_name, trt)

  meta_sub <- meta[sel, , drop = FALSE]
  meta_sub <- as.data.frame(meta_sub)

  dist_sub <- as.dist(as.matrix(dist_aitch)[sel, sel])

  pairwise_results[[paste(trt, "vs", control_name)]] <-
    adonis2(dist_sub ~ Treatment, data = meta_sub, permutations = 999)
}

pairwise_results


write.csv(pairwise_results$B, file = "pairwis_MC2B.csv")
```



## MC3


```{r}
sample_variables(MC_3)
rank_names(MC_3)
```
#Filter low abundnance ASVs

```{r}
# Filter low-abundance taxa (≥10 total reads) and remove T0
MC_3_filt <- subset_samples(MC_3, Replicate != "T0" )

ps_filt_3 <- MC_3_filt  %>%
  prune_taxa(taxa_sums(.) >= 10, .)


```


```{r}
# Extract count matrix (taxa x samples)
counts <- as(otu_table(ps_filt_3), "matrix")

# Ensure taxa are rows
if (!taxa_are_rows(ps_filt_3)) {
  counts <- t(counts)
}

# Extract condition vector
meta <- data.frame(sample_data(ps_filt_3))
conds <- meta$Treatment
conds <- factor(conds)
conds
```


```{r}
conds <- as.character(meta$Treatment)

set.seed(123)
x <- aldex.clr(counts,
               conds,
               mc.samples = 512,
               denom = "all",
               verbose = TRUE)
```


```{r}
class(conds)
unique(conds)
length(conds) == ncol(counts)
```

#Overall test to check for differential abundance (Kruskal–Wallis, >2 groups)
```{r}
#Overall test (Kruskal–Wallis, >2 groups)
kw <- aldex.kw(x)

```

```{r}
sig <- kw %>%
  as_tibble(rownames = "otu") %>%
  filter(kw.eBH < 0.05)

```

```{r}
# Effect size
res_pw <- list()

for (trt in setdiff(unique(conds), "A")) {

  sel <- conds %in% c("A", trt)

  x_pw <- aldex.clr(counts[, sel],
                    conds[sel],
                    mc.samples = 512)

  res_pw[[trt]] <- data.frame(
    aldex.ttest(x_pw),
    aldex.effect(x_pw)
  )
}



```

```{r}
sig <- res_pw$D %>%
  as_tibble(rownames = "ASV") %>%
  filter(wi.eBH < 0.05,
         abs(effect) > 1)



#use wilkinson due to small n and compositional data
```

```{r}
write.csv(res_pw$D, file = "res_pw_MC3D.csv")

```


```{r}
ggplot(res_pw$C, aes(x = effect, y = -log10(wi.eBH))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(title = "ALDEx2 Differential Abundance (18S)",
       x = "Effect size (CLR)",
       y = "-log10(FDR)")

#ALDEx2 results are log-ratios, not fold-changes

#Significant taxa are changing relative to the community, not in absolute abundance

#With n=3 per group:

#Expect few taxa to pass FDR (false discovery rate)

#Effect sizes are more informative than p-values
```
```{r}
#CLR + Aitchison distance + PERMANOVA FOR MC3 using counts (ps_filt)
#install.packages("compositions")
library(compositions)
clr_mat <- t(apply(counts, 2, function(x) clr(x + 1)))

dist_aitch <- dist(clr_mat, method = "euclidean")

adonis2(dist_aitch ~ Treatment,
        data = meta,
        permutations = 999)

bd <- betadisper(dist_aitch, meta$Treatment)
permutest(bd)
```


```{r}
# Your filtered phyloseq object
ps <- ps_filt_3   # rename if needed

# Extract count matrix (taxa x samples)
counts <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) counts <- t(counts)

# Metadata
meta <- as.data.frame(sample_data(ps))
meta$Treatment <- as.character(meta$Treatment)

```


```{r}
meta <- sample_data(ps)
meta <- data.frame(meta, check.names = FALSE, stringsAsFactors = FALSE)

meta$Treatment <- as.factor(meta$Treatment)
```


```{r}
class(meta)
str(meta$Treatment)
```


```{r}
meta <- meta[rownames(clr_mat), , drop = FALSE]

all(rownames(meta) == rownames(clr_mat))
```


```{r}
sel <- meta$Treatment %in% c("A", "B")

meta_sub <- meta[sel, , drop = FALSE]
meta_sub <- as.data.frame(meta_sub)

dist_sub <- as.dist(as.matrix(dist_aitch)[sel, sel])


adonis2(dist_sub ~ Treatment, data = meta_sub, permutations = 999)

```

```{r}
control_name <- "A"
pairwise_results <- list()

for (trt in setdiff(levels(meta$Treatment), control_name)) {

  sel <- meta$Treatment %in% c(control_name, trt)

  meta_sub <- meta[sel, , drop = FALSE]
  meta_sub <- as.data.frame(meta_sub)

  dist_sub <- as.dist(as.matrix(dist_aitch)[sel, sel])

  pairwise_results[[paste(trt, "vs", control_name)]] <-
    adonis2(dist_sub ~ Treatment, data = meta_sub, permutations = 999)
}

pairwise_results

```
```{r}
write.csv(pairwise_results$C, file = "pairwisew_MC3C.csv")
```


## Diversity
```{r}
##Alpha diversIty

Richness <- estimate_richness(ps_MC, split = TRUE, measures=c("Shannon", "Chao1", "observed"))
write.csv(Richness, file = "Richness.csv")
```

#Diveristy was added to the metadata
# M1
```{r}
#Diversity was added to the metadata
Richeness_edited <- read_xlsx("Metadata.xlsx")

M1_MT <- filter(Richeness_edited, Experiment=="E1")
```
#Chao
```{r}
M1_diversity_Chao <- M1_MT %>% ggplot( aes(x=Number_2, y=Chao1, fill=Treatment)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.5, alpha=0.9) +
   ggtitle("MC1") + ylim(40, 100) + theme(text = element_text(size = 15))

M1_diversity_Chao <- M1_diversity_Chao + scale_fill_manual(values =c("#ffc100", "#00688b", "#e0301e","#602320","#999999"))

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


#Shannon
```{r}
M1_diversity_shannon <- M1_MT %>% ggplot( aes(x=Number_2, y=Shannon, fill=Treatment)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.5, alpha=0.9) +
   ggtitle("MC1") + ylim(1.5, 3) + theme(text = element_text(size = 15))

M1_diversity_shannon <- M1_diversity_shannon + scale_fill_manual(values =c("#ffc100", "#00688b", "#e0301e","#602320","#999999"))
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
