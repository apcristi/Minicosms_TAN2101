library(dada2)
library(Biostrings)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(seqinr)

setwd("/projet/externe/univ/agutierrez/Antonia/dada2/TAN2101_MC")

project_name <- "TAN2101_Minicosm"

fastq_dir <-    "/projet/externe/univ/agutierrez/Antonia/dada2/TAN2101_MC/fastq"
filtered_dir <- "./fastq_filtered/"
qual_dir <-     "./qual_pdf/"      
dada2_dir <-    "./dada2/"     
blast_dir <-    "./blast/"           
database_dir <- "/projet/externe/univ/agutierrez/Antonia/dada2/TAN1802/database/pr2_version_4.12.0_18S_dada2.fasta.gz" 

dir.create(filtered_dir)
dir.create(qual_dir)
dir.create(dada2_dir)
dir.create(blast_dir)

primer_set_fwd = c("CTTTCCCTACACGACGCTCTTCCGATCTCCAGCASCYGCGGTAATTCC")
primer_set_rev = c("GGAGTTCAGACGTGTGCTCTTCCGATCTACTTTCGTTCTTGATYRATGA")

primer_length_fwd <- str_length(primer_set_fwd[1]) 
primer_length_rev <- str_length(primer_set_rev[1])

DB_levels <- c("Kingdom", "Supergroup","Division", "Class", "Order", "Family", "Genus", "Species")

fns <- sort(list.files(fastq_dir, full.names = TRUE)) 
fns <- fns[str_detect( basename(fns),".fastq.gz")]
fns_R1 <- fns[str_detect( basename(fns),"R1")]
fns_R2 <- fns[str_detect( basename(fns),"R2")]


sample.names <- str_split(basename(fns_R1), n=2, pattern = "_", simplify = TRUE) 
sample.names <- sample.names[,1]

filt_dir <- file.path(fastq_dir, "filtered")
filt_R1 <- str_c(filtered_dir, sample.names, "_R1_filt.fastq")
filt_R2 <- str_c(filtered_dir, sample.names, "_R2_filt.fastq")


df <- data.frame()  
for(i in 1:length(fns_R1)) {
    geom <- fastq.geometry(fns_R1[i])
    df_one_row <- data.frame (n_seq=geom[1], file_name=basename(fns[i]))
    df <- bind_rows(df, df_one_row)
    print(paste("Finished with file", fns[i], ". ", round(i/length(fns)*100, 2), "%", sep=""))
    if(i == length(fns)) {print(paste("Finished"))}
} 
df

for(i in 1:length(fns)) {
    p1 <- plotQualityProfile(fns[i])
    if (i <= 2) {print(p1)}
    p1_file <- paste0(qual_dir, basename(fns[i]),".qual.pdf")
    ggsave( plot=p1, filename= p1_file, device = "pdf", width = 15, height = 15, scale=1, units="cm")
    print(paste("Finished with file", fns[i], ". ", round(i/length(fns)*100, 2), "%", sep=""))
    if(i == length(fns)) {print(paste("Finished"))}
}   


out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, trimLeft = c(primer_length_fwd, primer_length_rev), maxN=0, maxEE=c(2, 2), truncQ=10, rm.phix=TRUE, compress=FALSE, multithread=FALSE) 

err_R1 <- learnErrors(filt_R1, multithread=FALSE)
plotErrors(err_R1, nominalQ=TRUE)
pdf("LearnErrors_R1.pdf", width=11.69, height=8.27, paper='special')
plot(plotErrors(err_R1, nominalQ=TRUE))
dev.off()

err_R2 <- learnErrors(filt_R2, multithread=FALSE)
plotErrors(err_R2, nominalQ=TRUE)
pdf("LearnErrors_R2.pdf", width=11.69, height=8.27, paper='special')
plot(plotErrors(err_R2, nominalQ=TRUE))
dev.off()


derep_R1 <- derepFastq(filt_R1, verbose=FALSE)
derep_R2 <- derepFastq(filt_R2, verbose=FALSE)


names(derep_R1) <- sample.names
names(derep_R2) <- sample.names 


dada_R1 <- dada(derep_R1, err=err_R1, multithread=FALSE)
dada_R2 <- dada(derep_R2, err=err_R2, multithread=FALSE)

dada_R1[[1]]
dada_R2[[1]]


mergers <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2, verbose=TRUE)

 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
t_seqtab <- t(seqtab)

table(nchar(getSequences(seqtab)))

write_delim(data.frame(seqtab), str_c(dada2_dir, "table_seq.txt"))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)

paste0("% of non chimeras : ",sum(seqtab.nochim)/sum(seqtab))
paste0("total number of sequences : ",sum(seqtab.nochim))

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada_R1, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim), round((rowSums(seqtab.nochim)/out[,1])*100,2))

colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim", "final_pct")
rownames(track) <- sample.names

write.csv(track, str_c(dada2_dir, "track.txt"))


DB_file <- paste0(database_dir)
taxa <- assignTaxonomy(seqtab.nochim, refFasta=DB_file, taxLevels = DB_levels, minBoot = 0, outputBootstraps = TRUE, verbose = TRUE)

saveRDS(taxa, str_c(dada2_dir, "TAN2101_MC.taxa.rds"))


print(paste("__##__assignTaxonomy finished__##__"))


taxa <-  readRDS(str_c(dada2_dir, "TAN2101_MC.taxa.rds"))  
write.csv(as.tibble(taxa$tax), str_c(dada2_dir, "taxa.txt"))  
write.csv(as.tibble(taxa$boot), str_c(dada2_dir, "taxa_boot.txt"))
write.csv(seqtab.nochim, str_c(dada2_dir, "TAN2101_MC_seqtab2.txt"))

taxa_tax <- as.data.frame(taxa$tax)
taxa_boot <- as.data.frame(taxa$boot)

taxa_tax <- taxa_tax %>% rownames_to_column(var = "sequence") %>%
    rowid_to_column(var = "OTUNumber") %>%
    mutate(OTUNumber = sprintf("otu%04d", OTUNumber))
row.names(taxa_tax) <- taxa_tax$OTUNumber
row.names(taxa_boot) <- taxa_tax$OTUNumber  


seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim))
row.names(seqtab.nochim_trans) <- taxa_tax$OTUNumber

bootstrap_min <- 80
taxa_tax_18S <- taxa_tax[taxa_boot$Supergroup >= bootstrap_min,]
taxa_boot_18S <- taxa_boot[taxa_boot$Supergroup >= bootstrap_min,]
seqtab.nochim_18S <- seqtab.nochim_trans[taxa_boot$Supergroup >= bootstrap_min,]
write_tsv(as.tibble(seqtab.nochim_18S, header=T, row.names=T), str_c(dada2_dir, "seqtab_nonchim_18S.txt"))
 
dada2_database <-   bind_cols(taxa_tax_18S, seqtab.nochim_18S)
write_tsv(dada2_database, str_c(dada2_dir, "TAN2101_MC.database.tsv"))
	
df <-  dada2_database %>%  mutate(sequence = str_replace_all(sequence, "(-|\\.)",""))

seq_out <- Biostrings::DNAStringSet(df$sequence)

names(seq_out) <- str_c(df$OTUNumber,
                                df$Supergroup,
                                df$Division,
                                df$Class,
                                df$Order,
                                df$Family,
                                df$Genus,
                                df$Species,
                                sep="|")

Biostrings::writeXStringSet(seq_out, str_c(blast_dir, "TAN2101_MC_ASV.fasta"), compress=FALSE, width = 20000)
  



