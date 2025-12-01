library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(sctransform)
library(patchwork)
library(GenomeInfoDb)
library(tidyverse)
library(grid)
library(readr)
library(harmony)
library(GenomicRanges)
library(SeuratData)
source("/data/lib14/R/Rscripts/utilities.R")
setwd("/data/Choi_lung/scLongreads/Multiomics_integration/")
peaks_annot <- read.csv("/data/Choi_lung/ChiaHan/Peak_and_motifBreaker_results/peaks_annotation_by_ChIPseeker_CHL_chranno.csv", header = TRUE, row.names = 1)
peaks_annot$peak <- paste(peaks_annot$seqnames, peaks_annot$start, peaks_annot$end, sep = "-")
# 
# NB.nominal.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_list.rds")
# 
# SNP.list <- lapply(NB.nominal.list, function(x) unique(x$snp))

# extract all snps located in RBP motif regions

peaks_region <-peaks_annot[, c(1,2,3,13)]
write.table(peaks_region, file = "/data/Choi_lung/scLongreads/Multiomics_integration/Multiome_peak_region.txt",
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

arg <- paste0("/usr/local/apps/plink/1.9.0-beta4.4/plink --bfile /data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype --extract Multiome_peak_region.txt --range --make-bed --out Multiome_peak --write-snplist")

system(arg)
SNPs <- read.table("/data/Choi_lung/scLongreads/Multiomics_integration/Multiome_peak.snplist")
SNPs <- SNPs$V1

library(snpStats)
fam <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.fam"
bim <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bim"
bed <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bed"
sample <- read.plink(bed, bim, fam)

snp_info <- sample$map[SNPs,c(1,2,4,5,6)]
colnames(snp_info) <- c("chr", "rsid", "pos", "alt", "ref")
cCRE.list <- list()
for (snp in SNPs) {
  cCRE <- peaks_annot$peak[which(peaks_annot$seqnames == snp_info[snp,]$chr & 
                                   as.numeric(peaks_annot$start) <= as.numeric(snp_info[snp,]$pos) & 
                                   as.numeric(peaks_annot$end) >= as.numeric(snp_info[snp,]$pos))]
  cCRE.list[[snp]] <- cCRE
}

save(list = c("snp_info","cCRE.list"),file = "/data/Choi_lung/scLongreads/Multiomics_integration/QTL_SNPs_integration.RData")


# summarize peak results
df <- data.frame(variants = names(cCRE.list), 
                 peak = unlist(cCRE.list))
NB.nominal.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_sig_list.rds")
isoQTL.list <- lapply(NB.nominal.sig.list, function(x) unique(x$snp))
isoQTLs <- unique(do.call(c, isoQTL.list))
sum(isoQTLs %in% df$variants)
# 6269080, 238678, 441277, 27555
enrich_tb <- matrix(c(6269080, 238678, 441277, 27555),
                    nrow = 2, byrow = TRUE,
                    dimnames = list(Total = c("all", "isoQTL"),
                                    peak = c("all", "isoQTL")))
ft <- fisher.test(enrich_tb)
ft <- chisq.test(enrich_tb)
log(ft$estimate)
ft$p.value


# 
promoter_peaks <- peaks_annot[grepl("Promoter", peaks_annot$annotation),]
TALON_afterqc_orf_secondpass2$promoter_peak <- NA

# Generate isoform ratio matrix
gff = as.data.frame(rtracklayer::import('/data/Choi_lung/scLongreads/TALON_workspace/test10s/NCI_lung_secondpass.isoforms_w_known.gtf'))
table(gff$type)
transcript_gff <- subset(gff, type == "transcript")
rm(gff)
pos <- data.frame(X = transcript_gff$transcript_id,
                  chr = transcript_gff$seqnames, 
                  start = transcript_gff$start, 
                  end = transcript_gff$end,
                  gene_id = transcript_gff$gene_id)
tss <- read.table("/data/Choi_lung/scLongreads/TALON_workspace/test10s/tss.bed", sep = "\t", header = FALSE)
rownames(tss) <- tss$V4
rownames(pos) <- pos$X


rownames(pos) <- pos$X
pos <- pos[TALON_afterqc_orf_secondpass2$annot_transcript_id,]
TALON_afterqc_orf_secondpass2$start <- pos$start
TALON_afterqc_orf_secondpass2$end <- pos$end
TALON_afterqc_orf_secondpass2$promoter_peak <- NA
for (i in 1:nrow(TALON_afterqc_orf_secondpass2)) {
  cCRE <- peaks_annot$peak[which(peaks_annot$seqnames == TALON_afterqc_orf_secondpass2$chrom.x[i] & 
                                   as.numeric(peaks_annot$start) <= as.numeric(TALON_afterqc_orf_secondpass2$start[i]) & 
                                   as.numeric(peaks_annot$end) >= as.numeric(TALON_afterqc_orf_secondpass2$start[i]))]
  if (length(cCRE) == 1) {
    TALON_afterqc_orf_secondpass2$promoter_peak[i] <- cCRE
  }
  
}
saveRDS(TALON_afterqc_orf_secondpass2, file = "/data/Choi_lung/scLongreads/TALON_workspace/test10s/TALON_afterqc_orf_secondpass_overlap_peaks.rds")
df_bg <- as.data.frame(table(TALON_afterqc_orf_secondpass2$transcript_novelty))
df_overlap_peak <- as.data.frame(table(TALON_afterqc_orf_secondpass2[!is.na(TALON_afterqc_orf_secondpass2$promoter_peak),]$transcript_novelty))
df_bg$percentage <- round((df_overlap_peak$Freq/df_bg$Freq)*100, 2)
colnames(df_bg) <- c("Structure category", "Number of transcripts", "TSS overlaps Multiome peaks")
pos_sub <- subset(pos, X %in% unique(final_table_AS_annot$phenotype_id))
pos_sub$promoter_peak <- NA
for (i in 1:nrow(pos_sub)) {
  cCRE <- peaks_annot$peak[which(peaks_annot$seqnames == pos_sub$chr[i] & 
                                   as.numeric(peaks_annot$start) <= as.numeric(pos_sub$TSS[i]) & 
                                   as.numeric(peaks_annot$end) >= as.numeric(pos_sub$TSS[i]))]
  if (length(cCRE) == 1) {
    pos_sub$promoter_peak[i] <- cCRE
  }
  
}
sum(is.na(pos_sub$promoter_peak))
saveRDS(TALON_afterqc_orf_secondpass2, file = "/data/Choi_lung/scLongreads/TALON_workspace/test10s/TALON_afterqc_orf_secondpass_overlap_peaks.rds")

peak_callin_by_ct <- read.csv("/data/Choi_lung/ChiaHan/Peak_and_motifBreaker_results/peaks_by_cell_types_CHL_chranno.csv", header = TRUE, row.names = 1)
peak_callin_by_ct$promoter_peak <- paste(peak_callin_by_ct$seqnames, peak_callin_by_ct$start, peak_callin_by_ct$end, sep = "-")
peak_called_in.list <- strsplit(peak_callin_by_ct$peak_called_in,",")
celltypes <- levels(share_seur$CellType)
peak_called_in.list <- lapply(peak_called_in.list, function(x) {
  i <- which(celltypes %in% x)
  idxs <- rep(FALSE, length(celltypes))
  idxs[i] <- TRUE
  return(idxs)
})
peak_called_in <- do.call(rbind, peak_called_in.list)
rownames(peak_called_in) <- peak_callin_by_ct$promoter_peak
colnames(peak_called_in) <- celltypes
ct_spec_peak <- peak_callin_by_ct$promoter_peak[which(rowSums(peak_called_in) == 1)]
peak_ct_only <- peak_called_in[ct_spec_peak,]
colSums(peak_ct_only)
tmp <- left_join(TALON_afterqc_orf_secondpass2, peak_callin_by_ct, by = "promoter_peak")
tmp_sub <- subset(tmp, promoter_peak %in% ct_spec_peak)

tmp_sub_AT2 <- subset(tmp_sub, peak_called_in == "AT2")
DEI_AT2 <- DEI_sig.list[["AT2"]]
DEI_AT2_unique<- DEI_unique.list[["AT2"]]
tmp_sub_AT2$transcript_name_unique[which(tmp_sub_AT2$transcript_name_unique %in% DEI_AT2$X)]
tmp_sub_AT2$transcript_name_unique[which(tmp_sub_AT2$transcript_name_unique %in% DEI_AT2_unique)]


tmp_sub_AT1 <- subset(tmp_sub, peak_called_in == "AT1")
DEI_AT1 <- DEI_sig.list[["AT1"]]
DEI_AT1_unique<- DEI_unique.list[["AT1"]]
tmp_sub_AT1$transcript_name_unique[which(tmp_sub_AT1$transcript_name_unique %in% DEI_AT1$X)]
tmp_sub_AT1$transcript_name_unique[which(tmp_sub_AT1$transcript_name_unique %in% DEI_AT1_unique)]

tmp_sub_NK <- subset(tmp_sub, peak_called_in == "NK")
DEI_NK <- DEI_sig.list[["NK_cells"]]
DEI_NK_unique<- DEI_unique.list[["NK_cells"]]
tmp_sub_NK$transcript_name_unique[which(tmp_sub_NK$transcript_name_unique %in% DEI_NK$X)]
tmp_sub_NK$transcript_name_unique[which(tmp_sub_NK$transcript_name_unique %in% DEI_NK_unique)]

tmp_sub_SMC <- subset(tmp_sub, peak_called_in == "SMC")
DEI_SMC <- DEI_sig.list[["SMC"]]
DEI_SMC_unique<- DEI_unique.list[["SMC"]]
tmp_sub_SMC$transcript_name_unique[which(tmp_sub_SMC$transcript_name_unique %in% DEI_SMC$X)]
tmp_sub_SMC$transcript_name_unique[which(tmp_sub_SMC$transcript_name_unique %in% DEI_SMC_unique)]

save(list = ls(),file = "/data/Choi_lung/scLongreads/Multiomics_integration/Multiple_comparison.RData")
