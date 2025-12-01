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

load("/data/Choi_lung/scLongreads/eisoQTL.RData")
# Obtain cCRE called-in matrix
share_seur <- readRDS("/data/Choi_lung/ChiaHan/CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds")
rm(CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB)
peaks_by_ct <- read.csv("/data/Choi_lung/ChiaHan/peaks_by_cell_types_CHL_chranno.csv", header = TRUE, row.names = 1)
peaks_annot <- read.csv("/data/Choi_lung/ChiaHan/Peak_and_motifBreaker_results/peaks_annotation_by_ChIPseeker_CHL_chranno.csv", header = TRUE, row.names = 1)
peak_called_in.list <- strsplit(peaks_by_ct$peak_called_in,",")
celltypes <- as.character(unique(share_seur$CellType))
peak_called_in.list <- lapply(peak_called_in.list, function(x) {
  i <- which(celltypes %in% x)
  idxs <- rep(FALSE, length(celltypes))
  idxs[i] <- TRUE
  return(idxs)
})
peak_called_in <- do.call(rbind, peak_called_in.list)

colnames(peak_called_in) <- celltypes
rownames(peak_called_in) <- paste(peaks_by_ct$seqnames,peaks_by_ct$start, peaks_by_ct$end, sep = "-")
# peak_called_in_sub <- peak_called_in[unique(GWAS_subset$colocalized_cCRE),]
# colSums(peak_called_in_sub)
# 
# ct_specific_ccre <- peak_called_in_sub[which(rowSums(peak_called_in_sub) == 1),]
# colSums(ct_specific_ccre)

# Load NB significant isoQTLs
NB.nominal.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_sig_list.rds")
NB.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")
NB.nominal.sig.list <- NB.nominal.sig.list[names(NB.sig.list)]
category <- c(rep("EPI",8),rep("IMMUNE",15),rep("EC", 6), rep("STROMA",4))
df.list <- list()
for (i in 1:length(NB.nominal.sig.list)) {
  nominal.rst <- NB.nominal.sig.list[[i]]
  dist <- as.numeric(nominal.rst$tss_distance/1000)
  df <- data.frame(dist = dist,
                   celltype = (names(NB.nominal.sig.list)[i]),
                   category = category[i])
  df.list[[i]] <- df
}
df_final <- Reduce(rbind, df.list)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)

fc <- colorRampPalette(c("#FA877533", "#610C04"))
fc_yellow <- colorRampPalette(c("#8B4000", "#FFD70033"))
fc_purple <- colorRampPalette(c("lightpink","#4B0082"))
fc_blue <- colorRampPalette(c("lightblue", "darkblue"))
fc_yellow(9)
cellcolors <- c(fc_yellow(9)[c(3,2,1,4,8,9,7,5,6)], fc(17)[c(1,17,15,16,14,2,13,4,11,
                                                             3,12,5,
                                                             6,7,9,10,8)],
                fc_purple(6)[c(4,2,6,3,5,1)], fc_blue(5)[c(5,1,2,4,3)])
plot(rep(1, 37),col = cellcolors, pch = 19, cex = 3)
cols <- cellcolors[-c(9,25,26,37)]
df_final$celltype <- factor(df_final$celltype, levels = names(NB.nominal.sig.list))
df_final$category <- factor(df_final$category, levels = c("EPI", "IMMUNE", "EC", "STROMA"))
ggplot(data=df_final, aes(x=dist, group=celltype, colour = celltype)) +
  geom_density(lwd = .75, linetype = 1) + xlab("Distance from TSS") +
  scale_color_manual(values = cols) + theme_bw()+
  facet_wrap(~category)

load("/data/Choi_lung/scLongreads/eisoQTL.RData")
final_table <- read.table("/data/lib14/project/scLongread/Tables/NB_eIsoform_sum.txt", sep = "\t",
                          header = TRUE)
length(unique(final_table$variant_id))
head(snp_info)

peaks_annot$peak <- paste(peaks_annot$seqnames, peaks_annot$start, peaks_annot$end, sep = "-")
peaks_annot$peak_called_id <- peaks_by_ct$peak_called_in

leadisoQTL_int_multiome.list <- list()
for (celltype in names(NB.sig.list)) {
  leadisoQTL <- NB.sig.list[[celltype]]
  leadisoQTL$peak <- NA
  leadisoQTL$Celltype <- celltype
  leadisoQTL <- leadisoQTL %>% group_by(phenotype_id) %>% mutate(transcript_name = Search_transcript_name2(phenotype_id),
                                                                 .after = phenotype_id)
  leadisoQTL <- leadisoQTL %>% group_by(variant_id) %>% mutate(chr = paste0("chr", chrom),
                                                               pos = snp_info[variant_id,]$pos,
                                                               ref = snp_info[variant_id,]$ref,
                                                               alt = snp_info[variant_id,]$alt,
                                                               .after = variant_id)
  for(i in 1:nrow(leadisoQTL)){
    cCRE <- peaks_annot$peak[which(peaks_annot$seqnames == leadisoQTL$chr[i] & 
                                     as.numeric(peaks_annot$start) <= as.numeric(leadisoQTL$pos[i]) & 
                                     as.numeric(peaks_annot$end) >= as.numeric(leadisoQTL$pos[i]))]
    if (length(cCRE) == 1) {
      leadisoQTL$peak[i] <- cCRE
    }
  }
  leadisoQTL_subset <- leadisoQTL[which(!is.na(leadisoQTL$peak)),]
  leadisoQTL_subset <- left_join(leadisoQTL_subset, peaks_annot, by = "peak")
  leadisoQTL_int_multiome.list[[celltype]] <- leadisoQTL_subset
  
}
leadisoQTL_int_multiome <- Reduce(rbind, leadisoQTL_int_multiome.list)
test <- (leadisoQTL_int_multiome[which(leadisoQTL_int_multiome$variant_id %in% ct_spec_isoQTLs),])

AllisoQTL_int_multiome.list <- list()
for (celltype in names(NB.nominal.sig.list)) {
  leadisoQTL <- NB.nominal.sig.list[[celltype]]
  leadisoQTL$peak <- NA
  leadisoQTL$chr <- paste0("chr", leadisoQTL$chrom)
  leadisoQTL$Celltype <- celltype
  leadisoQTL <- leadisoQTL %>% group_by(phenotype_id) %>% mutate(transcript_name = Search_transcript_name2(phenotype_id),
                                                                 .after = phenotype_id)
  
  for(i in 1:nrow(leadisoQTL)){
    cCRE <- peaks_annot$peak[which(peaks_annot$seqnames == leadisoQTL$chr[i] & 
                                     as.numeric(peaks_annot$start) <= as.numeric(leadisoQTL$pos[i]) & 
                                     as.numeric(peaks_annot$end) >= as.numeric(leadisoQTL$pos[i]))]
    if (length(cCRE) == 1) {
      leadisoQTL$peak[i] <- cCRE
    }

  }
  leadisoQTL_subset <- leadisoQTL[which(!is.na(leadisoQTL$peak)),]
  leadisoQTL_subset <- left_join(leadisoQTL_subset, peaks_annot, by = "peak")
  AllisoQTL_int_multiome.list[[celltype]] <- leadisoQTL_subset
  
}

write.table(unique(final_table$variant_id), file = "/data/Choi_lung/scLongreads/TWAS/1000Genome/leadisoQTL_rsids.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# plink2 --bfile 1kG_EUR --freq --allow-extra-chr --extract leadisoQTL_rsids.txt --out leadisoQTL_EUR
# plink2 --bfile 1kG_EAS --freq --allow-extra-chr --extract leadisoQTL_rsids.txt --out leadisoQTL_EAS
# plink2 --bfile 1kG_AFR --freq --allow-extra-chr --extract leadisoQTL_rsids.txt --out leadisoQTL_AFR

NB.nominal.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_list.rds")

SNP.list <- lapply(NB.nominal.list, function(x) unique(x$snp))

library(snpStats)
fam <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.fam"
bim <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bim"
bed <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bed"
sample <- read.plink(bed, bim, fam)

SNPs <- unique(do.call(c, SNP.list))
snp_info <- sample$map[SNPs,c(1,2,4,5,6)]
colnames(snp_info) <- c("chr", "rsid", "pos", "alt", "ref")
cCRE.list <- list()
for (snp in SNPs) {
  cCRE <- peaks_annot$peak[which(peaks_annot$seqnames == snp_info[snp,]$chr & 
                                   as.numeric(peaks_annot$start) <= as.numeric(snp_info[snp,]$pos) & 
                                   as.numeric(peaks_annot$end) >= as.numeric(snp_info[snp,]$pos))]
  cCRE[[snp]] <- cCRE
}



cCRE_mtx <- apply(snp_info[c(1000:1500),], 1, function(x){
  cCRE <- peaks_annot$peak[which(peaks_annot$seqnames == x & 
                                   as.numeric(peaks_annot$start) <= as.numeric(x[,3]) & 
                                   as.numeric(peaks_annot$end) >= as.numeric(x[,3]))]
})




AllSNPs_int_multiome.list <- list()
for (celltype in names(NB.nominal.list)) {
  leadisoQTL <- NB.nominal.list[[celltype]]
  leadisoQTL$peak <- NA
  leadisoQTL$chr <- paste0("chr", leadisoQTL$chrom)
  leadisoQTL$Celltype <- celltype
  leadisoQTL <- leadisoQTL %>% group_by(phenotype_id) %>% mutate(transcript_name = Search_transcript_name2(phenotype_id),
                                                                 .after = phenotype_id)
  
  for(i in 1:nrow(leadisoQTL)){
    cCRE <- peaks_annot$peak[which(peaks_annot$seqnames == leadisoQTL$chr[i] & 
                                     as.numeric(peaks_annot$start) <= as.numeric(leadisoQTL$pos[i]) & 
                                     as.numeric(peaks_annot$end) >= as.numeric(leadisoQTL$pos[i]))]
    if (length(cCRE) == 1) {
      leadisoQTL$peak[i] <- cCRE
    }
    
  }
  leadisoQTL_subset <- leadisoQTL[which(!is.na(leadisoQTL$peak)),]
  leadisoQTL_subset <- left_join(leadisoQTL_subset, peaks_annot, by = "peak")
  AllSNPs_int_multiome.list[[celltype]] <- leadisoQTL_subset
  
}

snp_subset <- snp_info[unique(final_table$variant_id),]
for (chr in 1:22) {
  chrchr <- paste0("chr",chr)
  chr_snps <- subset(snp_subset, chr == chrchr)
  snps <- chr_snps$rsid
  write.table(snps, file = paste0("/data/Choi_lung/scLongreads/TWAS/FUSION/1000genome_phase3/1000G_EUR_Phase3_plink/", chr,"_rsids.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  arg <- paste0("/usr/local/apps/plink/1.9.0-beta4.4/plink --bfile 1000G.EUR.QC.", chr,
                " --keep-allele-order --freq --extract ", 
                              chr, "_rsids.txt --out chr_",
                              chr, "_rsids_1kg_EUR")
  setwd("/data/Choi_lung/scLongreads/TWAS/FUSION/1000genome_phase3/1000G_EUR_Phase3_plink/")
  # plink --bfile FLCCA_controls_R05_MAF001 --keep-allele-order --r square --extract 6q22.1_snps.txt --out locus_6q22.1_snp_flcca_ctrl --write-snplist
  system(arg)
  arg <- paste0("/usr/local/apps/plink/1.9.0-beta4.4/plink --bfile /data/Choi_lung/scLongreads/TWAS/FUSION/1000genome_phase3/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.", chr,
                " --keep-allele-order --freq --extract ", 
                chr, "_rsids.txt --out chr_",
                chr, "_rsids_1kg_EAS")
  system(arg)
}
EAS <- NULL
EUR <- NULL
for (chr in 1:22) {
  EAS_path <- paste0("/data/Choi_lung/scLongreads/TWAS/FUSION/1000genome_phase3/1000G_EUR_Phase3_plink/chr_", chr,"_rsids_1kg_EAS.frq")
  EUR_path <- paste0("/data/Choi_lung/scLongreads/TWAS/FUSION/1000genome_phase3/1000G_EUR_Phase3_plink/chr_", chr,"_rsids_1kg_EUR.frq")
  EAS_frq <- read.table(EAS_path, sep = "", header = TRUE)
  EUR_frq <- read.table(EUR_path, sep = "", header = TRUE)
  EAS <- rbind(EAS, EAS_frq)
  EUR <- rbind(EUR, EUR_frq)
}
EAS$MAF_EAS <- EAS$MAF
rownames(EAS) <- EAS$SNP
rownames(EUR) <- EUR$SNP
EAS[rownames(EUR),]$MAF_EUR <- EUR$MAF

test <- left_join(EAS,EUR, by = "SNP", suffix = c(".EAS", ".EUR"))
test <- test[!is.na(test$CHR.EUR),]
test$MAF_EAS[which(test$A1.EAS != test$A1.EUR)] <- (1 - test$MAF.EAS[which(test$A1.EAS != test$A1.EUR)])
test$MAF_diff <- (test$MAF_EAS - test$MAF.EUR)
idx1 <- which(test$MAF_EAS > 0.5 & test$MAF_diff > 0.2) # minor allele is different in EAS and EUR
idx2 <- which(test$MAF.EUR < 0.05 & test$MAF_diff > 0.2)
idx <- union(idx1, idx2)
tmp <- test[idx,]
final_table_EAS <- subset(final_table, variant_id %in% tmp$SNP)
View(final_table_EAS[is.na(final_table_EAS$If_sGenes),])
write.table(test, file = "/data/Choi_lung/scLongreads/TWAS/FUSION/1000genome_phase3/LeadisoQTL_1kG_EAS_EUR_compare.txt", row.names = FALSE, quote = FALSE, sep = "\t")


EUR <- read.table("/data/Choi_lung/scLongreads/TWAS/1000Genome/leadisoQTL_EUR.afreq", 
                  sep = "", header = TRUE, comment.char = " ")
EAS <- read.table("/data/Choi_lung/scLongreads/TWAS/1000Genome/leadisoQTL_EAS.afreq", 
                  sep = "", header = TRUE, comment.char = " ")
AFR <- read.table("/data/Choi_lung/scLongreads/TWAS/1000Genome/leadisoQTL_AFR.afreq", 
                  sep = "", header = TRUE, comment.char = " ")
EUR$ALT_FREQS_EUR <- EUR$ALT_FREQS
EUR$ALT_FREQS_EAS <- EAS$ALT_FREQS
EUR$ALT_FREQS_AFR <- AFR$ALT_FREQS
EUR$id_ref_alt <- paste(EUR$ID, EUR$REF, EUR$ALT, sep = "_")

EUR <- EUR[!grepl("chr", EUR$X.CHROM),]
rownames(EUR) <- EUR$id_ref_alt
snp_subset$id_ref_alt <- paste(snp_subset$rsid,snp_subset$ref, snp_subset$alt, sep = "_")
rownames(snp_subset) <- snp_subset$id_ref_alt
length(intersect(EUR$id_ref_alt, snp_subset$id_ref_alt))
EUR$ALT_FREQS <- NULL
snp_info_w_1kG <- left_join(snp_subset, EUR, by = "id_ref_alt")
write.table(snp_info_w_1kG, file = "/data/Choi_lung/scLongreads/TWAS/1000Genome/LeadisoQTL_1kG_EUR_EAS_AFR_compare.txt", row.names = FALSE, quote = FALSE, sep = "\t")
