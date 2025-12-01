# Prepare RBP motifs from eCLIP-seq (Van Nostrand et al 2020 Nature)
# For enrichment analysis

# Bolun Li
# Feb 27th, 2025

library(dplyr)
library(ggplot2)
library(stringi)
library(stringr)
source("/data/lib14/R/Rscripts/utilities.R")

setwd("/data/Choi_lung/scLongreads/Enrichment/")
files <- list.files("/data/Choi_lung/scLongreads/Enrichment/", pattern = "qced.hg38.bed")
files <- list.files("/data/Choi_lung/scLongreads/Enrichment/", pattern = ".hg38.bed")
colnames <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
bed.list <- list()
bed.list1 <- list()
RBP_names <- NULL
for (file in files) {
  bed <- read.table(file, sep = "\t", header = FALSE)
  colnames(bed) <- colnames
  RBP_names <- c(RBP_names, unique(bed$name))
  bed.list1[[str_split_fixed(file, "\\.", 3)[,1]]] <- bed
}

names(bed.list) <- RBP_names
bed.list <- bed.list[-names(bed.list1)]
# QC with stringent cutoff
bed.qc.list <- lapply(bed.list, function(x){
  x$RBP <- str_split_fixed(x$name, "_", 3)[,1]
  x$cell_line <- str_split_fixed(x$name, "_", 3)[,2]
  x$sample_info <- str_split_fixed(x$name, "_", 3)[,3]
  x$peak <- paste(x$chr, x$start, x$end, sep = "-")
  # x <- subset(x, (signalValue >= 3) & (pValue >=5))
  return(x)
})

K562 <- bed.qc.list[grepl("K562", RBP_names)]
names(K562) <- str_split_fixed(names(K562), "_", 3)[,1]
K562_mtx <- do.call(rbind, K562)
HepG2 <- bed.qc.list[!grepl("K562", RBP_names)]
names(HepG2) <- str_split_fixed(names(HepG2), "_", 3)[,1]
HepG2_mtx <- do.call(rbind, HepG2)
tmp <- sapply(HepG2, function(x)nrow(x))
sum(tmp)
tmp2 <- lapply(HepG2, function(x)(x$peak))
tmp2 <- unique(Reduce(c,tmp2))


RBPs <- unique(str_split_fixed(RBP_names, "_", 3)[,1])
RBPs <- names(K562)
RBP_info <- read.csv("/data/Choi_lung/scLongreads/Enrichment/RBP_annotation.csv", header = TRUE)
rownames(RBP_info) <- RBP_info$RBP.name
RBP_info <- RBP_info[RBPs,c(1:21)]
Exprement_id <- read.table("/data/Choi_lung/scLongreads/Enrichment/ENCSR456FVU_metadata.tsv", header = TRUE, sep = "\t")
Exprement_id <- subset(Exprement_id, File.type == "bed")
Exprement_id <- Exprement_id[c(1:223)*3,]
# write.table(Exprement_id$File.accession, "/data/Choi_lung/scLongreads/Enrichment/IDR.list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
library(snpStats)
fam <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.fam"
bim <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bim"
bed <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bed"
sample <- read.plink(bed, bim, fam)


NB.nominal.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_sig_list.rds")
isoQTL.list <- lapply(NB.nominal.sig.list, function(x) unique(x$snp))
isoQTLs <- unique(do.call(c, isoQTL.list))
snp_info <- sample$map[isoQTLs,c(1,2,4,5,6)]
colnames(snp_info) <- c("chr", "rsid", "pos", "alt", "ref")
isoQTL_win_RBP_motif_peaks <- lapply(isoQTLs, function(snp){
  cCRE <- K562_mtx$peak[which(K562_mtx$chr == snp_info[snp,]$chr & 
                                as.numeric(K562_mtx$start) <= as.numeric(snp_info[snp,]$pos) & 
                                as.numeric(K562_mtx$end) >= as.numeric(snp_info[snp,]$pos))]
  RBP <- K562_mtx$RBP[which(K562_mtx$chr == snp_info[snp,]$chr & 
                              as.numeric(K562_mtx$start) <= as.numeric(snp_info[snp,]$pos) & 
                              as.numeric(K562_mtx$end) >= as.numeric(snp_info[snp,]$pos))]
  return(paste(cCRE,RBP, sep = "_"))
})
names(isoQTL_win_RBP_motif_peaks) <- isoQTLs
save(list = c("snp_info", "isoQTL_win_RBP_motif_peaks"), file = "/data/Choi_lung/scLongreads/Enrichment/isoQTL_RBP_motif.RData")
save(list = c("K562_mtx", "K562"), file = "/data/Choi_lung/scLongreads/Enrichment/K562_qced.RData")
save(list = c("HepG2_mtx", "HepG2"), file = "/data/Choi_lung/scLongreads/Enrichment/HepG2_qced.RData")


load("/data/Choi_lung/scLongreads/Enrichment/isoQTL_RBP_motif.RData")
load("/data/Choi_lung/scLongreads/Enrichment/K562_qced.RData")
load("/data/Choi_lung/scLongreads/Enrichment/HepG2_qced.RData")
isoQTLs <- names(isoQTL_win_RBP_motif_peaks)
test.list <- lapply(isoQTLs, function(isoQTL){
  if (length(isoQTL_win_RBP_motif_peaks[[isoQTL]])==0) {
    return(isoQTL_win_RBP_motif_peaks[[isoQTL]])
  }else{
    x <- isoQTL_win_RBP_motif_peaks[[isoQTL]]
    df <- data.frame(RBP_motif = str_split_fixed(x, "_", 2)[,1],
                     RBP = str_split_fixed(x, "_", 2)[,2],
                     snp = isoQTL)
    return(df)
  }
})
test <- do.call(rbind, test.list)


# extract all snps located in RBP motif regions

K562_RBP_region <-K562_mtx[, c(1,2,3,11)]
table(K562_RBP_region$chr)
K562_RBP_region <- K562_RBP_region[!(K562_RBP_region$chr %in% c("chr1_KI270706v1_random","chr14_GL000009v2_random", "chr22_KI270879v1_alt","chr4_GL000008v2_random", "chrUn_KI270742v1")),]
write.table(K562_RBP_region, file = "/data/Choi_lung/scLongreads/Enrichment/K562/K562_RBP_region.txt",
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

arg <- paste0("/usr/local/apps/plink/1.9.0-beta4.4/plink --bfile /data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype --extract ./K562/K562_RBP_region.txt --range --make-bed --out K562_RBP_motif --write-snplist")

system(arg)
RBP_bg_snps <- read.table("/data/Choi_lung/scLongreads/Enrichment/K562_RBP_motif.snplist")
library(snpStats)
fam <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.fam"
bim <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bim"
bed <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bed"
sample <- read.plink(bed, bim, fam)

snp_info <- sample$map[RBP_bg_snps$V1,c(1,2,4,5,6)]
colnames(snp_info) <- c("chr", "rsid", "pos", "alt", "ref")
RBP_bg_snps_peaks.list <- list()
RBP_bg_snps_peaks.list <- lapply(RBP_bg_snps$V1, function(snp){
  cCRE <- K562_mtx$peak[which(K562_mtx$chr == snp_info[snp,]$chr & 
                                as.numeric(K562_mtx$start) <= as.numeric(snp_info[snp,]$pos) & 
                                as.numeric(K562_mtx$end) >= as.numeric(snp_info[snp,]$pos))]
  RBP <- K562_mtx$RBP[which(K562_mtx$chr == snp_info[snp,]$chr & 
                              as.numeric(K562_mtx$start) <= as.numeric(snp_info[snp,]$pos) & 
                              as.numeric(K562_mtx$end) >= as.numeric(snp_info[snp,]$pos))]
  return(paste(cCRE,RBP, sep = "_"))
})
names(RBP_bg_snps_peaks.list) <- RBP_bg_snps$V1
RBP_bg_snps_peaks_df.list <- lapply(RBP_bg_snps$V1, function(snp){
  if (length(RBP_bg_snps_peaks.list[[snp]])==0) {
    return(RBP_bg_snps_peaks.list[[snp]])
  }else{
    x <- RBP_bg_snps_peaks.list[[snp]]
    df <- data.frame(RBP_motif = str_split_fixed(x, "_", 2)[,1],
                     RBP = str_split_fixed(x, "_", 2)[,2],
                     snp = snp)
    return(df)
  }
})
RBP_bg_snps_peaks_df <- do.call(rbind, RBP_bg_snps_peaks_df.list)
RBP_bg_sum <- as.data.frame(table(RBP_bg_snps_peaks_df$RBP))
RBP_isoQTL_sum <- as.data.frame(table(test$RBP))


# spliceosome
RBP_info$RBP.name[RBP_info$Spliceosome == 1]

# Try U2AF1 # 238678, 17, 6269080, 111
length(isoQTLs) # 238678 # total snps with MAF > 0.05: 6269080
length(RBP_bg_snps$V1) # 111
enrich_tb <- matrix(c(6269080, 238678, 111, 17),
                       nrow = 2, byrow = TRUE,
                       dimnames = list(Population = c("isoQTL", "background"),
                                       Allele = c("Total", "RBP")))
ft <- fisher.test(enrich_tb)
log(ft$estimate)

RBP_isoQTL_enrich_tb <- left_join(RBP_isoQTL_sum, RBP_bg_sum, by = "Var1")
colnames(RBP_isoQTL_enrich_tb) <- c("RBP", "isoQTL", "background")
RBP_isoQTL_enrich_tb <- RBP_isoQTL_enrich_tb %>% group_by(RBP) %>% mutate(pval = fisher.test(matrix(c(6269080, 238678, background, isoQTL),
                                                                            nrow = 2, byrow = TRUE))$p.value,
                                                  logOR = log(fisher.test(matrix(c(6269080, 238678, background, isoQTL),
                                                                            nrow = 2, byrow = TRUE))$estimate))
RBP_isoQTL_enrich_tb$FDR <- p.adjust(RBP_isoQTL_enrich_tb$pval, method = "BH" )
colnames(RBP_info)[1] <- "RBP"
RBP_isoQTL_enrich_tb <- left_join(RBP_isoQTL_enrich_tb, RBP_info, by = "RBP")
colSums(RBP_isoQTL_enrich_tb[which(RBP_isoQTL_enrich_tb$FDR < 0.05),c(8:26)])
write.table(RBP_isoQTL_enrich_tb[,c(1:26)], file = "/data/Choi_lung/scLongreads/Enrichment/HepG2_RBP_motif_enrichment_rst.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
RBP_info$functions <- NA
functions <- colnames(RBP_info)[c(3:21)]
for (i in 1:nrow(RBP_info)) {
  RBP_info$functions[i] <- paste(functions[(RBP_info[i,c(3:21)]==1)], sep = ",",collapse = ",")
}
RBP_info$functions <- paste(colnames(RBP_info)[(RBP_info[RBP,]==1)], sep = ",",collapse = ",")
test$functions <- NULL
test <- test %>% group_by(RBP) %>% mutate(functions = RBP_info[RBP,"functions"])
write.table(test, file = "/data/Choi_lung/scLongreads/Enrichment/HepG2_RBP_motif_enrichment.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

K562_rst <- read.table("/data/Choi_lung/scLongreads/Enrichment/K562_RBP_motif_enrichment_rst.txt",
                       sep = "\t", header = TRUE)
HepG2_rst <- read.table("/data/Choi_lung/scLongreads/Enrichment/HepG2_RBP_motif_enrichment_rst.txt",
                       sep = "\t", header = TRUE)

solid_RBP <- intersect(K562_rst$RBP[which(K562_rst$FDR < 0.05)], HepG2_rst$RBP[which(HepG2_rst$FDR < 0.05)])
tmp <- subset(K562_rst, RBP %in% solid_RBP)
colSums(tmp[,c(8:26)])

load("/data/Choi_lung/scLongreads/Enrichment/K562_enrichment_rst.RData")
# NB.nominal.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_list.rds")
all_qtl <-  readRDS(file = "/data/Choi_lung/scLongreads/mashr/all.qtl.rds")
features <- unique(all_qtl$feature_id)
ceiling(length(features)/20)
feature_chunks <- split(features, rep(1:20,ceiling(length(features)/20))[-42940])
cis.SNPs.list <- list()
# error.list <- list()
for (i in 1:20) {
  # betas1 <- subset(all_qtl, feature_id %in% feature_chunks[[i]]) %>% pivot_wider(names_from = celltype, values_from = betas,
  #                                                                                id_cols = id) %>%
  #   tibble::column_to_rownames(var = "id") %>% qDF()
  tmp <- subset(all_qtl, feature_id %in% feature_chunks[[i]]) 
  snps <- str_split_fixed(rownames(tmp), "\\|",2)[,2]
  cis.SNPs.list[[i]] <- snps
  
}

cis.SNPs <- unique(all_qtl$variant_id)
NB.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")
snp.pruned <- read.table("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/prune.prune.in")
# 2000kb
# snp.pruned <- read.table("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/prune_2kk.prune.in")

RBP_bg_snps_peaks_df_cisSNP <- subset(RBP_bg_snps_peaks_df, snp %in% cis.SNPs)
RBP_bg_snps_peaks_df_pruned_cisSNP <- subset(RBP_bg_snps_peaks_df, snp %in% snp.pruned$V1)
test_pruned <- subset(test, snp %in% snp.pruned$V1)
isoQTL_pruned <- isoQTLs[which(isoQTLs%in%snp.pruned$V1)]
functions <- RBP_info[,c(3:20)]
func_category <- colnames(functions)
enrich_tb.list <- list()
for (func in func_category) {
  func_RBPs <- rownames(functions)[(functions[[func]]==1)]
  enriched_tb <- subset(test_pruned, RBP %in% func_RBPs)
  if(nrow(enriched_tb) == 0){
    enrich_tb.list[[func]] <- NULL
  }else{
    enriched_tb$functional <- func
    enrich_tb.list[[func]] <- enriched_tb
  }
  
}
enrich_tb <- Reduce(rbind, enrich_tb.list)

bg.list <- list()
for (func in func_category) {
  func_RBPs <- rownames(functions)[(functions[[func]]==1)]
  bg_tb <- subset(RBP_bg_snps_peaks_df_pruned_cisSNP, RBP %in% func_RBPs)
  if(nrow(bg_tb) == 0){
    bg.list[[func]] <- NULL
  }else{
    bg_tb$functional <- func
    bg.list[[func]] <- bg_tb
  }
  
}
bg_tb <- Reduce(rbind, bg.list)
RBP_bg_sum <- as.data.frame(table(bg_tb$functional))
RBP_isoQTL_sum <- as.data.frame(table(enrich_tb$functional))

RBP_isoQTL_enrich_tb <- left_join(RBP_bg_sum, RBP_isoQTL_sum, by = "Var1")
colnames(RBP_isoQTL_enrich_tb) <- c("Funtional", "background", "isoQTL")
RBP_isoQTL_enrich_tb$non_isoQTL <- (RBP_isoQTL_enrich_tb$background - RBP_isoQTL_enrich_tb$isoQTL)
RBP_isoQTL_enrich_tb[is.na(RBP_isoQTL_enrich_tb)] <- 0
RBP_isoQTL_enrich_tb <- RBP_isoQTL_enrich_tb %>% group_by(Funtional) %>% mutate(pval = fisher.test(matrix(c((nrow(snp.pruned)-length(isoQTL_pruned)), length(isoQTL_pruned), non_isoQTL, isoQTL),
                                                                                                    nrow = 2, byrow = TRUE))$p.value,
                                                                          logOR = log(fisher.test(matrix(c((nrow(snp.pruned)-length(isoQTL_pruned)), length(isoQTL_pruned), non_isoQTL, isoQTL),
                                                                                                         nrow = 2, byrow = TRUE))$estimate))
RBP_isoQTL_enrich_tb$FDR <- p.adjust(RBP_isoQTL_enrich_tb$pval, method = "BH" )
write.table(RBP_isoQTL_enrich_tb, "/data/Choi_lung/scLongreads/Enrichment/K562/K562_functional_annotation_enrichment_pruned_rst.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(RBP_isoQTL_enrich_tb, "/data/Choi_lung/scLongreads/Enrichment/K562/K562_functional_annotation_enrichment_pruned_2kk_rst.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
save(list = ls(), file = "/data/Choi_lung/scLongreads/Enrichment/K562_enrichment_isoQTL_pruned_rst.RData")


isoQTL_all <- RBP_isoQTL_enrich_tb
eQTL_all <- RBP_isoQTL_enrich_tb
isoQTL_new <- RBP_isoQTL_enrich_tb
colnames(eQTL_all) <- colnames(isoQTL_all)
colnames(isoQTL_new) <- colnames(isoQTL_all)
isoQTL_all$QTL_type <- "isoQTL_all"
eQTL_all$QTL_type <- "eQTL"
isoQTL_new$QTL_type <- "isoQTL_new"
df <- do.call(rbind, list(isoQTL_all, isoQTL_new, eQTL_all))
df <- df[-which(df$logOR == -Inf),]
df$QTL_type <- factor(df$QTL_type, levels = c("isoQTL_all", "isoQTL_new", "eQTL"))
df$Funtional <- factor(df$Funtional, levels = rev(rownames))
df2 <- arrange(df, Funtional) |> distinct(Funtional)
n <- c(1:8)
df2 <- df2[2*n,]

ggplot(df, aes(y = Funtional, x = -log10(FDR), size = isoQTL, shape = QTL_type)) + 
  geom_point(aes(color = logOR), alpha = 1.0) + 
  geom_tile(data = df2,inherit.aes = FALSE,
            aes(y = Funtional, x = 0, width = Inf), fill = "gray95", height = 0.99)+
  theme(axis.text = element_text(size = 12, color = "black"), 
        axis.title  = element_text(size = 14, face = "bold"),
        legend.position='top', 
        panel.grid  = element_line(color = "grey"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = NA)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")+
  ylab("Functional annotation") +
  scale_colour_viridis_c()
ggsave("/data/lib14/project/scLongread/Fig5D_V2.pdf",width = 6,height = 7)
## get index of geom rect layer 
rect_ind <- which(lapply(g$layers, function(x) class(x$geom)[1]) == "GeomTile")

## manually change the order to put geom rect first
g$layers <- c(g$layers[rect_ind], g$layers[-rect_ind])
g


df_sum_details <- as.matrix.data.frame(table(enrich_tb$functional, enrich_tb$RBP))
colnames(df_sum_details) <- names(table(enrich_tb$RBP))
rownames(df_sum_details) <- names(table(enrich_tb$functional))
library(pheatmap)
library(ComplexHeatmap)
p <- ComplexHeatmap::pheatmap(df_sum_details, scale = "none", 
                         cluster_rows = TRUE, cluster_cols = TRUE, border_color = "black",
                         color = colorRampPalette(c("white", "yellow", brewer.pal(7,"RdYlBu")[1]))(100), 
                         number_format = "%.0f",number_color = "white",
                         display_numbers = T, fontsize_number = 9, 
                         column_names_side = c("top"), angle_col = "45")

pdf("/data/lib14/project/scLongread/FigS7A.pdf", width = 10)
print(p)
dev.off()
# ggsave("/data/lib14/project/scLongread/FigS7A.pdf", p,width = 12,height = 7)
rownames <- rownames(df_sum_details)[pheatmap::pheatmap(df_sum_details)$tree_row$order]
rownames(RBP_isoQTL_enrich_tb) <- RBP_isoQTL_enrich_tb$Funtional
RBP_isoQTL_enrich_tb <- RBP_isoQTL_enrich_tb[rownames,]
RBP_isoQTL_enrich_tb$Funtional <- factor(RBP_isoQTL_enrich_tb$Funtional, levels = rownames)

# Draw enrichment bubble plot
library(ggplot2)
library(forcats)

ggplot(RBP_isoQTL_enrich_tb, aes(y = reorder(Funtional, rev(as.numeric(Funtional))), x = logOR, size = isoQTL)) + 
  geom_point(aes(color = -log10(FDR)), alpha = 1.0) + theme(legend.position='bottom') +
  theme_bw(base_size = 12)+ylab("Functional annotation") +
  scale_colour_viridis_c()

ggplot(RBP_isoQTL_enrich_tb, aes(y = reorder(Funtional, rev(as.numeric(Funtional))), x = logOR, size = isoQTL)) + 
  geom_point(aes(color = -log10(FDR)), alpha = 1.0) + 
  theme(axis.text = element_text(size = 12, color = "black"), 
        axis.title  = element_text(size = 14, face = "bold"),
        legend.position='top', 
        panel.grid  = element_line(color = "grey"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")) +
  ylab("Functional annotation") +
  scale_colour_viridis_c()
ggsave("/data/lib14/project/scLongread/Fig5D.pdf",width = 6,height = 7)
write.table(RBP_isoQTL_enrich_tb, file = "/data/Choi_lung/scLongreads/Enrichment/K562/K562_functional_annotation_enrichment_rst.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


Nominal_rst <- do.call(rbind,NB.nominal.sig.list)
Nominal_rst <- subset(Nominal_rst, snp %in% enrich_tb$snp)
Nominal_rst_sub <- Nominal_rst[!duplicated(Nominal_rst$snp),]
rownames(Nominal_rst_sub) <- Nominal_rst_sub$snp
enrich_tb <- enrich_tb %>% group_by(snp) %>% mutate(eIsoform = Nominal_rst_sub[snp,]$phenotype_id,
                                                    transcript_name = Search_transcript_name2(Nominal_rst_sub[snp,]$phenotype_id),
                                                    Celltype = Nominal_rst_sub[snp,]$Celltype
                                                    )
write.table(enrich_tb, file = "/data/Choi_lung/scLongreads/Enrichment/K562/LDpruned_isoQTL_located_RBP_motif.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

NB.cis.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_list.rds")
tested_phenotypes <- lapply(NB.cis.list, function(x) x$phenotype_id)
tested_phenotypes <- Reduce(c, tested_phenotypes)
tested_phenotypes <- unique(tested_phenotypes)
tmp <- subset(TALON_afterqc_orf_secondpass2, annot_transcript_id %in% tested_phenotypes)
table(tmp$subcategory)
