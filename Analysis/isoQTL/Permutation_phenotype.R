# Generate phenotype files for permutation

# Bolun Li
# Jun 30 2024


.libPaths()
library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)
library(future)
library(stringr)
library(stringi)
library(forcats)
library(Matrix)


# load transcript metadata generated from TALON and integrated with TranDecoder results
TALON_afterqc_orf_secondpass <- readRDS("/data/Choi_lung/scLongreads/Seurat/TALON_afterqc_orf_secondpass_w_known.rds")
# load the final version with cell type annotation
lr <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final_v2_wRNAexpr.RDS")

# load the original isoform level expression profiles
lr.isoform <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/final/lr.sc.data.Isoform.22samples.combined.rds")
ncol(lr)
lr.isoform.sub <- subset(lr.isoform, cell = colnames(lr))
lr.isoform.sub$Sample <- lr$Sample
lr.isoform.sub$Sample_NCI <- lr$Sample_NCI
lr.isoform.sub$Celltype <- lr$Celltype
# prepare splicing QTL analysis
gff = as.data.frame(rtracklayer::import('/data/Choi_lung/scLongreads/TALON_workspace/test10s/NCI_lung_secondpass.isoforms_w_known.gtf'))
table(gff$type)
transcript_gff <- subset(gff, type == "transcript")
rm(gff)
pos <- data.frame(X = transcript_gff$transcript_id,
                  chr = transcript_gff$seqnames, 
                  start = transcript_gff$start, 
                  end = transcript_gff$start,
                  gene_id = transcript_gff$gene_id)
tss <- read.table("/data/Choi_lung/scLongreads/TALON_workspace/test10s/tss.bed", sep = "\t", header = FALSE)
rownames(tss) <- tss$V4
rownames(pos) <- pos$X

library("tidyverse")
`%nin%` <- Negate(`%in%`)

celltypes <- levels(lr.isoform.sub$Celltype)
number_iso <- c()
tensor_path <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level"

celltypes <- celltypes[1:8]
celltypes <- celltypes[-2]
for(celltype in celltypes){
  # pilot analysis
  seur <- subset(lr.isoform.sub, Celltype == celltype)
  celltype_for_path <- gsub(" ", "_", celltype)
  cell_number <- table(seur$Sample)
  ind_filtered <- names(cell_number)[which(cell_number < 5)]
  table(seur$Sample)
  seur <- subset(seur, Sample %nin% ind_filtered)
  count_mtx <- seur@assays$Isoform@counts
  count_mtx <- count_mtx[TALON_afterqc_orf_secondpass$transcript_name_unique,]
  count_mtx@Dimnames[[1]] <- TALON_afterqc_orf_secondpass$annot_transcript_id
  seur$Sample <- as.character(seur$Sample)
  group <- seur$Sample %>% fct_inorder()
  group_mat <- sparse.model.matrix(~ 0 + group) %>% t
  # Adjust row names to get the correct final row names
  rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
  count_mtx_T <- Matrix::t(count_mtx)
  count_mtx_sum <- group_mat %*% count_mtx_T  
  dim(count_mtx_sum)
  sum(count_mtx_sum)
  count_mtx_sum <- Matrix::t(count_mtx_sum)
  dim(count_mtx_sum)
  pos <- data.frame(X = transcript_gff$transcript_id,
                    chr = transcript_gff$seqnames, 
                    start = transcript_gff$start, 
                    end = transcript_gff$start,
                    gene_id = transcript_gff$gene_id)
  length(unique(pos$X))
  
  rownames(pos) <- pos$X
  pos <- pos[TALON_afterqc_orf_secondpass$annot_transcript_id,]
  gene_group <- pos$gene_id %>% fct_inorder()
  gene_group_mat <- sparse.model.matrix(~ 0 + gene_group) %>% t
  # Adjust row names to get the correct final row names
  rownames(gene_group_mat) <- rownames(gene_group_mat) %>% str_extract("(?<=^gene_group).+")
  dim(gene_group_mat)
  count_mtx_gene <- gene_group_mat %*% count_mtx_sum
  sum(count_mtx_gene)
  dim(count_mtx_gene)
  N <- 129-length(ind_filtered)
  message(celltype_for_path, ": N = ", N)
  idx <- order(pos$gene_id)
  pos <- pos[idx,]
  non_zero_mtx <- ((count_mtx_sum[idx,]) != 0)
  dim(non_zero_mtx)
  non_zero_transcripts <- rowSums(non_zero_mtx)
  idx_trancripts_.1 <- which(non_zero_transcripts > N*0.2)
  count_mtx_sum = NormalizeData(count_mtx_sum[idx,])
  
  centered_mtx <- scale(t(count_mtx_sum[idx_trancripts_.1,]))
  dim(centered_mtx)
  which(centered_mtx=="NaN")
  library(limma)
  quantile.norm.mtx <- normalizeQuantiles(t(centered_mtx))
  dim(quantile.norm.mtx)
  samples <- colnames(count_mtx_sum)
  transcripts <- rownames(count_mtx_sum)
  transcripts <- transcripts[idx_trancripts_.1]
  
  transcripts_ids <- rownames(tss)[which(rownames(tss) %in% transcripts)]
  colnames(tss) <- c("#chr", "start", "end", "trascript_id")
  tss_sub <- tss[transcripts_ids,c(1:4)]
  quantile.norm.mtx_permutated <- apply(quantile.norm.mtx, 1, function(x){
    seed <- sample(100000, 1)
    message(seed)
    set.seed(seed)
    idxs <- sample(N, N, replace = FALSE)
    x <- x[idxs]
  })
  tmp1 <- sapply(1:nrow(quantile.norm.mtx), function(i){
    cor((quantile.norm.mtx[i,]), (quantile.norm.mtx_permutated[,i]))
  })
  message("Median: ", median(tmp1), "; Mean: ", mean(tmp1))
  rownames(quantile.norm.mtx_permutated) <- colnames(count_mtx_sum)
  # bed <- cbind(tss_sub,quantile.norm.mtx[rownames(tss_sub),])
  bed <- cbind(tss_sub,t(quantile.norm.mtx_permutated)[rownames(tss_sub),]) # permutation
  dir.create(file.path(tensor_path, celltype_for_path), showWarnings = FALSE)
  # write.table(samples, file = paste0(file.path(tensor_path, celltype_for_path), "/samples.txt"), quote = F, row.names = F, col.names = F)
  write.table(bed, file = paste0(file.path(tensor_path, celltype_for_path), "/phenotype_expr_permutated.bed"), sep='\t', quote=F, row.names=FALSE, col.names = T, eol = '\n')
  # peer <- t(quantile.norm.mtx[rownames(tss_sub),])
  peer <- t(t(quantile.norm.mtx_permutated)[rownames(tss_sub),]) # permutation
  write.table(peer, file = paste0(file.path(tensor_path, celltype_for_path), "/peer_expr_permutated.csv"), sep = ',', quote = F, row.names = F, col.names = F)
}



