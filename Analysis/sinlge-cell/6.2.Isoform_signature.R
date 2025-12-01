lr.list <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr.list.byCT.rds")

seur <- lr.list$AT2
seur <- lr.list$AT1
rst.list <- list()
celltypes <- levels(lr$Celltype)
for (i in 1:9) {
  celltype <- celltypes[i]
  seur <- lr.list[[celltype]]
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
  idx <- order(pos$gene_id)
  pos <- pos[idx,]
  transcript_number <- table(pos$gene_id)
  count_mtx_gene_rep <- apply(as.matrix(count_mtx_gene)[unique(pos$gene_id),], 2, function(x){
    rep(x, transcript_number)
  })
  
  non_zero_mtx <- ((count_mtx_sum[idx,]) != 0)
  dim(non_zero_mtx)
  non_zero_transcripts <- rowSums(non_zero_mtx)
  count_mtx_per <- as.matrix(count_mtx_sum)[idx,]/count_mtx_gene_rep
  count_mtx_per[count_mtx_per == "NaN"] <- 0
  tmp <- count_mtx_per[grep("ENSG00000204305",rownames(count_mtx_gene_rep)),]
  library(reshape2)
  tmp_long <- melt(tmp)
  tmp_long$Celltype <- celltype
  rst.list[[i]] <- tmp_long
  
}
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
idx <- order(pos$gene_id)
pos <- pos[idx,]
transcript_number <- table(pos$gene_id)
count_mtx_gene_rep <- apply(as.matrix(count_mtx_gene)[unique(pos$gene_id),], 2, function(x){
  rep(x, transcript_number)
})

non_zero_mtx <- ((count_mtx_sum[idx,]) != 0)
dim(non_zero_mtx)
non_zero_transcripts <- rowSums(non_zero_mtx)
count_mtx_per <- as.matrix(count_mtx_sum)[idx,]/count_mtx_gene_rep
idx_trancripts_.1 <- which(non_zero_transcripts > N*0.3)
number_iso <- c(number_iso, length(idx_trancripts_.1))
count_mtx_per[count_mtx_per == "NaN"] <- 0
tmp <- count_mtx_per[grep("ENSG00000168484",rownames(count_mtx_gene_rep)),]
tmp <- count_mtx_per[grep("ENSG00000000003",rownames(count_mtx_gene_rep)),]
library(reshape2)
tmp_long1 <- melt(tmp)
tmp_long1$Celltype <- "AT1"
tmp_long <- melt(tmp)
tmp_long$Celltype <- "AT2"
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
df <- Reduce(rbind, rst.list)
ggplot(data=df, aes(x=value, group=Var1, fill=Var1)) +
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()+facet_wrap(~Celltype)
ggplot(data=df, aes(x=value, group=Celltype, fill=Celltype)) +
  geom_density(adjust=1.5, alpha=.4) +
  facet_wrap( ~ Var1, scales="free_x")
ggsave("/data/lib14/project/scLongread/ENSG00000204305.pdf", width = 20, height = 20)
ggplot(data=subset(tmp_long, Var1 %in% tmp_long$Var1[1:10]), aes(x=value, group=Var1, fill=Var1)) +
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()+NoLegend()

