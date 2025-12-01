# generate colocalization heatmap
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
library(arrow)
transcripts <- c("HLA-DRB5-TALONT001746251","SECISBP2L-TALONT001050043",
                 "HLA-DRB1-TALONT001751909","CHMP3-203","RSPH4A-TALONT003150695",
                 "PPIL6-207","HLA-A-TALONT001165898","PSMA4-201","HSPA1B-201")
celltypes <- c("AT2", "Multiciliated", "Alveolar_transitional_cells",
               "NK_cells", "Non_classical_monocytes", "Lymphatic_EC")
NB.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")
NB.sig.list <- lapply(names(NB.sig.list), function(x){
  rst <- NB.sig.list[[x]]
  rst$Celltype <- x
  return(rst)
})
isoQTL_sig <- do.call(rbind,NB.sig.list)

isoQTL_sig$sig_iso_ct <- paste(isoQTL_sig$phenotype_id, isoQTL_sig$Celltype, sep = "-")
Byun_total <- read.table("/data/Choi_lung/scLongreads/colocalization/Coloc_Byun_total_lung_PP_sum.txt", header = TRUE)
Byun_luad <- read.table("/data/Choi_lung/scLongreads/colocalization/Byun_coloc/ADE_sig_coloc_table.txt", header = TRUE)
Byun_sqc <- read.table("/data/Choi_lung/scLongreads/colocalization/Byun_coloc/SQC_sig_coloc_table.txt", header = TRUE)
Shi_EAS <- read.table("/data/Choi_lung/scLongreads/colocalization/Coloc_Shi_EAS_PP_sum.txt", header = TRUE)
Shi_meta <- read.table("/data/Choi_lung/scLongreads/colocalization/Coloc_Shi_meta_PP_sum.txt", header = TRUE)
coloc.list <- list(Shi_EAS=Shi_EAS, Shi_meta=Shi_meta, Byun_lung=Byun_total,Byun_luad= Byun_luad[,c(1:9)],Byun_sqc= Byun_sqc[,c(1:9)])

coloc.list <- lapply(1:5,function(i){
  coloc_rst <- coloc.list[[i]]
  coloc_rst$GWAS <- names(coloc.list)[i]
  return(coloc_rst)
})

coloc_rst <- do.call(rbind, coloc.list)
coloc_rst <- subset(coloc_rst, PP.H4.abf > 0.7)

output.fev1 <- readRDS("/data/Choi_lung/scLongreads/colocalization/lung_function/Colocalization_output_table_fev1.rds")
output.ratio <- readRDS("/data/Choi_lung/scLongreads/colocalization/lung_function/Colocalization_output_table.rds")
output.ratio$GWAS <- "FEV1_FVC"
output.fev1$GWAS <- "FEV1"
output.fev1 <- output.fev1[!duplicated(output.fev1$sig_iso_ct),]
output.ratio <- output.ratio[!duplicated(output.ratio$sig_iso_ct),]
coloc_rst$GWAS <- "Lung_cancer"
coloc_rst$isoform_ct <- paste(coloc_rst$transcript_id, coloc_rst$Celltype,sep = "_")
coloc_rst <- coloc_rst[!duplicated(coloc_rst$isoform_ct),]
tmp <- rbind(coloc_rst[,c(1:10)], output.ratio[,c(1:9,12)])
coloc_final <- rbind(tmp, output.fev1[,c(1:9,12)])

coloc_rst_PPH4 <- as.matrix.data.frame(table(coloc_final$Celltype, coloc_final$transcript_name))
rownames <- names(table(coloc_final$Celltype))
rownames(coloc_rst_PPH4) <- rownames
colnames <- names(table(coloc_final$transcript_name))
colnames(coloc_rst_PPH4) <- colnames

transcripts <- unique(coloc_final$transcript_name)
celltypes <- unique(coloc_final$Celltype)
celltype_levels <- levels(lr$Celltype)
celltypes <- intersect(gsub(" ", "_", celltype_levels), celltypes)
coloc_rst_PPH4 <- coloc_rst_PPH4[celltypes,]
p <- pheatmap::pheatmap(coloc_rst_PPH4, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
                   color = c(NA,"yellow", "orange","red"), border_color = NA,cellwidth = 15, cellheight = 15)
ggsave("/data/lib14/project/scLongread/ASHG_coloc_count.pdf", p,width = 15, height = 7)
count_mtx <- lr.isoform.sub@assays$Isoform@counts
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass$transcript_name_unique,]
TALON_afterqc_orf_secondpass$total_counts <- Matrix::rowSums(count_mtx)

TALON_afterqc_orf_secondpass <- TALON_afterqc_orf_secondpass %>% group_by(annot_gene_id) %>%
  mutate(percentage = total_counts/sum(total_counts))
rownames(TALON_afterqc_orf_secondpass) <- TALON_afterqc_orf_secondpass$transcript_name_unique

hm <- matrix(t(TALON_afterqc_orf_secondpass2[colnames,"percentage"]), nrow = 1)
colnames(hm) <- colnames
rownames(hm) <- "Percentage"
p <- pheatmap::pheatmap(hm, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
                        cellwidth = 15, cellheight = 15)
ggsave("/data/lib14/project/scLongread/ASHG_coloc_expr_percentage.pdf", p,width = 15, height = 3.2)

tmp <- data.frame(Freq = rowSums(coloc_rst_PPH4),
                  Celltype = celltypes)

tmp$Celltype <- factor(tmp$Celltype, levels = celltypes)
x.axis.bar.df <- tmp
x.axis.bar.df$Var1 <- factor(x.axis.bar.df$Celltype, levels = celltypes)
x.axis.bar.df$Category <- c(rep("Epithelial",8), rep("Immune", 7), rep("Endothelial", 2))
x.axis.bar.df$Category <- factor(x.axis.bar.df$Category , levels = c("Epithelial", "Immune","Endothelial"))
bp.x <- ggplot(data = x.axis.bar.df, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity", aes(fill = Category)) + theme_classic() +
  theme(axis.title.y = element_blank(), axis.text.y =  element_text(size = 12), 
        axis.ticks.y = element_blank(), axis.text.x = element_text(size = 15), 
        axis.title.x = element_text(size = 20, margin = margin(10,0,0,0)),
        legend.position = "none") +
  scale_fill_manual(values = c("#ffd700","#fa8775","#cd34b5"))
bp.x
bp.x.clean <- bp.x +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x  = element_blank(),
        legend.position="none")
bp.x.clean
ggsave("/data/lib14/project/scLongread/ASHG_coloc_count_bar.pdf", bp.x.clean, width = 7, height = 5)

coloc_final$ID_celltype <- paste(coloc_final$transcript_name, coloc_final$Celltype, sep = "_")
DEI_sig$ID_celltype <- paste(DEI_sig$X, DEI_sig$Celltype, sep = "_")
unique(coloc_final$transcript_id[which(coloc_final$ID_celltype %in% DEI_sig$ID_celltype)])
unique(coloc_final$transcript_id)
df_coloc_sum <- data.frame(transcript_id = unique(coloc_final$transcript_id),
                           transcript_name = unique(coloc_final$transcript_name),
                           percentage = hm[1,]
                           )
df_coloc_sum$abundant <- (df_coloc_sum$percentage > 0.144)
df_coloc_sum$gene_name <- NA
tmp_gene_names <- (TALON_afterqc_orf_secondpass2[df_coloc_sum$transcript_name,"annot_gene_name"])
df_coloc_sum[["gene_name"]] <- tmp_gene_names$annot_gene_name
df_coloc_sum$gene_name[39] <- "HLA-DRB1+HLA-DRB5"
df_coloc_sum$gene_name[15] <- "HLA-DQA1+HLA-DQA2"
df_coloc_sum$replicated <- (df_coloc_sum$gene_name %in% c("SECISBP2L", "RBM23", "MUC1","DBNDD2",
                                                          "MRPL21", "MAP1LC3A", "HLA-DRB5","CD55"))
df_coloc_sum$eQTL_coloc <- (df_coloc_sum$transcript_id %in% unique(eQTL_coloc_sub1$phenotype_id.x))
df_coloc_sum$ct_spec <- (df_coloc_sum$transcript_id %in% ct_spec_isoforms)
df_coloc_sum$novel <- grepl("TALON",df_coloc_sum$transcript_id)


lung_func_gwas <- read.table("Need_clean/Lung_func_gwas_association.txt", header = TRUE, sep = "\t")
lung_func_gwas_fev1 <- read.table("Need_clean/Lung_func_fev1_gwas_association.txt", header = TRUE, sep = "\t")
lung_func_gwas <- lung_func_gwas[,c("STRONGEST.SNP.RISK.ALLELE", "REGION")]
lung_func_gwas <- lung_func_gwas[!duplicated(lung_func_gwas),]
colnames(lung_func_gwas)[1] <- "leadSNP"
lung_func_gwas_fev1 <- lung_func_gwas_fev1[,c("STRONGEST.SNP.RISK.ALLELE", "REGION")]
lung_func_gwas_fev1 <- lung_func_gwas_fev1[!duplicated(lung_func_gwas_fev1),]
colnames(lung_func_gwas_fev1)[1] <- "leadSNP"
output.fev1 <- left_join(output.fev1, lung_func_gwas_fev1, by = "leadSNP")
output.ratio <- left_join(output.ratio, lung_func_gwas, by = "leadSNP")
tmppp1 <- output.fev1[,c("transcript_id", "REGION")]
tmppp1 <- tmppp1[!duplicated(tmppp1),]
tmppp2 <- output.ratio[,c("transcript_id", "REGION")]
tmppp2 <- tmppp2[!duplicated(tmppp2),]
tmppp <- rbind(tmppp2,tmppp1)
tmppp <- tmppp[!duplicated(tmppp),]
df_coloc_sum <- left_join(df_coloc_sum, tmppp, by = "transcript_id")
df_coloc_sum <- df_coloc_sum[!duplicated(df_coloc_sum$transcript_id),]
df_coloc_sum$REGION[2:8] <- c("2p11.2","6q22.1","6p21.33","6p21.32","6p21.33","6q21","6p21.33")
df_coloc_sum$eQTL_coloc[which(!(df_coloc_sum$transcript_id %in% eQTL_coloc$phenotype_id.x))] <- "No_matched"
table(df_coloc_sum$eQTL_coloc)
table(df_coloc_sum$eQTL_coloc, df_coloc_sum$replicated)
df_coloc_sum$eQTL_coloc_replicated <- paste(df_coloc_sum$eQTL_coloc, df_coloc_sum$replicated, sep = "-")
table(df_coloc_sum$eQTL_coloc_replicated, df_coloc_sum$abundant)



TALON_afterqc_orf_secondpass <- TALON_afterqc_orf_secondpass %>% group_by(annot_gene_id) %>%
  mutate(rank = (rank(-total_counts))/length(total_counts))
rownames(TALON_afterqc_orf_secondpass) <- TALON_afterqc_orf_secondpass$transcript_name_unique
df_coloc_sum$rank <- TALON_afterqc_orf_secondpass[df_coloc_sum$transcript_name,]$rank


# TWAS
colnames(hm) <- hm["transcript_name_unique",]
hm <- (TALON_afterqc_orf_secondpass2[unique(TWAS_tb$transcript_name),"percentage"])
names(hm) <- unique(TWAS_tb$transcript_name)
df_twas_sum <- data.frame(transcript_id = unique(TWAS_tb$ID),
                           transcript_name = unique(TWAS_tb$transcript_name),
                           percentage = hm$percentage
)
df_twas_sum$abundant <- (df_twas_sum$percentage > ninetyth_percentile)

tmp_gene_names <- (TALON_afterqc_orf_secondpass[df_twas_sum$transcript_name,"annot_gene_name"])

df_twas_sum$gene_name <-  tmp_gene_names

rownames(df_twas_sum) <- c(1:74)
df_twas_sum$gene_name[15] <- "HLA-DRB1+HLA-DRB5"
df_twas_sum$gene_name[24] <- "HLA-DRB1+HLA-DRB5"
df_twas_sum$gene_name[11] <- "HLA-DQA1+HLA-DQA2"
df_twas_sum$gene_name[59] <- "ATP5MG+KMT2A"
TWAS_previous <- readxl::read_xlsx("TWAS_sum_previous.xlsx")

TWAS_previous$gene_name <- str_split_fixed(TWAS_previous$Phenptype_id, ":", 2)[,1]
TWAS_previous_iso <- intersect(unique(TWAS_previous$gene_name), unique(df_twas_sum$gene_name))
df_twas_sum$replicated <- (df_twas_sum$gene_name %in% TWAS_previous_iso)
eQTL_coloc_sub1 <- subset(eQTL_coloc, PP.H4.abf > 0.7)
df_twas_sum$eQTL_coloc <- (df_twas_sum$transcript_id %in% unique(eQTL_coloc_sub1$phenotype_id.x))
sum(df_twas_sum$eQTL_coloc)
df_twas_sum$ct_spec <- (df_twas_sum$transcript_id %in% ct_spec_isoforms)
df_twas_sum$novel <- grepl("TALON",df_twas_sum$transcript_id)
tmppppp <- TWAS_tb[,c("ID", "Locus")]
tmppppp <- tmppppp[!duplicated(tmppppp),]
colnames(tmppppp)[1]  <- "transcript_id"
df_twas_sum <- left_join(df_twas_sum, tmppppp, by = "transcript_id")
df_twas_sum <- df_twas_sum[!duplicated(df_twas_sum$transcript_id),]
df_twas_sum$Locus_novel <- "Replicated"
df_twas_sum$Locus_novel[which(df_twas_sum$Locus == "")] <- "New"

Novel_sum <- read.table("Novel_isoform_sum.txt", sep = "\t", header = TRUE)

Novel_sum$pct <- Novel_sum$Novel.isoforms/Novel_sum$Total.number
Novel_sum$pct <- round(Novel_sum$pct * 100, 2)
Novel_sum_long <- melt(Novel_sum[,c(1:3)])

ggplot(Novel_sum_long, aes(x=variable, y=value, fill = variable)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("")+ylab("Number of total/novel isoforms") +
  theme(axis.text.x = element_text(color = "black", size = 12,face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        strip.background = element_blank(), 
        legend.justification="right")+facet_wrap(~Study)
ggsave("LR_studies_novel_sum.pdf", width = 10, height = 7)


TWAS_tb$match_CT_id <- paste(TWAS_tb$transcript_name, TWAS_tb$Celltype, sep = ":")
DEI_sig$Celltype <- str_split_fixed(rownames(DEI_sig), "\\.", 2)[,1]
DEI_sig$match_CT_id <- paste(DEI_sig$X, DEI_sig$Celltype, sep = ":")

ct_specific_TWAS <- TWAS_tb$match_CT_id[which(TWAS_tb$match_CT_id %in% DEI_sig$match_CT_id)]
ct_specific_TWAS <- unique(str_split_fixed(ct_specific_TWAS, ":", 2)[,1])
df_twas_sum$ct_spec_DEI <- (df_twas_sum$transcript_name %in% ct_specific_TWAS)
df_twas_sum$eIsoform <- (df_twas_sum$transcript_id %in% unique(final_table$phenotype_id))
table(df_twas_sum$eIsoform, df_twas_sum$ct_spec_DEI)

twas_isoform_nonHLA_loci <- df_twas_sum$transcript_id[!grepl("6p21", df_twas_sum$Locus)]
TWAS_tb_nonHLA_loci <- subset(TWAS_tb, ID %in% twas_isoform_nonHLA_loci)
table(TWAS_tb_nonHLA_loci$Celltype) # 11 in immune, 9 in epithelial, 5 in endothelial, 3 in stromal


count_mtx <- readRDS("~/Desktop/scLongread/sc_files/count_mtx.rds")
metadata <- readRDS("~/Desktop/scLongread/sc_files/metadata.rds")

count_mtx <- count_mtx[TALON_afterqc_orf_secondpass2$transcript_name_unique,]
metadata$orig.ident <- as.character(metadata$orig.ident)
group <- metadata$orig.ident %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx)
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
sum(count_mtx_sum)
count_mtx_sum <- Matrix::t(count_mtx_sum)
dim(count_mtx_sum)
batch_level_expr_mtx <- count_mtx_sum

# Individual level
count_mtx <- readRDS("~/Desktop/scLongread/sc_files/count_mtx.rds")
metadata <- readRDS("~/Desktop/scLongread/sc_files/metadata.rds")

count_mtx <- count_mtx[TALON_afterqc_orf_secondpass2$transcript_name_unique,]
table(metadata$Sample)
metadata$Sample <- as.character(metadata$Sample)
group <- metadata$Sample %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx)
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
sum(count_mtx_sum)
count_mtx_sum <- Matrix::t(count_mtx_sum)
dim(count_mtx_sum)
Idv_level_expr_mtx <- count_mtx_sum

DEI_sig <- do.call(rbind, DEI_sig.list)
DEG_sig <- do.call(rbind, DEG_sig.list)
DEI_sig$gene <- TALON_afterqc_orf_secondpass[DEI_sig$X,"annot_gene_name"]
length(unique(DEI_sig$gene)) # 
HLCA_markers <- c(Epi_sig_list, Endo_sig_list, Imm_sig_list, Stroma_sig_list)
length(which(HLCA_markers %in% unique(DEI_sig$gene)))
length(unique(DEI_sig$X))
length(unique(DEG_sig$X))
length(which(unique(DEI_sig$gene) %in% unique(DEG_sig$X)))

TALON_afterqc_orf_secondpass2$Total_count <- rowSums(count_mtx)
TALON_afterqc_orf_secondpass2 <- TALON_afterqc_orf_secondpass2 %>% 
  group_by(annot_gene_id) %>%
  mutate(Abundance_rank = rank(-Total_count))

most_abundant_isoforms <- TALON_afterqc_orf_secondpass2$transcript_name_unique[which(TALON_afterqc_orf_secondpass2$Abundance_rank == 1)]

length(which(most_abundant_isoforms %in% unique(DEI_sig$X)))

most_abundant_DEI <- most_abundant_isoforms[which(most_abundant_isoforms %in% unique(DEI_sig$X))]
names(most_abundant_DEI) <- TALON_afterqc_orf_secondpass[most_abundant_DEI,"annot_gene_name"]
length(most_abundant_DEI[intersect(HLCA_markers, names(most_abundant_DEI) )])

# lead sQTLs
GTEx_sGene <- read.table(gzfile("../../Downloads/GTEx_Analysis_v8_sQTL/Lung.v8.sgenes.txt.gz"),sep="\t", header = TRUE)
GTEx_sGene$gene_id_new <- str_split_fixed(GTEx_sGene$gene_id, "\\.", 2)[,1]
set3 <- unique(GTEx_sGene$gene_id_new[which(GTEx_sGene$qval < 0.05)])
GTEx_sGene_sig <- GTEx_sGene[which(GTEx_sGene$qval < 0.05),]
final_table$If_sGenes <- FALSE
final_table$If_sGenes[which(final_table$gene_id %in% set3)] <- TRUE
length(unique(final_table$gene_id[final_table$If_sGenes]))

AF_diff_sum <- final_table[,c("gene_id", "variant_id", "If_sGenes")]
AF_diff_sum <- AF_diff_sum[!duplicated(AF_diff_sum),]
colnames(AF_diff_sum)[2] <- "rsid"
AF_diff_final <- left_join(AF_diff, AF_diff_sum, by = "rsid")
length(unique(AF_diff_final$rsid[AF_diff_final$If_sGenes]))
length(unique(AF_diff_final$rsid[!AF_diff_final$If_sGenes]))
length(unique(AF_diff_new$rsid[AF_diff_new$If_sGene & AF_diff_new$ALT_FREQS.EUR < 0.01 & AF_diff_new$ALT_FREQS.EAS > 0.01]))

variants_EAS_only <- unique(AF_diff_new$rsid[AF_diff_new$If_sGene & AF_diff_new$ALT_FREQS.EUR < 0.01 & AF_diff_new$ALT_FREQS.EAS > 0.01])

AF_diff_new <- read.table("Need_clean/AF_isoQTL_1KGenome_updated.txt", header = TRUE, sep = "\t")
AF_diff_new_sGenes <- subset(AF_diff_new, If_sGene == TRUE)
colnames(AF_diff_new_sGenes)[13] <- "gene_id_new"
AF_diff_new_sGenes_sQTL <- left_join(AF_diff_new_sGenes, GTEx_sGene_sig, by = "gene_id_new")

mytoken = "149dc339070c"

LDpair("rs61818036", 
       "rs3118122", 
       pop = "EAS", 
       token = mytoken, 
       output = "table", 
       genome_build = "grch38"
)

tb <- LDpair("rs61818036", 
       "rs3118122", 
       pop = "EUR", 
       token = mytoken, 
       output = "table", 
       genome_build = "grch38"
)

tb.list <- lapply(c(1:1334), function(i){
  tryCatch({
    tb <- LDpair(AF_diff_new_sGenes_sQTL$rsid[i], 
                 AF_diff_new_sGenes_sQTL$rs_id_dbSNP151_GRCh38p7[i], 
                 pop = "EAS", 
                 token = mytoken, 
                 output = "table", 
                 genome_build = "grch38"
    )
  }, error=function(e){})
  
})

tb.list.EUR <- lapply(c(1:1334), function(i){
  tryCatch({
    tb <- LDpair(AF_diff_new_sGenes_sQTL$rsid[i], 
                 AF_diff_new_sGenes_sQTL$rs_id_dbSNP151_GRCh38p7[i], 
                 pop = "EUR", 
                 token = mytoken, 
                 output = "table", 
                 genome_build = "grch38"
    )
  }, error=function(e){})
  
})
df_tmp <- as.data.frame(matrix(as.character(NA), nrow = 1, ncol = 18))

colnames(df_tmp) <- colnames(tb.list[[1]])
tb.list.tmp <- lapply(tb.list, function(x){
  if(is.null(x)){return(df_tmp)}else{return(x[,c(1:18)])}
})
library(data.table)
tb_EAS <- rbindlist(tb.list, fill = TRUE)
tb_EAS <- do.call(rbind, tb.list.tmp)
sum(grepl("linkage equilibrium", tb_EAS$ld))
tb_EAS_only <- subset(tb_EAS, var1 %in% variants_EAS_only)
sum(grepl("linkage equilibrium", tb_EAS_only$ld))



tb_EUR <- rbindlist(tb.list.EUR, fill = TRUE)
sum(grepl("linkage equilibrium", tb_EUR$ld))

# Figure 2F
length(which(df_niso_per_gene$total_n == 1 & df_niso_per_gene$known_n == 1))
length(which(df_niso_per_gene$total_n == 1))
df_niso_per_gene$novel_n[is.na(df_niso_per_gene$novel_n)]

load("/Users/lib14/Downloads/Additional_seq_data.RData")

Idv_level_if_expr <- Idv_level_expr_mtx > 0 
batch_level_if_expr <- batch_level_expr_mtx > 0 
TALON_afterqc_orf_secondpass2$batch_level <- rowSums(batch_level_if_expr)
TALON_afterqc_orf_secondpass2$Indv_level <- rowSums(Idv_level_if_expr)
NCI4_4SMRT_isoforms <- (unique(intersect(iso.info[[1]]$annot_transcript_id, TALON_afterqc_orf_secondpass2$annot_transcript_id))) # 83512
NCI4_2SMRT_isoforms <- (unique(intersect(iso.info1[[1]]$annot_transcript_id, TALON_afterqc_orf_secondpass2$annot_transcript_id))) # 73880
NCI4_added_isoforms <- setdiff(NCI4_4SMRT_isoforms, NCI4_2SMRT_isoforms)
NCI4_added_isoforms_info <- subset(TALON_afterqc_orf_secondpass2, annot_transcript_id %in% NCI4_added_isoforms)
NCI4_added_isoforms_info <- subset(NCI4_added_isoforms_info, batch_level != 0)
hist(NCI4_added_isoforms_info$batch_level, breaks = 20)
hist(NCI4_added_isoforms_info$Indv_level, breaks = 200)
sum(NCI4_added_isoforms_info$Indv_level == 1)
sum(NCI4_added_isoforms_info$batch_level == 1)/nrow(NCI4_added_isoforms_info)
NCI19_4SMRT_isoforms <- (unique(intersect(iso.info[[2]]$annot_transcript_id, TALON_afterqc_orf_secondpass2$annot_transcript_id))) # 89569
NCI19_2SMRT_isoforms <- (unique(intersect(iso.info1[[2]]$annot_transcript_id, TALON_afterqc_orf_secondpass2$annot_transcript_id))) # 78183
NCI19_added_isoforms <- setdiff(NCI19_4SMRT_isoforms, NCI19_2SMRT_isoforms)

NCI19_added_isoforms_info <- subset(TALON_afterqc_orf_secondpass2, annot_transcript_id %in% NCI19_added_isoforms)
NCI19_added_isoforms_info <- subset(NCI19_added_isoforms_info, batch_level != 0)
hist(NCI19_added_isoforms_info$batch_level, breaks = 20)
hist(NCI19_added_isoforms_info$Indv_level, breaks = 100)
sum(NCI19_added_isoforms_info$Indv_level == 1)
sum(NCI19_added_isoforms_info$batch_level == 1)
NCI19_added_isoforms_info$unique_dataset[NCI19_added_isoforms_info$batch_level == 0]
sum(NCI19_added_isoforms_info$batch_level == 1)/nrow(NCI19_added_isoforms_info)
library(ggplot2)
NCI19_added_isoforms_info$Novelty <- "Novel"
NCI19_added_isoforms_info$Novelty[which(NCI19_added_isoforms_info$transcript_catalog == "Known")] <- "Known"
NCI19_added_isoforms_df <- as.data.frame(table(NCI19_added_isoforms_info$Novelty, NCI19_added_isoforms_info$batch_level))
NCI19_added_isoforms_df$pct <- NCI19_added_isoforms_df$Freq/sum(NCI19_added_isoforms_df$Freq)
NCI19_added_isoforms_df$Var1 <- factor(NCI19_added_isoforms_df$Var1, levels = c("Novel","Known"))
n = 2
cols = gg_color_hue(n)
cat.palette = c( "Known"="#6BAED6", 
                 "Novel"="#F8766D")
ggplot(NCI19_added_isoforms_df) +
  geom_bar(aes(x=Var2, y=pct, fill=Var1), position="stack", stat="identity")+
  scale_fill_manual(values = cat.palette)+ ylab("Fraction")+
  xlab("Number of single-cell batches")+
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, face = "plain"),
        axis.title.x = element_text(color = "black", size = 16,  face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        # legend.position="none",
        legend.justification="right")
ggsave("figures/FigS2C_NCI19_new.pdf", width = 10, height = 4)

NCI4_added_isoforms_info$Novelty <- "Novel"
NCI4_added_isoforms_info$Novelty[which(NCI4_added_isoforms_info$transcript_catalog == "Known")] <- "Known"
NCI4_added_isoforms_df <- as.data.frame(table(NCI4_added_isoforms_info$Novelty, NCI4_added_isoforms_info$batch_level))
NCI4_added_isoforms_df$pct <- NCI4_added_isoforms_df$Freq/sum(NCI4_added_isoforms_df$Freq)
NCI4_added_isoforms_df$Var1 <- factor(NCI4_added_isoforms_df$Var1, levels = c("Novel","Known"))
n = 2
cols = gg_color_hue(n)
cat.palette = c( "Known"="#6BAED6", 
                 "Novel"="#F8766D")
ggplot(NCI4_added_isoforms_df) +
  geom_bar(aes(x=Var2, y=pct, fill=Var1), position="stack", stat="identity")+
  scale_fill_manual(values = cat.palette)+ ylab("Fraction")+
  xlab("Number of single-cell batches")+
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, face = "plain"),
        axis.title.x = element_text(color = "black", size = 16,  face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        # legend.position="none",
        legend.justification="right")
ggsave("figures/FigS2C_NCI4_new.pdf", width = 10, height = 4)

save(list = ls(), file = "data_for_figure_1016.RData")



# Cell type level
count_mtx <- readRDS("~/Desktop/scLongread/sc_files/count_mtx.rds")
metadata <- readRDS("~/Desktop/scLongread/sc_files/metadata.rds")

count_mtx <- count_mtx[TALON_afterqc_orf_secondpass2$transcript_name_unique,]
table(metadata$Celltype)
metadata$Sample <- as.character(metadata$Celltype)
group <- metadata$Celltype %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx)
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
sum(count_mtx_sum)
count_mtx_sum <- Matrix::t(count_mtx_sum)
dim(count_mtx_sum)
CT_level_expr_mtx <- count_mtx_sum

CT_level_patowary <- rowSums(CT_level_expr_mtx < 3)
keep <- !(CT_level_patowary == 37)

keep2 <- which(rowMeans(count_mtx > 0 ) >= 0.1)
count_mtx_average <- rowMeans(count_mtx)
count_mtx_keep <- count_mtx[keep,]
count_mtx_average["BAMBI-201"]
count_mtx_average["UBA7-TALONT001920798"]
dim(count_mtx_keep)
library(Seurat)
lr <- CreateSeuratObject(counts = count_mtx_keep, meta.data = metadata)
count_mtx_true <- count_mtx > 0
count_mtx_T <- Matrix::t(count_mtx_true)
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
CT_level_true_mtx <- count_mtx_sum
celltype_numbers <- as.vector(table(metadata$Celltype))
celltype_orders <- names(table(metadata$Celltype))
CT_level_true_mtx <- CT_level_true_mtx[celltype_orders,]
CT_level_pct_mtx <- CT_level_true_mtx/celltype_numbers
test_mtx <- (CT_level_pct_mtx > 0.9)
idx1 <- which(colSums(CT_level_pct_mtx > 0.9) > 0)


# check min reads of isoforms tested in DEseq2 and isoQTL
isoform_tested <- do.call(rbind, isoQTL_list)$phenotype_id
isoform_tested <- unique(isoform_tested)
isoform_tested_info <- subset(TALON_afterqc_orf_secondpass, annot_transcript_id %in% isoform_tested)
min(isoform_tested_info$total_counts)
isoform_tested_indv <- Idv_level_expr_mtx[isoform_tested_info$transcript_name_unique,]
isoform_tested_indv_num <- rowSums(isoform_tested_indv > 0)
min(isoform_tested_info$total_counts[(isoform_tested_indv_num > 8)])
Idv_level_expr_mtx["SMN1-TALONT001426912",]
dim(Idv_level_expr_mtx)



linear_model <- lm(df_niso_per_gene$novel_n ~ df_niso_per_gene$known_n)
summary(linear_model)

coloc_final_lungcancer <- subset(coloc_final, GWAS == "Lung_cancer")
lungfun_sus_isoform <- unique((coloc_final_lungfun$transcript_name))
lungcancer_sus_isoform <- unique(c(coloc_final_lungcancer$transcript_name, TWAS_tb$transcript_name))
intersect(lungcancer_sus_isoform,lungfun_sus_isoform)

df_ct_num <- as.data.frame(table(metadata$Celltype))
colnames(df_ct_num) <- c("Celltype", "cell_number")
df_ct_num$Celltype <- gsub(" ", "_", df_ct_num$Celltype)
colnames(DEIs_sc_wilcox)
DEIs_sc_sig <- subset(DEIs_sc_wilcox, avg_log2FC > 1 & p_val_adj < 0.05)
df_DEI_num_sc <- as.data.frame(table(DEIs_sc_sig$cluster))
colnames(df_DEI_num_sc) <- c("Celltype", "DEI_sc")
df_DEI_num_sc$Celltype <- gsub(" ", "_", df_DEI_num_sc$Celltype)
colnames(df_DEI) <- c("Celltype", "DEI_pseudo")
df_DEI_sum <- left_join(df_DEI, df_DEI_num_sc, by = "Celltype")
df_DEI_sum <- left_join(df_DEI_sum, df_ct_num, by = "Celltype")
DEIs_sc_sig$Celltype <- gsub(" ", "_", DEIs_sc_sig$cluster)
celltypes <- levels(metadata$Celltype)
celltypes <- gsub(" ", "_", celltypes)
rownames(df_DEI_sum) <- df_DEI_sum$Celltype
DEI_sig$Celltype <- str_split_fixed(rownames(DEI_sig), "\\.", 2)[,1]
df_DEI_sum$overlaps <- NA
for (celltype in celltypes) {
  tmp_sc <- subset(DEIs_sc_sig, Celltype == celltype)
  tmp_ps <- subset(DEI_sig, Celltype == celltype)
  overlaps <- intersect(tmp_sc$gene, tmp_ps$X)
  df_DEI_sum[celltype,"overlaps"] <- length(overlaps)
}
PPIL6_GWAS <- read.table("gwas_fn.tsv", sep = "\t", header = TRUE)
PPIL6_eQTL <- subset(PPIL6_eQTL,  rsnum %in% isoQTL_fn$rsid)
PPIL6_isoQTL <- subset(PPIL6_isoQTL, rsnum %in% isoQTL_fn$rsid)
PPIL6_sub <- PPIL6_sub[,c(1,2,5,3,4)]
rs12528822_LD <- LDmatrix(PPIL6_isoQTL$rsnum,
                              pop = "EAS", 
                              token = mytoken,
                              genome_build = "grch38"
)
PPIL6_CADD <- read.table(("GRCh38-v1.7_anno_788a2238917383d7a4716982389ed291.tsv"), header = TRUE, 
                         sep = "\t", quote = "")
PPIL6_CADD$variant_id <- paste0(PPIL6_CADD)
write_vcf(PPIL6_sub, "PPIL6_variant.vcf")
colnames(PPIL6_sub) <- c("#CHROM", "POS", "ID", "REF", "ALT")
write.table(PPIL6_sub, "PPIL6_leadisoQTL.vcf", sep = "\t", row.names = FALSE, quote = FALSE)

PPIL6_isoQTL <- left_join(PPIL6_isoQTL, PPIL6_eQTL[,c("rsnum", "pval_nominal")], by = "rsnum")


PPIL6_isoQTL_MPRA <- intersect(PPIL6_isoQTL$rsnum[which(PPIL6_isoQTL$pvalue < 10^-6)], MPRA_V2$RSNUM)
MPRA_V2_isoQTL <- subset(MPRA_V2, RSNUM %in% PPIL6_isoQTL_MPRA)

PPIL6_eQTL_MPRA <- intersect(PPIL6_eQTL_tmp$rsnum[which(PPIL6_eQTL_tmp$pval_nominal < 10^-12)], MPRA_V2$RSNUM)
MPRA_V2_eQTL <- subset(MPRA_V2, RSNUM %in% PPIL6_eQTL_MPRA)

# MUC1
MUC1_coloc_leadSNP <- "rs4072037"
isoQTL_plot_pub("AT2", MUC1_coloc_leadSNP, transcript = "ENST00000368392")
isoQTL_plot_pub("Alveolar_transitional_cells", MUC1_coloc_leadSNP, transcript = "ENST00000368392")
isoQTL_plot_pub("Secretory_transitional_cells", MUC1_coloc_leadSNP, transcript = "ENST00000368392")
isoQTL_plot_pub("Club", MUC1_coloc_leadSNP, transcript = "ENST00000368392")

LDpair("rs12411216", 
       "rs4072037", 
       pop = "EAS", 
       token = mytoken, 
       output = "table", 
       genome_build = "grch38"
)
LDpair("rs12411216", 
       "rs4072037", 
       pop = "EUR", 
       token = mytoken, 
       output = "table", 
       genome_build = "grch38"
)ld_insampe_PPIL6 <- read.table("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/locus_6q21.ld", header = FALSE, sep = "\t")
snp_ld <- read.table("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/locus_6q21.snplist")
colnames(ld_insampe_PPIL6) <- snp_ld$V1
rownames(ld_insampe_PPIL6) <- snp_ld$V1

# 
df_insample_ld <- data.frame(rsnum = rownames(ld_insampe_PPIL6),
                             insample_ld = ld_insampe_PPIL6[,"rs12528822"])
PPIL6_isoQTL <- left_join(PPIL6_isoQTL, df_insample_ld, by = "rsnum")
Y <- read.table(gzfile("/data/Choi_lung/scLongreads/jaxqtl/Multiciliated/phenotype.bed.gz"), header = TRUE, comment.char = "", 
                sep = "\t")
which(Y$trascript_id == "ENST00000521072")
library(susieR)
fitted_rss1 <- susie_rss(bhat = PPIL6_isoQTL$effect, shat = PPIL6_isoQTL$se, n = 127, R = ld_insampe_PPIL6, var_y = var(as.numeric(Y[10015,c(5:131)])), L = 10,
                         estimate_residual_variance = TRUE,coverage = 0.9)

summary(fitted_rss1)$cs
PPIL6_isoQTL$rsnum[fitted_rss1$sets$cs$L1]
susie_plot(fitted_rss1, y="PIP", b=PPIL6_isoQTL$effect)
isoQTL_plot_pub("Multiciliated", "rs12528822", transcript = "ENST00000521072")
isoQTL_plot_pub("Multiciliated", "rs12528822", transcript = "TALONT003040002")
isoQTL_plot_pub("Multiciliated", "rs10080925", transcript = "ENST00000521072")
isoQTL_plot_pub("Multiciliated", "rs10080925", transcript = "TALONT003040002")
isoQTL_plot_pub("Multiciliated", "rs10080925", transcript = "ENST00000424445")
isoQTL_plot_pub("Multiciliated", "rs12528822", transcript = "ENST00000424445")
isoQTL_plot_pub("Multiciliated", "rs12528822", transcript = "TALONT003040049")
rownames(PPIL6_isoQTL) <- PPIL6_isoQTL$rsnum
K562_mtx <- subset(K562_mtx, chr == "chr6")
PPIL6_isoQTL_win_RBP_motif_peaks <- lapply(PPIL6_isoQTL$rsnum, function(snp){
  cCRE <- K562_mtx$peak[which(as.numeric(K562_mtx$start) <= as.numeric(PPIL6_isoQTL[snp,]$pos) & 
                                as.numeric(K562_mtx$end) >= as.numeric(PPIL6_isoQTL[snp,]$pos))]
  RBP <- K562_mtx$RBP[which(as.numeric(K562_mtx$start) <= as.numeric(PPIL6_isoQTL[snp,]$pos) & 
                              as.numeric(K562_mtx$end) >= as.numeric(PPIL6_isoQTL[snp,]$pos))]
  message("Finish ", snp)
  return(paste(cCRE,RBP, sep = "_"))
})
names(PPIL6_isoQTL_win_RBP_motif_peaks) <- PPIL6_isoQTL$rsnum
test.list <- lapply(PPIL6_isoQTL$rsnum, function(isoQTL){
  if (length(PPIL6_isoQTL_win_RBP_motif_peaks[[isoQTL]])==0) {
    return(PPIL6_isoQTL_win_RBP_motif_peaks[[isoQTL]])
  }else{
    x <- PPIL6_isoQTL_win_RBP_motif_peaks[[isoQTL]]
    df <- data.frame(RBP_motif = str_split_fixed(x, "_", 2)[,1],
                     RBP = str_split_fixed(x, "_", 2)[,2],
                     snp = isoQTL)
    return(df)
  }
})
test <- do.call(rbind, test.list)
PPIL6_isoQTL_RBP <- subset(PPIL6_isoQTL, rsnum %in% test$snp)
isoQTL_plot_pub("Multiciliated", "rs3778475", transcript = "ENST00000521072")
isoQTL_plot_pub("Multiciliated", "rs3778475", transcript = "TALONT003040002")


DotPlot(lr.isoform.sub, features = c("CAV1-201","CAV1-202","CAV1-203","CAV1-204"))
DotPlot(lr.isoform.sub, features = c("CAV1-203","CAV1-204", "CAV1-TALONT003062900", "CAV1+CAV2-TALONT003056734"))
DotPlot(lr.isoform.sub, features = c("FGF14-202"))

dim(lr.isoform.sub)
count_mtx <- lr.isoform.sub@assays$Isoform@counts
metadata <- lr.isoform.sub@meta.data # please check this line if it is correct
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass2$transcript_name_unique,]
pct <- rowSums(count_mtx > 0)/ncol(count_mtx)

hist(pct[unique(DEI_sig$X)], breaks = 200)

hist(pct[unique(DEIs_sc_sig$gene)], breaks = 200)
length(intersect(unique(DEI_sig$X), unique(DEIs_sc_sig$gene)))
df_DEI_dist1 <- data.frame(pct = pct[unique(DEI_sig$X)], method = "Indv")
df_DEI_dist2 <- data.frame(pct = pct[unique(DEIs_sc_sig$gene)], method = "sc")
df_DEI_dist <- rbind(df_DEI_dist1, df_DEI_dist2)
sum(df_DEI_dist1$pct <= 0.015, na.rm = TRUE)/nrow(df_DEI_dist1)
sum(df_DEI_dist2$pct <= 0.015, na.rm = TRUE)/nrow(df_DEI_dist2)
ggplot(df_DEI_dist) +
  geom_histogram(aes(x=pct,y=..density..,fill=method),position="identity",alpha=.5,binwidth = 0.005)
rownames(TALON_afterqc_orf_secondpass)

AT1_sc$gene_name <- TALON_afterqc_orf_secondpass[AT1_sc$gene,"annot_gene_name"]
df_niso_perDEI1 <- as.data.frame(table(AT1_sc$gene_name))
df_niso_perDEI2 <- as.data.frame(table(AT1_idv$gene))
hist(df_niso_perDEI1$Freq, breaks = 50)
hist(df_niso_perDEI2$Freq, breaks = 50)
DEIs_sc_sig$gene_name <- TALON_afterqc_orf_secondpass[DEIs_sc_sig$gene,"annot_gene_name"]
df_niso_perDEI1 <- as.data.frame(table(DEIs_sc_sig$gene_name))
df_niso_perDEI2 <- as.data.frame(table(DEI_sig$gene))
df_DEIperGene_dist1 <- data.frame(niso = df_niso_perDEI2$Freq, method = "Indv")
df_DEIperGene_dist2 <- data.frame(niso = df_niso_perDEI1$Freq, method = "sc")
df_DEIperGene_dist <- rbind(df_DEIperGene_dist1, df_DEIperGene_dist2)
ggplot(df_DEIperGene_dist) +
  geom_histogram(aes(x=niso,y=..density..,fill=method),position="identity",alpha=.5,binwidth = 1)
df_DEIperGene_dist1_ct <- DEI_sig %>% group_by(Celltype) %>% summarise(mean = mean(table(gene)))
df_DEIperGene_dist2_ct <- DEIs_sc_sig %>% group_by(cluster) %>% summarise(mean = mean(table(gene_name)))
colnames(df_DEIperGene_dist2_ct)[1] <- "Celltype"
df_DEIperGene_dist2_ct$Celltype <- gsub(" ", "_",df_DEIperGene_dist2_ct$Celltype)
df_DEIperGene_ct <- left_join(df_DEIperGene_dist1_ct, df_DEIperGene_dist2_ct, by = "Celltype", suffix = c("_idv","_sc"))
library(reshape2)
df_DEIperGene_ct_rs <- melt(df_DEIperGene_ct)
library(ggpubr)
ggboxplot(df_DEIperGene_ct_rs, x = "variable", y = "value", fill = "variable") + ylab("Average number of DEIs per gene within each cell type")+
  xlab("Method")

# 1031 on Biowulf
ld_insampe_PPIL6 <- read.table("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/locus_6q21.ld", header = FALSE, sep = "\t")
snp_ld <- read.table("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/locus_6q21.snplist")
colnames(ld_insampe_PPIL6) <- snp_ld$V1
rownames(ld_insampe_PPIL6) <- snp_ld$V1

# 
df_insample_ld <- data.frame(rsnum = rownames(ld_insampe_PPIL6),
                             insample_ld = ld_insampe_PPIL6[,"rs12528822"])
PPIL6_isoQTL <- left_join(PPIL6_isoQTL, df_insample_ld, by = "rsnum")
Y <- read.table(gzfile("/data/Choi_lung/scLongreads/jaxqtl/Multiciliated/phenotype.bed.gz"), header = TRUE, comment.char = "", 
                sep = "\t")
which(Y$trascript_id == "ENST00000521072")
library(susieR)
fitted_rss1 <- susie_rss(bhat = PPIL6_isoQTL$effect, shat = PPIL6_isoQTL$se, n = 127, R = ld_insampe_PPIL6, var_y = var(as.numeric(Y[10015,c(5:131)])), L = 10,
                         estimate_residual_variance = TRUE,coverage = 0.9)

summary(fitted_rss1)$cs
PPIL6_isoQTL$rsnum[fitted_rss1$sets$cs$L1]
susie_plot(fitted_rss1, y="PIP", b=PPIL6_isoQTL$effect)
isoQTL_plot_pub("Multiciliated", "rs12528822", transcript = "ENST00000521072")
isoQTL_plot_pub("Multiciliated", "rs12528822", transcript = "TALONT003040002")
isoQTL_plot_pub("Multiciliated", "rs10080925", transcript = "ENST00000521072")
isoQTL_plot_pub("Multiciliated", "rs10080925", transcript = "TALONT003040002")
isoQTL_plot_pub("Multiciliated", "rs10080925", transcript = "ENST00000424445")
isoQTL_plot_pub("Multiciliated", "rs12528822", transcript = "ENST00000424445")
isoQTL_plot_pub("Multiciliated", "rs12528822", transcript = "TALONT003040049")
rownames(PPIL6_isoQTL) <- PPIL6_isoQTL$rsnum
K562_mtx <- subset(K562_mtx, chr == "chr6")
PPIL6_isoQTL_win_RBP_motif_peaks <- lapply(PPIL6_isoQTL$rsnum, function(snp){
  cCRE <- K562_mtx$peak[which(as.numeric(K562_mtx$start) <= as.numeric(PPIL6_isoQTL[snp,]$pos) & 
                                as.numeric(K562_mtx$end) >= as.numeric(PPIL6_isoQTL[snp,]$pos))]
  RBP <- K562_mtx$RBP[which(as.numeric(K562_mtx$start) <= as.numeric(PPIL6_isoQTL[snp,]$pos) & 
                              as.numeric(K562_mtx$end) >= as.numeric(PPIL6_isoQTL[snp,]$pos))]
  message("Finish ", snp)
  return(paste(cCRE,RBP, sep = "_"))
})
names(PPIL6_isoQTL_win_RBP_motif_peaks) <- PPIL6_isoQTL$rsnum
test.list <- lapply(PPIL6_isoQTL$rsnum, function(isoQTL){
  if (length(PPIL6_isoQTL_win_RBP_motif_peaks[[isoQTL]])==0) {
    return(PPIL6_isoQTL_win_RBP_motif_peaks[[isoQTL]])
  }else{
    x <- PPIL6_isoQTL_win_RBP_motif_peaks[[isoQTL]]
    df <- data.frame(RBP_motif = str_split_fixed(x, "_", 2)[,1],
                     RBP = str_split_fixed(x, "_", 2)[,2],
                     snp = isoQTL)
    return(df)
  }
})
test <- do.call(rbind, test.list)
PPIL6_isoQTL_RBP <- subset(PPIL6_isoQTL, rsnum %in% test$snp)
isoQTL_plot_pub("Multiciliated", "rs3778475", transcript = "ENST00000521072")
isoQTL_plot_pub("Multiciliated", "rs3778475", transcript = "TALONT003040002")


DotPlot(lr.isoform.sub, features = c("CAV1-201","CAV1-202","CAV1-203","CAV1-204"))
DotPlot(lr.isoform.sub, features = c("CAV1-203","CAV1-204", "CAV1-TALONT003062900", "CAV1+CAV2-TALONT003056734"))
DotPlot(lr.isoform.sub, features = c("FGF14-202"))

dim(lr.isoform.sub)
count_mtx <- lr.isoform.sub@assays$Isoform@counts
metadata <- lr.isoform.sub@meta.data # please check this line if it is correct
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass2$transcript_name_unique,]
pct <- rowSums(count_mtx > 0)/ncol(count_mtx)

hist(pct[unique(DEI_sig$X)], breaks = 200)

hist(pct[unique(DEIs_sc_sig$gene)], breaks = 200)
length(intersect(unique(DEI_sig$X), unique(DEIs_sc_sig$gene)))
df_DEI_dist1 <- data.frame(pct = pct[unique(DEI_sig$X)], method = "Indv")
df_DEI_dist2 <- data.frame(pct = pct[unique(DEIs_sc_sig$gene)], method = "sc")
df_DEI_dist <- rbind(df_DEI_dist1, df_DEI_dist2)
sum(df_DEI_dist1$pct <= 0.015, na.rm = TRUE)/nrow(df_DEI_dist1)
sum(df_DEI_dist2$pct <= 0.015, na.rm = TRUE)/nrow(df_DEI_dist2)
ggplot(df_DEI_dist) +
  geom_histogram(aes(x=pct,y=..density..,fill=method),position="identity",alpha=.5,binwidth = 0.005)
rownames(TALON_afterqc_orf_secondpass)

AT1_sc$gene_name <- TALON_afterqc_orf_secondpass[AT1_sc$gene,"annot_gene_name"]
df_niso_perDEI1 <- as.data.frame(table(AT1_sc$gene_name))
df_niso_perDEI2 <- as.data.frame(table(AT1_idv$gene))
hist(df_niso_perDEI1$Freq, breaks = 50)
hist(df_niso_perDEI2$Freq, breaks = 50)
DEIs_sc_sig$gene_name <- TALON_afterqc_orf_secondpass[DEIs_sc_sig$gene,"annot_gene_name"]
df_niso_perDEI1 <- as.data.frame(table(DEIs_sc_sig$gene_name))
df_niso_perDEI2 <- as.data.frame(table(DEI_sig$gene))
df_DEIperGene_dist1 <- data.frame(niso = df_niso_perDEI2$Freq, method = "Indv")
df_DEIperGene_dist2 <- data.frame(niso = df_niso_perDEI1$Freq, method = "sc")
df_DEIperGene_dist <- rbind(df_DEIperGene_dist1, df_DEIperGene_dist2)
ggplot(df_DEIperGene_dist) +
  geom_histogram(aes(x=niso,y=..density..,fill=method),position="identity",alpha=.5,binwidth = 1)
df_DEIperGene_dist1_ct <- DEI_sig %>% group_by(Celltype) %>% summarise(mean = mean(table(gene)))
df_DEIperGene_dist2_ct <- DEIs_sc_sig %>% group_by(cluster) %>% summarise(mean = mean(table(gene_name)))
colnames(df_DEIperGene_dist2_ct)[1] <- "Celltype"
df_DEIperGene_dist2_ct$Celltype <- gsub(" ", "_",df_DEIperGene_dist2_ct$Celltype)
df_DEIperGene_ct <- left_join(df_DEIperGene_dist1_ct, df_DEIperGene_dist2_ct, by = "Celltype", suffix = c("_idv","_sc"))
library(reshape2)
df_DEIperGene_ct_rs <- melt(df_DEIperGene_ct)
library(ggpubr)
ggboxplot(df_DEIperGene_ct_rs, x = "variable", y = "value", fill = "variable") + ylab("Average number of DEIs per gene within each cell type")+
  xlab("Method")


feature_chunks <- split(PPIL6_isoQTL$rsnum, 
                        cut(seq_along(PPIL6_isoQTL$rsnum), 20, labels=FALSE))
for (i in 1:20) {
  output <- feature_chunks[[i]]
  output_tb <- subset(PPIL6_sub, ID %in% output)
  write.table(PPIL6_sub, file = paste0("PPIL6_leadisoQTL_chunk",i ,".vcf"), sep = "\t", row.names = FALSE, quote = FALSE)
}

PPIL6_CADD <- read.table(("GRCh38-v1.7_anno_6c9e1e4292660f582bfe7367dadd6931.tsv.gz"), header = TRUE, 
                         sep = "\t",comment.char = "", skip = 1)

write.table(colnames(PPIL6_CADD), file = "CADD_categories.txt", quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = FALSE)


PPIL6_sub$variant_id <- paste(PPIL6_sub$POS, PPIL6_sub$REF,PPIL6_sub$ALT, sep = "_")
PPIL6_CADD$variant_id <- paste(PPIL6_CADD$Pos, PPIL6_CADD$Ref,PPIL6_CADD$Alt, sep = "_")
length(unique(PPIL6_CADD$variant_id))
PPIL6_CADD <- left_join(PPIL6_CADD, PPIL6_sub[,c("variant_id","ID")], by = "variant_id")

# DEI
DEI_sc_info <- subset(TALON_afterqc_orf_secondpass2, transcript_name_unique %in% DEIs_sc_wilcox$gene)
DEI_idv_info <- subset(TALON_afterqc_orf_secondpass2, transcript_name_unique %in% DEI_sig$X)
sum(DEI_idv_info$unique_dataset == 1)
sum(DEI_sc_info$unique_dataset == 1)

isoQTL_plot_pub("Club", MUC1_coloc_leadSNP, transcript = "ENST00000368392")

isoQTL_plot_pub_v2("Multiciliated", "rs12528822", transcript = "ENST00000521072")
ggsave("", p,width = 15, height = 7)
isoQTL_plot_pub_v2("Multiciliated", "rs12528822", transcript = "TALONT003040002")


unique(TWAS_tb_imm_epi$transcript_name[which(TWAS_tb_imm_epi$Celltype %in% celltypes[1:9])])
unique(TWAS_tb_imm_epi$transcript_name[which(TWAS_tb_imm_epi$Celltype %in% celltypes[10:26])])

tb.list <- lapply(c(1:1134), function(i){
  tryCatch({
    tb <- LDpair(AF_diff_new_sGenes_sQTL_tmp$rsid[i], 
                 AF_diff_new_sGenes_sQTL_tmp$rs_id_dbSNP151_GRCh38p7[i], 
                 pop = "EAS", 
                 token = mytoken, 
                 output = "table", 
                 genome_build = "grch38"
    )
  }, error=function(e){})
  
})
library(data.table)
tb_EAS <- rbindlist(tb.list, fill = TRUE)
tb_EUR <- rbindlist(tb.list.EUR, fill = TRUE)
sum(grepl("linkage equilibrium", tb_EAS$ld)) # 437
isoQTL_nonsQTL <- unique(AF_diff_new_sGenes_sQTL_tmp[!AF_diff_new_sGenes_sQTL_tmp$replicated_V8,"rsid"])
tb_EAS_nosQTL <- tb_EAS[(tb_EAS$var1 %in% isoQTL_nonsQTL),]
sum(grepl("linkage equilibrium", tb_EAS_nosQTL$ld)) # 415

tb_EUR_nosQTL <- tb_EUR[(tb_EUR$var1 %in% isoQTL_nonsQTL),]
sum(grepl("linkage equilibrium", tb_EUR_nosQTL$ld)) # 415
