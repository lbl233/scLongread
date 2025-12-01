setwd("/Users/lib14/Desktop/scLongread/")


library(locuscomparer)
isoQTL_fn = 'AT2.15q21.2.tsv'
eQTL_fn = 'Shrine.ratio.15q21.2.tsv'
gwas_fn = 'Shi.EAS.15q21.2.tsv'
locuscompare(in_fn1 = eQTL_fn, in_fn2 = isoQTL_fn, title = 'eQTL', title2 = 'isoQTL',
             population = "EAS", snp = "rs7169304")
library(ggplot2)

p1 <- locuscompare(in_fn1 = gwas_fn, in_fn2 = isoQTL_fn, title = 'GWAS', title2 = 'isoQTL',
             population = "EAS", snp = "rs7169304")
p2 <- locuscompare(in_fn1 = gwas_fn, in_fn2 = eQTL_fn, title = 'GWAS', title2 = 'eQTL',
             population = "EAS", snp = "rs7169304")
p3 <- locuscompare(in_fn1 = isoQTL_fn, in_fn2 = gwas_fn, title = 'isoQTL', title2 = 'GWAS',
                   population = "EAS", snp = "rs7169304")
ggsave("./figures/Fig7B1.pdf", p1, width = 12, height = 8)
ggsave("./figures/Fig7B2.pdf", p2, width = 12, height = 8)
ggsave("./figures/Fig7B3.pdf", p3, width = 12, height = 8)

mytoken = "149dc339070c"
PPIL6_isoQTL <- read.table("PPIL6_isoQTL.txt", header = TRUE, sep = "\t")
PPIL6_eQTL <- read.table("PPIL6_eQTL.txt", header = TRUE, sep = "\t")
PPIL6_eQTL_tmp <- read.table("PPIL6_eQTL.txt", header = TRUE, sep = "\t")
PPIL6_eQTL <- read.table("", header = TRUE, sep = "\t")
test_tmp <- LDmatrix(head(PPIL6_isoQTL$rsnum[order(PPIL6_isoQTL$pvalue)] , n =100), pop = "EAS",token = mytoken, genome_build = "grch38")

df_test_tmp <- test_tmp[,c("RS_number","rs12528822")]
colnames(df_test_tmp)[1] <- "rsnum"
rs_list <- intersect(colnames(genotypes), PPIL6_isoQTL$rsnum)
PPIL6_isoQTL <- left_join(PPIL6_isoQTL, df_test_tmp, by = "rsnum")

geno_tmp <- genotypes[,PPIL6_isoQTL$rsnum]

rs12528822_LDproxy <- LDproxy("rs12528822", 
       pop = "EAS", 
       token = mytoken,
       genome_build = "grch38"
)
# isoQTL lead LD
tb.list <- lapply(PPIL6_isoQTL$rsnum, function(x){
  tryCatch({
    tb <- LDpair("rs12528822", x, pop = "EAS", token = mytoken, genome_build = "grch38")
  }, error=function(e){})
})

df_tmp <- as.data.frame(matrix(as.character(NA), nrow = 1, ncol = 18))

colnames(df_tmp) <- colnames(tb.list[[1]])
tb.list.tmp <- lapply(tb.list, function(x){
  if(is.null(x)){return(df_tmp)}else{return(x[,c(1:18)])}
})

# eQTL lead
library(LDlinkR)
tb.list <- lapply(PPIL6_eQTL$rsnum, function(x){
  tryCatch({
    tb <- LDpair("rs736830", x, pop = "EAS", token = mytoken, genome_build = "grch38")
  }, error=function(e){})
})

df_tmp <- as.data.frame(matrix(as.character(NA), nrow = 1, ncol = 18))

colnames(df_tmp) <- colnames(tb.list[[1]])
tb.list.tmp <- lapply(tb.list, function(x){
  if(is.null(x)){return(df_tmp)}else{return(x[,c(1:18)])}
})
library(data.table)
PPIL6_isoQTL_LD <- rbindlist(tb.list, fill = TRUE)
PPIL6_isoQTL_LD_R2 <- PPIL6_isoQTL_LD[,c("var2", "r2")]
colnames(PPIL6_isoQTL_LD_R2)[1] <- "rsnum"
PPIL6_isoQTL_LD_Dprime <- PPIL6_isoQTL_LD[,c("var2", "d_prime")]
colnames(PPIL6_isoQTL_LD_Dprime)[1] <- "rsnum"
PPIL6_isoQTL <- left_join(PPIL6_isoQTL, PPIL6_isoQTL_LD_R2, by = "rsnum")
PPIL6_isoQTL <- left_join(PPIL6_isoQTL, PPIL6_isoQTL_LD_Dprime, by = "rsnum")

PPIL6_isoQTL_LD <- rbindlist(tb.list, fill = TRUE)
PPIL6_eQTL_LD_R2 <- PPIL6_isoQTL_LD[,c("var2", "r2")]
colnames(PPIL6_eQTL_LD_R2)[1] <- "rsnum"
PPIL6_eQTL <- left_join(PPIL6_eQTL, PPIL6_eQTL_LD_R2, by = "rsnum")


PPIL6_isoQTL_lead <- subset(PPIL6_isoQTL, r2 > 0.8 & pvalue < 0.0001)
PPIL6_eQTL_lead <- subset(PPIL6_eQTL, r2 > 0.8 & pval_nominal < 0.0001)
PPIL6_eQTL_lead_CADD <- subset(PPIL6_CADD, ID %in% PPIL6_eQTL_lead$rsnum)
PPIL6_isoQTL_lead_CADD <- subset(PPIL6_CADD, ID %in% PPIL6_isoQTL_lead$rsnum)
PHRED_eQTL <- PPIL6_eQTL_lead_CADD$PHRED[!duplicated(PPIL6_eQTL_lead_CADD$ID)]
PHRED_isoQTL <- PPIL6_isoQTL_lead_CADD$PHRED[!duplicated(PPIL6_isoQTL_lead_CADD$ID)]
mean(PHRED_isoQTL)
mean(PHRED_eQTL)
wilcox.test(PHRED_isoQTL, PHRED_eQTL )
library(readxl)
CADD_category <- read_xlsx("CADD_annotation.xlsx")
idx <- grep("splicing",CADD_category$`Functional type`)
PPIL6_eQTL_lead_CADD_splicing_mtx <- PPIL6_eQTL_lead_CADD[!duplicated(PPIL6_eQTL_lead_CADD$ID),idx[-c(1:2)]]
PPIL6_isoQTL_lead_CADD_splicing_mtx <- PPIL6_isoQTL_lead_CADD[!duplicated(PPIL6_isoQTL_lead_CADD$ID),idx[-c(1:2)]]
PPIL6_eQTL_lead_CADD_splicing_mtx_nonzero <- (PPIL6_eQTL_lead_CADD_splicing_mtx > 0 | PPIL6_eQTL_lead_CADD_splicing_mtx < 0)
PPIL6_isoQTL_lead_CADD_splicing_mtx_nonzero <- (PPIL6_isoQTL_lead_CADD_splicing_mtx > 0 | PPIL6_isoQTL_lead_CADD_splicing_mtx < 0)
rowSums(PPIL6_eQTL_lead_CADD_splicing_mtx_nonzero, na.rm = TRUE)
sum(rowSums(PPIL6_eQTL_lead_CADD_splicing_mtx_nonzero, na.rm = TRUE) > 0 )
rowSums(PPIL6_isoQTL_lead_CADD_splicing_mtx_nonzero, na.rm = TRUE)
sum(rowSums(PPIL6_isoQTL_lead_CADD_splicing_mtx_nonzero, na.rm = TRUE) > 0 )

PPIL6_eQTL_lead_CADD$PHRED[!duplicated(PPIL6_eQTL_lead_CADD$ID)]

# PPIL6 Alpha fold
# 202 
# https://alphafold.ebi.ac.uk/search/sequence/MARPQPCGPPHARCGSPSLPERPLQVKVVGLFSCPNFQIAKSAAEELKNETWEYSSSVISFVNGQFLGDALDLQKWAHEVWDIVDIKPSALYDALTEDFSAKFLRDTKHDFVFLDICIDSSPIGRLIFELYCDVCPKTCKNFQVLCTGKAGFSQRGIRLHYKNSIFHRIVQNGWIQGGDIVYGKGDNGESIYGPTFEDENFSVPHNKRGVLGMANKGRHSNGSQFYITLQATPYLDRKFVAFGQLIEGTEVLKQLELVPTQNERPIHMCRITDSGDPYA

# T0002
# https://alphafold.ebi.ac.uk/search/sequence/MARPQPCGPPHARCGSPSLPERPLQVKVVGLFSCPNFQIAKSAAENLKNNHPSKFEDPILVPLQEFAWHQYLQEKKRELKNETWEYSSSVISFVNGQFLGDALDLQKWAHEVWDIVDIKPSALYDALTEDFSAKFLRDTKHDFVFLDICIDSSPIGRLIFELYCDVCPKTCKNFQVLCTGKAGFSQRGIRLHYKNSIFHRIVQNGWIQGGDIVYGKGDNGESIYGPTFEGN

# 207
# https://alphafold.ebi.ac.uk/search/sequence/MARPQPCGPPHARCGSPSLPERPLQVKVVGLFSCPNFQIAKSAAENLKNNHPSKFEDPILVPLQEFAWHQYLQEKKRELKNETWEYSSSVISFVNGQFLGDALDLQKWAHEVWDIVDIKPSALYDALTEDFSAKFLRDTKHDFVFLDICIDSSPIGRLIFELYCDVCPKTCKNFQVLCTGKAGFSQRGIRLHYKNSIFHRIVQNGWIQGGDIVYGKGDNGESIYGPTFEDENFSVPHNKRGVLGMANKGRHSNGSQFYITLQATPYLDRKFVAFGQLIEGTEVLKQLELVPTQNERPIHMCRITDSGDPYA
ENST00000424445
ENST00000521072
TALONT003040002
isoQTL_plot_pub("Multiciliated", "rs12528822", transcript = "ENST00000424445")
ggsave("./figures/Fig11_1_D.pdf",width = 6, height = 6)
isoQTL_plot_pub("Multiciliated", "rs12528822", transcript = "ENST00000521072")
ggsave("./figures/Fig11_1_B.pdf",width = 6, height = 6)
isoQTL_plot_pub("Multiciliated", "rs12528822", transcript = "TALONT003040002")
ggsave("./figures/Fig11_1_F.pdf",width = 6, height = 6)

PPIL6_lead_cs_eQTL <- read.table("PPIL6_cs_for_BL.tsv", sep = "\t", header = TRUE)
PPIL6_lead_cs_eQTL <- str_split_fixed(PPIL6_lead_cs_eQTL$rsid, pattern = "_", 2)[,1]

PPIL6_isoQTL$rsnum[fitted_rss1$sets$cs$L1]
PPIL6_eQTL_lead_CADD <- subset(PPIL6_CADD, ID %in% PPIL6_lead_cs_eQTL)
PPIL6_isoQTL_lead_CADD <- subset(PPIL6_CADD, ID %in% PPIL6_isoQTL$rsnum[fitted_rss1$sets$cs$L1])
PPIL6_isoQTL_lead_CADD <- subset(PPIL6_CADD, ID %in% PPIL6_isoQTL_cs)
write.table(PPIL6_isoQTL_lead_CADD, file = "PPIL6_isoQTL_cs_CADD.txt", sep = "\t", row.names = FALSE, quote = FALSE)
PHRED_eQTL <- PPIL6_eQTL_lead_CADD$PHRED[!duplicated(PPIL6_eQTL_lead_CADD$ID)]
PHRED_isoQTL <- PPIL6_isoQTL_lead_CADD$PHRED[!duplicated(PPIL6_isoQTL_lead_CADD$ID)]
mean(PHRED_isoQTL)
mean(PHRED_eQTL)
wilcox.test(PHRED_isoQTL, PHRED_eQTL )
library(readxl)
CADD_category <- read_xlsx("CADD_annotation.xlsx")
idx <- grep("splicing",CADD_category$`Functional type`)
PPIL6_eQTL_lead_CADD_splicing_mtx <- PPIL6_eQTL_lead_CADD[!duplicated(PPIL6_eQTL_lead_CADD$ID),idx[-c(1:2)]]
PPIL6_isoQTL_lead_CADD_splicing_mtx <- PPIL6_isoQTL_lead_CADD[!duplicated(PPIL6_isoQTL_lead_CADD$ID),idx[-c(1:2)]]
PPIL6_eQTL_lead_CADD_splicing_mtx_nonzero <- (PPIL6_eQTL_lead_CADD_splicing_mtx > 0 | PPIL6_eQTL_lead_CADD_splicing_mtx < 0)
PPIL6_isoQTL_lead_CADD_splicing_mtx_nonzero <- (PPIL6_isoQTL_lead_CADD_splicing_mtx > 0 | PPIL6_isoQTL_lead_CADD_splicing_mtx < 0)
rowSums(PPIL6_eQTL_lead_CADD_splicing_mtx_nonzero, na.rm = TRUE)
sum(rowSums(PPIL6_eQTL_lead_CADD_splicing_mtx_nonzero, na.rm = TRUE) > 0 )
rowSums(PPIL6_isoQTL_lead_CADD_splicing_mtx_nonzero, na.rm = TRUE)
sum(rowSums(PPIL6_isoQTL_lead_CADD_splicing_mtx_nonzero, na.rm = TRUE) > 0 )

PPIL6_eQTL_lead_CADD$PHRED[!duplicated(PPIL6_eQTL_lead_CADD$ID)]


Multiciliated_idx <- which(metadata$Celltype == "Multiciliated")
sample_annot <- as.character(metadata$Sample[Multiciliated_idx])

group <- sample_annot %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx[,Multiciliated_idx])
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
count_mtx_sum <- Matrix::t(count_mtx_sum)
Multiciliated_idv_idv <- colSums(count_mtx_sum)
Multiciliated_PPIL6_idv <- count_mtx_sum[grepl("PPIL6", rownames(count_mtx_sum)),]
Multiciliated_PPIL6_idv_sum <- colSums(Multiciliated_PPIL6_idv)
Multiciliated_PPIL6_idv_count <- count_list$Multiciliated
Multiciliated_PPIL6_idv_total_norm <- log1p((Multiciliated_PPIL6_idv_sum/Multiciliated_idv_idv)*1000)
# PPIL6 eIsoforms
PPIL6_eIsoforms <- final_table$transcript_name[which(final_table$gene_name == "PPIL6")]
PPIL6_eIsoforms_id <- final_table$phenotype_id[which(final_table$gene_name == "PPIL6")]
matrix_pct <- Multiciliated_PPIL6_idv_count[PPIL6_eIsoforms_id,-c(1:4)]
matrix_pct_final <- t(t(matrix_pct)/Multiciliated_PPIL6_idv_sum[colnames(matrix_pct)])
isoQTL_plot_pct(celltype = "Multiciliated", matrix_pct = matrix_pct_final, 
                rs = "rs12528822", transcript = "ENST00000424445")
isoQTL_plot_pct(celltype = "Multiciliated", matrix_pct = matrix_pct_final, 
                rs = "rs12528822", transcript = "ENST00000521072")
isoQTL_plot_pct(celltype = "Multiciliated", matrix_pct = matrix_pct_final, 
                rs = "rs12528822", transcript = "TALONT003040002")
isoQTL_plot_pct(celltype = "Multiciliated", matrix_pct = matrix_pct_final, 
                rs = "rs12528822", transcript = "TALONT003040049")
samples <- colnames(matrix_pct_final)
rs = "rs12528822"
transcript = "TALONT003040002"
datainput2 <- data.frame(snp = genotypes[samples,rs],
                         isoform = as.numeric(matrix_pct_final[transcript,]))
datainput2$snp <- factor(datainput2$snp, levels = c("03","02","01"))
datainput2$snp <- as.numeric(datainput2$snp)
fit_TALONT003040002 <- lm(isoform~snp,data = datainput2)
summary(fit_TALONT003040002)

transcript = "ENST00000424445"
datainput2 <- data.frame(snp = genotypes[samples,rs],
                         isoform = as.numeric(matrix_pct_final[transcript,]))
datainput2$snp <- factor(datainput2$snp, levels = c("03","02","01"))
datainput2$snp <- as.numeric(datainput2$snp)
fit_ENST00000424445 <- lm(isoform~snp,data = datainput2)
summary(fit_ENST00000424445)

transcript = "ENST00000521072"
datainput2 <- data.frame(snp = genotypes[samples,rs],
                         isoform = as.numeric(matrix_pct_final[transcript,]))
datainput2$snp <- factor(datainput2$snp, levels = c("03","02","01"))
datainput2$snp <- as.numeric(datainput2$snp)
fit_ENST00000521072 <- lm(isoform~snp,data = datainput2)
summary(fit_ENST00000521072)

transcript = "TALONT003040049"
datainput2 <- data.frame(snp = genotypes[samples,rs],
                         isoform = as.numeric(matrix_pct_final[transcript,]))
datainput2$snp <- factor(datainput2$snp, levels = c("03","02","01"))
datainput2$snp <- as.numeric(datainput2$snp)
fit_TALONT003040049 <- lm(isoform~snp,data = datainput2)
summary(fit_TALONT003040049)

Multiciliated_PPIL6_idv_total_norm_mtx <- t(as.matrix(Multiciliated_PPIL6_idv_total_norm))
rownames(Multiciliated_PPIL6_idv_total_norm_mtx) <- "PPIL6"
samples <- colnames(Multiciliated_PPIL6_idv_total_norm_mtx)

datainput2 <- data.frame(snp = genotypes[samples,rs],
                         isoform = as.numeric(Multiciliated_PPIL6_idv_total_norm_mtx["PPIL6",]))
datainput2$snp_ra <- "0|0"
datainput2$snp_ra[which(datainput2$snp == "01")] <- "1|1"
datainput2$snp_ra[which(datainput2$snp == "02")] <- "0|1"
df <- as.vector(table(datainput2$snp_ra))
datainput2$snp_ra <- paste0("0|0\n", "n=", df[1])
datainput2$snp_ra[which(datainput2$snp == "01")] <- paste0("1|1\n", "n=", df[3])
datainput2$snp_ra[which(datainput2$snp == "02")] <- paste0("0|1\n", "n=", df[2])
p2 <- ggplot(datainput2, aes(x=snp_ra, y=isoform, color=snp_ra)) +
  geom_violin(trim=FALSE,linewidth = 0.75)+
  geom_boxplot(aes(middle=mean(isoform)),
               width=0.1, linewidth = 0.75)+
  geom_smooth(mapping = aes(x = snp_ra, y = isoform, group = 1),formula = y~x, color = "gray",
              method='lm', size = 1, se =TRUE,fill = alpha("gray", .5) ) +
  labs(title="",x=rs, y = paste0("Normalized expression of ", transcript))+
  scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1")) +
  ggtitle(paste0("Long-read data of ", celltype, "\n",transcript , ": ",rs, ": ", paste(snp_info[rs,1], snp_info[rs,3], snp_info[rs,5], snp_info[rs,4], "b38", sep = "_")))+
  theme(axis.text.x = element_text(color = "black", size = 12,face = "plain"),
        axis.text.y = element_text(color = "black", size = 12,face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.position="none")
p2

# short read data

samples <- colnames(gene_list$Multiciliated)[-c(1:4)]
which(gene_list$Multiciliated$gene_id == "ENSG00000185250")
datainput2 <- data.frame(snp = genotypes[samples,rs],
                         isoform = as.numeric(gene_list$Multiciliated[593,-c(1:4)]))
datainput2$snp_ra <- "0|0"
datainput2$snp_ra[which(datainput2$snp == "01")] <- "1|1"
datainput2$snp_ra[which(datainput2$snp == "02")] <- "0|1"
df <- as.vector(table(datainput2$snp_ra))
datainput2$snp_ra <- paste0("0|0\n", "n=", df[1])
datainput2$snp_ra[which(datainput2$snp == "01")] <- paste0("1|1\n", "n=", df[3])
datainput2$snp_ra[which(datainput2$snp == "02")] <- paste0("0|1\n", "n=", df[2])
p2 <- ggplot(datainput2, aes(x=snp_ra, y=isoform, color=snp_ra)) +
  geom_violin(trim=FALSE,linewidth = 0.75)+
  geom_boxplot(aes(middle=mean(isoform)),
               width=0.1, linewidth = 0.75)+
  geom_smooth(mapping = aes(x = snp_ra, y = isoform, group = 1),formula = y~x, color = "gray",
              method='lm', size = 1, se =TRUE,fill = alpha("gray", .5) ) +
  labs(title="",x=rs, y = paste0("Standardized and normalized expression of ", transcript))+
  scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1")) +
  ggtitle(paste0("Short-read data of ", celltype, "\n",transcript , ": ",rs, ": ", paste(snp_info[rs,1], snp_info[rs,3], snp_info[rs,5], snp_info[rs,4], "b38", sep = "_")))+
  theme(axis.text.x = element_text(color = "black", size = 12,face = "plain"),
        axis.text.y = element_text(color = "black", size = 12,face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.position="none")
p2
rst.list <- readRDS(file = "/data/Choi_lung/scLongreads/colocalization/Byun_total_coloc_rst_list.rds")
Multiciliated <- rst.list$Multiciliated
saveRDS(Multiciliated, file = "/your/path/Byun_total_coloc_rst_list_Multiciliated.rds")

ld_mtx <- read.table("locus_6q21_snp_flcca_ctrl_r.ld", header = FALSE, sep = "\t")
snp_ld <- read.table("locus_6q21_snp_flcca_ctrl_r.snplist")
colnames(ld_mtx) <- snp_ld$V1
rownames(ld_mtx) <- snp_ld$V1
ld_mtx <- as.matrix(ld_mtx)
dim(ld_mtx)
PPIL6_isoQTL_susie <- subset(PPIL6_isoQTL, rsnum %in% snp_ld$V1)
which(PPIL6_isoQTL_susie$rsnum != snp_ld$V1)
Y = as.numeric(count_list$Multiciliated[which(count_list$Multiciliated$trascript_id == "ENST00000521072"),-c(1:4)])
fitted_rss1 <- susie_rss(bhat = PPIL6_isoQTL_susie$effect, 
                         shat = PPIL6_isoQTL_susie$se, 
                         n = 4544, 
                         R = ld_mtx, 
                         var_y = var(Y), 
                         L = 10,
                         estimate_residual_variance = TRUE)
summary(fitted_rss1)$cs

fitted_rss2 = susie_rss(z = PPIL6_isoQTL_susie$zscore, R = ld_mtx, n = 4544, L = 10,
                        estimate_residual_variance = TRUE)
summary(fitted_rss2)$cs
fitted_rss3 <- susie_rss(bhat = PPIL6_isoQTL_susie$effect, 
                         shat = PPIL6_isoQTL_susie$se, 
                         n = 4544, 
                         R = ld_mtx, 
                         var_y = var(Y), 
                         L = 10,
                         coverage = 0.8,
                         estimate_residual_variance = TRUE)
summary(fitted_rss3)$cs

PPIL6_isoQTL_cs <- PPIL6_isoQTL_susie$rsnum[fitted_rss1$sets$cs$L1]
PPIL6_isoQTL_lead_CADD <- subset(PPIL6_CADD, ID %in% PPIL6_isoQTL_cs)

Byun_6q21 <- Byun_total_meta_gwas_31loci_1MB$`6q21`
GWAS_6q21_byun <- Byun_6q21 %>% group_by(variant_id) %>% summarise(pvalue = min(p_value))
colnames(GWAS_6q21_byun) <- c("rsid","pval")
write.table(GWAS_6q21_byun, file = "GWAS_6q21_byun.tsv", row.names = FALSE, quote = FALSE, sep = "\t")






# MUC1 Percentage
AT2_idx <- which(metadata$Celltype == "AT2")
sample_annot <- as.character(metadata$Sample[AT2_idx])

group <- sample_annot %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx[,AT2_idx])
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
count_mtx_sum <- Matrix::t(count_mtx_sum)
AT2_idv_idv <- colSums(count_mtx_sum)
AT2_MUC1_idv <- count_mtx_sum[grepl("MUC1", rownames(count_mtx_sum)),]
AT2_MUC1_idv_sum <- colSums(AT2_MUC1_idv)
AT2_MUC1_idv_count <- count_list$AT2

# MUC1 eIsoforms
MUC1_eIsoforms <- unique(final_table$transcript_name[which(final_table$gene_name == "MUC1")])
MUC1_eIsoforms_id <- unique(final_table$phenotype_id[which(final_table$gene_name == "MUC1")])
matrix_pct <- AT2_MUC1_idv_count[MUC1_eIsoforms_id,-c(1:4)]
matrix_pct_final <- t(t(matrix_pct)/AT2_MUC1_idv_sum[colnames(matrix_pct)])
isoQTL_plot_pct(celltype = "AT2", matrix_pct = matrix_pct_final, 
                rs = "rs4072037", transcript = "ENST00000368392")



# MUC1 Percentage
Secretory_transitional_cells_idx <- which(metadata$Celltype == "Secretory transitional cells")
sample_annot <- as.character(metadata$Sample[Secretory_transitional_cells_idx])

group <- sample_annot %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx[,Secretory_transitional_cells_idx])
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
count_mtx_sum <- Matrix::t(count_mtx_sum)
Secretory_transitional_cells_idv_idv <- colSums(count_mtx_sum)
Secretory_transitional_cells_MUC1_idv <- count_mtx_sum[grepl("MUC1", rownames(count_mtx_sum)),]
Secretory_transitional_cells_MUC1_idv_sum <- colSums(Secretory_transitional_cells_MUC1_idv)
Secretory_transitional_cells_MUC1_idv_count <- count_list$Secretory_transitional_cells

# MUC1 eIsoforms
MUC1_eIsoforms <- unique(final_table$transcript_name[which(final_table$gene_name == "MUC1")])
MUC1_eIsoforms_id <- unique(final_table$phenotype_id[which(final_table$gene_name == "MUC1")])
matrix_pct <- Secretory_transitional_cells_MUC1_idv_count[MUC1_eIsoforms_id,-c(1:4)]
matrix_pct_final <- t(t(matrix_pct)/Secretory_transitional_cells_MUC1_idv_sum[colnames(matrix_pct)])
isoQTL_plot_pct(celltype = "Secretory_transitional_cells", matrix_pct = matrix_pct_final, 
                rs = "rs4072037", transcript = "ENST00000368392")

# MUC1 Percentage
Alveolar_transitional_cells_idx <- which(metadata$Celltype == "Alveolar transitional cells")
sample_annot <- as.character(metadata$Sample[Alveolar_transitional_cells_idx])

group <- sample_annot %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx[,Alveolar_transitional_cells_idx])
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
count_mtx_sum <- Matrix::t(count_mtx_sum)
Alveolar_transitional_cells_idv_idv <- colSums(count_mtx_sum)
Alveolar_transitional_cells_MUC1_idv <- count_mtx_sum[grepl("MUC1", rownames(count_mtx_sum)),]
Alveolar_transitional_cells_MUC1_idv_sum <- colSums(Alveolar_transitional_cells_MUC1_idv)
Alveolar_transitional_cells_MUC1_idv_count <- count_list$Alveolar_transitional_cells

# MUC1 eIsoforms
MUC1_eIsoforms <- unique(final_table$transcript_name[which(final_table$gene_name == "MUC1")])
MUC1_eIsoforms_id <- unique(final_table$phenotype_id[which(final_table$gene_name == "MUC1")])
matrix_pct <- Alveolar_transitional_cells_MUC1_idv_count[MUC1_eIsoforms_id,-c(1:4)]
matrix_pct_final <- t(t(matrix_pct)/Alveolar_transitional_cells_MUC1_idv_sum[colnames(matrix_pct)])
isoQTL_plot_pct(celltype = "Alveolar_transitional_cells", matrix_pct = matrix_pct_final, 
                rs = "rs4072037", transcript = "ENST00000368392")




