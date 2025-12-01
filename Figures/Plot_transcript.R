# Plotting transcripts

library(magrittr)
library(dplyr)
library(ggtranscript)
library(ggplot2)
library(rtracklayer)
library(geomtextpath)
source("/data/lib14/R/Rscripts/utilities.R")
table(TALON_afterqc_orf_secondpass$structural_category)
idx1 <- which(TALON_afterqc_orf_secondpass2$structural_category == "fusion")
TALON_afterqc_orf_secondpass2$transcript_catalog <- TALON_afterqc_orf_secondpass2$transcript_novelty
TALON_afterqc_orf_secondpass2$transcript_catalog[idx1] <- "fusion"
table(TALON_afterqc_orf_secondpass2$transcript_catalog)
cat.palette = c( "Known"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "fusion"="goldenrod1", "Genomic"="#969696", "Antisense"="#66C2A4")
NB.nominal.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_sig_list.rds")
AT2_nominal <- NB.nominal.sig.list$AT2
Multiciliated_nominal <- NB.nominal.sig.list$Multiciliated
gtf <- rtracklayer::import("/data/Choi_lung/scLongreads/TALON_workspace/test10s/NCI_lung_secondpass.isoforms_w_known.gtf")
# gtf <- rtracklayer::import("/data/Choi_lung/lbl/refdata-gex-GRCh38-2020-A/genes/genes.gtf")

class(gtf)


gtf <- gtf %>% dplyr::as_tibble()

class(gtf)


# filter your gtf for the gene of interest, here "OAS1"
gene_of_interest <- "ENSG00000111671" #SPSB2
gene_of_interest <- "ENSG00000164265"
gene_of_interest <- "ENSG00000089127"
gene_of_interest <- "ENSG00000186501" # TMEM222
gene_of_interest <- "ENSG00000196735" # HLA-DQA1
gene_of_interest <- "ENSG00000185250" # PPIL6
SPSB2_annotation_from_gtf <- gtf %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_id == gene_of_interest
  ) %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_name
    #transcript_biotype
  )

# extract the required annotation columns
SPSB2_annotation_from_gtf <- SPSB2_annotation_from_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_name
    #transcript_biotype
  )

SPSB2_annotation_from_gtf %>% head()


# extract exons
SPSB2_exons <- SPSB2_annotation_from_gtf %>% dplyr::filter(type == "exon")


SPSB2_exons$transcript_name_unique <- SPSB2_exons$transcript_name
SPSB2_exons$transcript_name_unique[grepl("TALON", SPSB2_exons$transcript_name)] <- paste("PPIL6", SPSB2_exons$transcript_name[grepl("TALON", SPSB2_exons$transcript_name)], sep = "-")
SPSB2_exons$transcript_name_unique[grepl("TALON", SPSB2_exons$transcript_name)] <- paste("SPSB2", SPSB2_exons$transcript_name[grepl("TALON", SPSB2_exons$transcript_name)], sep = "-")
SPSB2_exons$transcript_name_unique[grepl("TALON", SPSB2_exons$transcript_name)] <- paste("TMEM222", SPSB2_exons$transcript_name[grepl("TALON", SPSB2_exons$transcript_name)], sep = "-")
SPSB2_exons$transcript_name_unique[grepl("TALON", SPSB2_exons$transcript_name)] <- paste("HLA-DQA1", SPSB2_exons$transcript_name[grepl("TALON", SPSB2_exons$transcript_name)], sep = "-")
# SPSB2_exons$transcript_name_unique[grepl("TALON", SPSB2_exons$transcript_name)] <- paste("SCGB3A2", SPSB2_exons$transcript_name[grepl("TALON", SPSB2_exons$transcript_name)], sep = "-")
SPSB2_exons$transcript_name_unique[grepl("TALON", SPSB2_exons$transcript_name)] <- paste("OAS1", SPSB2_exons$transcript_name[grepl("TALON", SPSB2_exons$transcript_name)], sep = "-")
colnames(SPSB2_exons)[7] <- "annot_transcript_name"
SPSB2_exons <- left_join(SPSB2_exons, TALON_afterqc_orf_secondpass2[,c("annot_transcript_name", "transcript_catalog")], by = "annot_transcript_name")
SPSB2_exons <- SPSB2_exons %>% dplyr::filter(transcript_name %in% c("SPSB2-205", "TALONT000636440", "TALONT000636472", "SPSB2-204"))
SPSB2_exons <- SPSB2_exons %>% filter(transcript_name %in% c("TMEM222-211", "TALONT001958422"))
# SPSB2_exons <- SPSB2_exons %>% filter(transcript_name %in% c("OAS1-201", "OAS1-202", "OAS1-203"))
SPSB2_exons <- SPSB2_exons %>% filter(transcript_name %in% c("HLA-DQA1-203", "TALONT001781255"))

# only for HLA-DQA1 with longer intron

# extract exons

pknox1_rescaled <- shorten_gaps(
  exons = SPSB2_exons, 
  introns = to_intron(SPSB2_exons, "transcript_name"), 
  group_var = "transcript_name"
)

# shorten_gaps() returns exons and introns all in one data.frame()
# let's split these for plotting 
pknox1_rescaled_exons <- pknox1_rescaled %>% dplyr::filter(type == "exon") 
pknox1_rescaled_introns <- pknox1_rescaled %>% dplyr::filter(type == "intron") 

SPSB2_base_plot <-  pknox1_rescaled_exons %>% 
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = transcript_catalog)
  ) +
  geom_intron(
    data = pknox1_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 300
  )




SPSB2_base_plot <-  SPSB2_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = annot_transcript_name
  )) +
  geom_range(
    aes(fill = transcript_catalog)
  ) +
  geom_intron(
    data = to_intron(SPSB2_exons, "annot_transcript_name"),
    aes(strand = strand), 
    arrow.min.intron.length = 3500
  )
SPSB2_base_plot + theme_bw() + 
  ggtitle("HLA-DQA1 gene") +  #theme_bw() + 
  xlab("chr6") + ylab("Isoforms") + 
  scale_fill_manual(values = cat.palette) +
  # geom_vline(xintercept = c(6872998), linetype = 1,
  #            colour = "darkgreen")+
  theme(legend.position = "bottom", legend.spacing.x  =  unit(0.01, 'cm'), legend.key.size = unit(0.25, 'cm')) 
SPSB2_base_plot + theme_bw() + 
  ggtitle("PPIL6 gene") +  #theme_bw() + 
  xlab("chr6") + ylab("Isoforms") + 
  scale_fill_manual(values = cat.palette) +
  # geom_vline(xintercept = c(6872998), linetype = 1,
  #            colour = "darkgreen")+
  theme(legend.position = "bottom", legend.spacing.x  =  unit(0.01, 'cm'), legend.key.size = unit(0.25, 'cm')) 

SPSB2_base_plot + theme_bw() + 
  ggtitle("TMEM222 gene") +  #theme_bw() + 
  xlab("chr1") + ylab("Isoforms") + 
  scale_fill_manual(values = cat.palette) +
  # geom_vline(xintercept = c(6872998), linetype = 1,
  #            colour = "darkgreen")+
  theme(legend.position = "bottom", legend.spacing.x  =  unit(0.01, 'cm'), legend.key.size = unit(0.25, 'cm')) 
SPSB2_base_plot + theme_bw() + 
  ggtitle("SPSB2 gene") +  #theme_bw() + 
  xlab("chr12") + ylab("Isoforms") + 
  scale_fill_manual(values = cat.palette) +
  geom_vline(xintercept = c(6872998), linetype = 1,
             colour = "darkgreen")+
  theme(legend.position = "top", legend.spacing.x  =  unit(0.01, 'cm'), legend.key.size = unit(0.25, 'cm')) 
SPSB2_base_plot + theme_bw() + 
  ggtitle("OAS1 gene") +  #theme_bw() + 
  xlab("chr12") + ylab("Isoforms") + 
  scale_fill_manual(values = cat.palette) +
  # geom_textvline(aes(xintercept=6872998, label = "rs11064437"), lty = 1, colour= 'red',)+
  # geom_vline(xintercept = c(112919388), linetype = 1,
  #            colour = "darkgreen")+
  theme(legend.position = "top", legend.spacing.x  =  unit(0.01, 'cm'), legend.key.size = unit(0.25, 'cm')) 
ggsave("/data/lib14/project/scLongread/Fig5F.pdf", width = 12,height = 4)
ggsave("/data/lib14/project/scLongread/Fig5H1.pdf", width = 12,height = 4)
ggsave("/data/lib14/project/scLongread/FigS8D.pdf", width = 12,height = 7)
ggsave("/data/lib14/project/scLongread/FigS11A.pdf", width = 10,height = 10)
ggsave("/data/lib14/project/scLongread/FigS13B.pdf", width = 10,height = 7)
# add_x_break <- function(plot, xval) {
#   
#   p2 <- ggplot_build(plot)
#   breaks <- p2$layout$panel_params[[1]]$x$breaks
#   breaks <- breaks[!is.na(breaks)]
#   
#   plot +
#     geom_vline(xintercept = xval) +
#     scale_x_continuous(breaks = sort(c(xval, breaks)))
# }
# 
# add_x_break(SPSB2_base_plot, )
# 
# transcript_gtf %>%
#     ggplot(aes(
#         xstart = start,
#         xend = end,
#         y = transcript_id
#     )) +
#     geom_range(
#         aes(fill = structural_category)) +
#     geom_intron(
#         data = to_intron(transcript_gtf, "transcript_id"),
#         aes(strand = strand)
#     )  + ggtitle("PICALM gene") +  theme_jh() + xlab("chr11") + ylab("Isoforms") + theme(legend.position = "bottom", legend.spacing.x  =  unit(0.01, 'cm'), axis.title.y=element_text(angle=0)) +
#   scale_fill_manual(values=c("#440154FF",  "#31688EFF", "#FDE725FF", "#35B779FF"))
# 
# 
# 
# p <- 
#   exons %>% 
#   arrange((regular_FDR)) %>%
#   mutate( transcript_id = factor(transcript_id, levels = unique(rev(transcript_id)))) %>%
#   ggplot(aes_string(
#     xstart = "start",
#     xend = "end",
#     y = "transcript_id",
#     fill = "structural_category"
#   )) +
#   
#   geom_intron(arrow.min.intron.length = 10000,
#               data = intron_data,
#               aes(strand = strand)) +
#   geom_range(height = 0.25) +
#   geom_range(data = cds) +
#   labs(y = "") +
#   #scale_fill_distiller(name = colourby, palette = "RdBu") + 
#   theme_jh()  +
#   scale_fill_viridis_d(name=NULL) + ggtitle("PICALM gene") + xlab("chr11")  + 
#   theme(legend.position = "bottom", legend.spacing.x  =  unit(0.01, 'cm'), legend.key.size = unit(0.25, 'cm')) 

ATC_nominal <- NB.nominal.sig.list$Alveolar_transitional_cells

gene_id = "ENSG00000164265"
idx <- grepl(gene_id, TALON_afterqc_orf_secondpass2$annot_gene_id)
isoform_list <- TALON_afterqc_orf_secondpass2$transcript_name_unique[idx]
head(isoform_list)
isoform_list <- TALON_afterqc_orf_secondpass2$annot_transcript_id[idx]
Nominal_combined <- do.call(rbind, NB.nominal.sig.list)
Nominal_combined_sub <- Nominal_combined %>% 
  dplyr::select(chrom,snp,pos,a0,a1)
Nominal_combined_sub <- Nominal_combined_sub[!duplicated(Nominal_combined_sub),]
Nominal_combined_sub$chrom <- paste0("chr", Nominal_combined_sub$chrom)
sjs <- read.table("/data/Choi_lung/scLongreads/TALON_workspace/test10s/NCI_lung_eIsoform_sjg_introns.tsv",
                  header = TRUE, sep = "\t")
sjs <- sjs[!duplicated(sjs),]
# if strand +
sjs$donor_start <- sjs$start-2
sjs$donor_end <- sjs$start+1
sjs$receptor_start <- sjs$stop-1
sjs$receptor_end <- sjs$stop+2
sjs$sjs <- paste(sjs$chrom, sjs$start, sjs$stop, sep = "-")

sjs$donor_start[which(sjs$strand == "-")] <- sjs$donor_start[which(sjs$strand == "-")] + 1
sjs$donor_end[which(sjs$strand == "-")] <- sjs$donor_end[which(sjs$strand == "-")] + 1
sjs$receptor_start[which(sjs$strand == "-")] <- sjs$receptor_start[which(sjs$strand == "-")] - 1
sjs$receptor_end[which(sjs$strand == "-")] <- sjs$receptor_end[which(sjs$strand == "-")] - 1
located_sjs.list <- list()
for (i in 1:nrow(Nominal_combined_sub)) {
  located_sjs <- sjs$sjs[which(( Nominal_combined_sub$chrom[i] == sjs$chrom & Nominal_combined_sub$pos[i] >= sjs$donor_start & Nominal_combined_sub$pos[i] <= sjs$donor_end) |
                                 (Nominal_combined_sub$chrom[i] == sjs$chrom & Nominal_combined_sub$pos[i] >= sjs$receptor_start & Nominal_combined_sub$pos[i] <= sjs$receptor_end) )]
  located_sjs.list[[i]] <- located_sjs
  print(i)
}
names(located_sjs.list) <- Nominal_combined_sub$snp
snps <- Nominal_combined_sub$snp
test.list <- lapply(snps, function(snp){
  if (length(located_sjs.list[[snp]])==0) {
    return(located_sjs.list[[snp]])
  }else{
    x <- located_sjs.list[[snp]]
    df <- data.frame(sjss = paste(x, collapse = ","),
                     snp = snp)
    return(df)
  }
})
tmp <- do.call(rbind,  test.list)
Nominal_combined_sjsloc <- subset(Nominal_combined, snp %in% unique(tmp$snp))
Nominal_combined_sjsloc <- Nominal_combined_sjsloc %>% group_by(phenotype_id) %>% mutate(transcript_name = Search_transcript_name2(phenotype_id), .after = phenotype_id)
Nominal_combined_sjsloc <- left_join(Nominal_combined_sjsloc, tmp, by = "snp")
length(unique(tmp$phenotype_id))
write.table(Nominal_combined_sjsloc, "/data/Choi_lung/scLongreads/DEI/isoQTLs_located_sjs_rst.txt",
            sep = "\t", row.names = FALSE,quote = FALSE)
write.table(Nominal_combined_sjsloc_final, "/data/Choi_lung/scLongreads/DEI/isoQTLs_located_sjs_rst_final.txt",
            sep = "\t", row.names = FALSE,quote = FALSE)
# # if strand -
# # change start and end orders for coding convenience 
# sjs$donor_end <- sjs$start+2
# sjs$donor_start <- sjs$start-1
# sjs$receptor_end <- sjs$stop+1
# sjs$receptor_start <- sjs$stop-2

length(set1)
write.table(set1, "/data/Choi_lung/scLongreads/TALON_workspace/test10s/NCI_lung_secondpass.isoGenes",
            row.names = FALSE, col.names = FALSE, quote = FALSE)



NB.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")
names(NB.sig.list)
Nominal_combined_sjsloc <- read.table("/data/Choi_lung/scLongreads/DEI/isoQTLs_located_sjs_rst_final.txt",
                                      sep = "\t", header = TRUE)
Nominal_combined_sjsloc <- subset(Nominal_combined_sjsloc, snp %in% unique(tmp$snp))
celltypes <- intersect(names(NB.sig.list),unique(Nominal_combined_sjsloc$Celltype))

library(reshape2)

tmp <- Nominal_combined_sjsloc[,c("snp", "Celltype")]
tmp <- tmp[!duplicated(tmp),]
sjs_snps <- as.matrix.data.frame(table(tmp$snp, tmp$Celltype))
rownames(sjs_snps) <- names(table(tmp$snp))
colnames(sjs_snps) <- names(table(tmp$Celltype))
sjs_snps <- sjs_snps[,celltypes]
library(UpSetR)
set_vars <- celltypes
main_bar_col <- "violetred4"
sets_bar_col <- "turquoise4"
matrix_col <- "slateblue4"
shade_col <- "wheat4"

library(ComplexHeatmap)
m <- make_comb_mat(sjs_snps)
cs = comb_size(m)
intersect(which(comb_degree(m) == 1), which(order(comb_degree(m), -cs) %in% 1:40))
m = m[union(which(comb_degree(m) == 1), which(order(comb_degree(m), -cs) %in% 1:40))] 
m = m[which(order(comb_degree(m), -cs) %in% 1:40)] 
m = m[which(comb_degree(m) == 1)]
ss = set_size(m)
cs = comb_size(m)
# pdf("/data/lib14/project/AT2_co_CS/peak_by_celltype.pdf", width = 10)
pdf("/data/Choi_lung/lbl/SHARE-seq/correction/peak_by_celltype.pdf", width = 10 )
ht = UpSet(m, 
           set_order = order(ss),pt_size =  unit(3, "mm"),comb_col = "slateblue4",lwd = 3,
           comb_order = order(comb_degree(m), -cs),
           top_annotation = HeatmapAnnotation(
             "Numbers of cell type specific peaks" = anno_barplot(cs, 
                                                                  ylim = c(0, max(cs)*1.1),
                                                                  border = FALSE, 
                                                                  gp = gpar(fill = c(rep("turquoise4",23), rep("grey",7)),fontfamily = "Times New Roman"), 
                                                                  height = unit(5, "cm"), width = 6
             ), 
             annotation_name_side = "left",
             annotation_name_rot = 0
           ),
           right_annotation = rowAnnotation(
             "Total numbers of peaks" = anno_barplot(ss, 
                                                     baseline = 0,
                                                     # axis_param = list(
                                                     #    at = c(0, -50000, -100000),
                                                     #    labels = c(0, 50000, 100000),
                                                     #   labels_rot = 0),
                                                     border = FALSE, 
                                                     gp = gpar(fill = "violetred4",fontfamily = "Times New Roman"), 
                                                     width = unit(4, "cm")
             ) 
             # set_name = anno_text(set_name(m), 
             #                       location = -12.7, 
             #                       just = "right"
             #                       #width = max_text_width(set_name(m)) + unit(2, "mm")
             #                     )
           ), 
           #right_annotation = NULL,
           show_row_names = TRUE)
ht = draw(ht)
od = column_order(ht)
decorate_annotation("Numbers of cell type specific peaks", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "black",fontfamily = "Times New Roman"), rot = 45)
})
dev.off()

Nominal_combined_sub <- Nominal_combined %>% 
  dplyr::select(snp,phenotype_id)
Nominal_combined_sub <- Nominal_combined_sub[!duplicated(Nominal_combined_sub),]
snp.pruned <- read.table("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/prune.prune.in")
Nominal_combined_sub <- subset(Nominal_combined_sub, snp %in% snp.pruned$V1)
tmp <- Nominal_combined_sjsloc[,c("snp", "phenotype_id")]
tmp <- tmp[!duplicated(tmp),]
hist(as.data.frame(table(Nominal_combined_sub$snp))$Freq)
hist(as.data.frame(table(tmp$snp))$Freq)

df1 <- as.data.frame(table(Nominal_combined_sub$snp))
df2 <- as.data.frame(table(tmp$snp))
df1$Type <- "ld_pruned_isoQTLs"
df2$Type <- "sjs_isoQTLs"
fisher.test(matrix(c((60-33), 33, (3044-2697), 2697),
                   nrow = 2, byrow = TRUE))
library(hrbrthemes)
library(viridis)
library(reshape2)
df <- rbind(df1, df2)
ggplot(subset(df, Freq < 6), aes(x=Freq, color=Type, fill=Type)) +
  geom_density(alpha = 0.6) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) + theme_bw()
  
ggplot(df, aes(x=Freq, color=Type, fill=Type)) +
  geom_histogram(aes(y = ..density..), binwidth = 1, alpha = 0.6)
df_merged <- as.matrix.data.frame(table(df$Freq, df$Type))
colnames(df_merged) <- c("ld_pruned_isoQTLs", "sjs_isoQTLs")
df_merged <- as.data.frame(df_merged)
df_merged$eIsoform_num <- as.numeric( names(table(df$Freq)))
df_merged$ld_pruned_isoQTLs_frac <- df_merged$ld_pruned_isoQTLs/3044
df_merged$sjs_isoQTLs_frac <- df_merged$sjs_isoQTLs/179

df_long <- melt(df_merged[,c(3,4,5)], id.vars = "eIsoform_num")
p <- ggplot(df_long, aes(x=eIsoform_num, y = value, fill=variable)) +
  geom_bar(position="dodge", stat="identity", alpha = 0.6, width=1) + scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) + ylab("Fraction") + xlab("Number of eIsoforms associated with SNP")+
  theme(axis.text.x = element_text(color = "black", size = 12,face = "plain"),
        axis.text.y = element_text(color = "black", size = 12,face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1))
ggsave("/data/lib14/project/scLongread/Fig5E.pdf", width = 8,height = 7)
print(p)



# isoQTL distribution on chromosome
gff = as.data.frame(rtracklayer::import('/data/Choi_lung/scLongreads/TALON_workspace/test10s/NCI_lung_secondpass.isoforms_w_known.gtf'))
Nominal_combined <- do.call(rbind, NB.nominal.sig.list)
colnames(Nominal_combined)[6] <- "transcript_id"
Nominal_combined <- left_join(Nominal_combined, transcript_gff[,c("transcript_id", "start", "end")], by = "transcript_id")
Nominal_combined_sub <- Nominal_combined %>% 
  dplyr::select(chrom,snp,pos,a0,a1, tss, tss_distance, start, end, Celltype)
Nominal_combined_sub <- Nominal_combined_sub[!duplicated(Nominal_combined_sub),]
idx <- which((Nominal_combined_sub$pos >= Nominal_combined_sub$start & Nominal_combined_sub$pos <= Nominal_combined_sub$end))
tmp <- Nominal_combined_sub[idx,]
hist(abs(tmp$tss_distance))
hist(abs(subset(tmp, Celltype == "AT2")$tss_distance))
hist(abs(subset(tmp, Celltype == "AT1")$tss_distance))
hist(abs(subset(tmp, Celltype == "Alveolar_macrophages")$tss_distance))
GTEx_sGene <- GTEx_sGene[which(GTEx_sGene$qval < 0.05),]
idx <- which((GTEx_sGene$variant_pos >= GTEx_sGene$gene_start & GTEx_sGene$variant_pos <= GTEx_sGene$gene_end))
hist(abs(GTEx_sGene[idx,]$tss_distance))
