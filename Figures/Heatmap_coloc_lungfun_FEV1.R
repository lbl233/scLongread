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

output.ratio <- readRDS("/data/Choi_lung/scLongreads/colocalization/lung_function/Colocalization_output_table_fev1.rds")
NB.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")
NB.sig.list <- lapply(names(NB.sig.list), function(x){
  rst <- NB.sig.list[[x]]
  rst$Celltype <- x
  return(rst)
})
isoQTL_sig <- do.call(rbind,NB.sig.list)

isoQTL_sig$sig_iso_ct <- paste(isoQTL_sig$phenotype_id, isoQTL_sig$Celltype, sep = "-")
output.ratio$sig_iso_ct <- paste(output.ratio$transcript_id, output.ratio$Celltype, sep = "-")
output.ratio <- subset(output.ratio, sig_iso_ct %in% isoQTL_sig$sig_iso_ct)

transcripts <- unique(output.ratio$transcript_name)
celltypes <- unique(output.ratio$Celltype)
celltype_levels <- levels(lr$Celltype)
celltypes <- intersect(gsub(" ", "_", celltype_levels), celltypes)
output.ratio <- output.ratio %>% group_by(transcript_id) %>% mutate(gene_name = Search_transcript_gene_name(transcript_id))
output.ratio <- output.ratio %>% group_by(transcript_id) %>% mutate(gene_id = Search_transcript_gene_id(transcript_id))
colnames(lung_fun_assoc )[21] <- "leadSNP"
output.ratio <- left_join(output.ratio, lung_fun_assoc[,c("leadSNP", "chr","pos")], by = "leadSNP")
coloc_rst_PPH4 <- output.ratio[,c(6,8,9,10)]
write.table(coloc_rst_PPH4, file = "/data/Choi_lung/scLongreads/colocalization/TableS5.txt", sep = "\t", row.names = FALSE, quote = FALSE)
coloc_rst_PPH4$Celltype <- factor(coloc_rst_PPH4$Celltype, levels = celltypes)
library(reshape2)
# coloc_rst_PPH4$value[is.na(coloc_rst_PPH4$value)] <- 0
hm <- ggplot(data = coloc_rst_PPH4, aes(x = Celltype, y = transcript_name, fill = as.numeric(PP.H4.abf))) + geom_tile() + 
  scale_fill_distiller(name = "Legend title", palette = "Reds", direction = 1, na.value = "white") +
  geom_text(aes(label = round(PP.H4.abf, 2)), size = 3) +
  scale_x_discrete(breaks = unique(coloc_rst_PPH4$Celltype), labels = unique(coloc_rst_PPH4$Celltype)) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 15), legend.key.size = unit(1,"cm"),
        legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))
hm
# Substract legend for heatmap
tmp <- ggplot_gtable(ggplot_build(hm))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
# Remove legend from heatmap
hm.clean <- hm +
  theme(
    # axis.title.y = element_blank(), axis.text.y = element_blank(),
    #     axis.ticks.y = element_blank(), axis.title.x = element_blank(),
    #     axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text  = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none")
hm.clean
tmp <- output.ratio[,c("Celltype", "transcript_id")]
tmp <- tmp[!duplicated(tmp),]
tmp$Celltype <- factor(tmp$Celltype, levels = celltypes)
x.axis.bar.df <- as.data.frame(table(tmp$Celltype))
x.axis.bar.df$Var1 <- factor(x.axis.bar.df$Var1, levels = celltypes)
x.axis.bar.df$Category <- c(rep("Epithelial",7), rep("Immune", 3), rep("Endothelial", 2))
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
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
pdf("/data/lib14/project/scLongread/Fig7B.pdf", width = 7.8 , height = 10)
print(grid.arrange(arrangeGrob(legend, bp.x.clean, ncol = 2,widths = c(10,20)), hm.clean,  nrow = 2, ncol = 1, heights = c(10, 80)))
dev.off()


transcripts_genes <- c("HLA-DRB5","SECISBP2L",
                 "HLA-DRB1","CHMP3","RSPH4A",
                 "PPIL6","HLA-A","PSMA4","HSPA1B")

write.table(output.ratio, file = "/data/Choi_lung/scLongreads/colocalization/lung_function/Colocalization_output_table.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(output.fev1, file = "/data/Choi_lung/scLongreads/colocalization/lung_function/Colocalization_output_table_fev1.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)


shrine_coloc <- read.table("/data/Choi_lung/scLongreads/colocalization/lung_function/Shrine_coloc.txt", header = TRUE, sep = "\t")
unique(shrine_coloc$gene[grepl("Lung", shrine_coloc$tissue)])
unique(output.ratio$gene_name)
length(setdiff(unique(output.ratio$gene_name), unique(shrine_coloc$gene[grepl("Lung", shrine_coloc$tissue)])))
intersect(unique(output.ratio$gene_name), unique(shrine_coloc$gene))
