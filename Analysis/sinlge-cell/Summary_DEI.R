library(stringr)
library(stringi)
library(dplyr)
library(ggplot2)
library(Seurat)

files <- list.files("/data/Choi_lung/Elelta/eQTLProjects/BolunTask/PostIsoformFiltering/WholeGroupCorrected/Outputfeb19/unfiltered/", pattern = ".csv")


celltypes <- str_split_fixed(files, "_", 4)[,4]
celltypes <- str_split_fixed(celltypes, "\\.", 2)[,1]
celltypes <- gsub(" ", "_", celltypes)

DEI.list <- list()
Sig_DEI <- NULL
for (i in 1:37) {
  DEI_rst <- read.csv(paste0("/data/Choi_lung/Elelta/eQTLProjects/BolunTask/PostIsoformFiltering/WholeGroupCorrected/Outputfeb19/unfiltered/",files[i]))
  DEI.list[[celltypes[i]]] <- DEI_rst
  Sig_DEI <- c(length(which(DEI_rst$padj < 0.05 & abs(DEI_rst$log2FoldChange) > 1)), Sig_DEI)
}
names(Sig_DEI) <- celltypes

DEI_sig.list <- list()
DEI_unique.list <- list()
for (i in 1:37) {
  DEI_df <- DEI.list[[i]]
  DEI_sub <- DEI_df[which(DEI_df$padj < 0.05 & DEI_df$log2FoldChange > 1),]
  DEI_sub$Celltype <- celltypes[i]
  DEI_sig.list[[i]] <- DEI_sub
  tmp <- do.call(rbind, DEI.list[-i])
  DEI_other <- unique(tmp$X[which(tmp$padj < 0.05 & tmp$log2FoldChange > 0)])
  DEI_unique <- setdiff(DEI_sub$X, DEI_other)
  DEI_unique.list[[i]] <- DEI_unique
}
names(DEI_sig.list) <- celltypes
names(DEI_unique.list) <- celltypes
save(list = ls(), file = "/data/Choi_lung/scLongreads/DEI/DEI_workspace.RData")
df_DEI <- data.frame(Celltype = celltypes,
                     Number = sapply(DEI_sig.list, nrow))
ct_levels <- c(
  'Secretory transitional cells', 'Club', 'AT2', 'Multiciliated', 'AT1', 'Alveolar transitional cells', 'Goblet', 'Basal', 'Neuroendocrine',  # Epithelial
  'Lymphatic EC', 'EC aerocyte capillary', 'EC arterial', 'EC venous systemic', 'EC venous pulmonary', 'EC general capillary',  # Endothelial
  'CD8 T cells', 'Alveolar macrophages', 'NK T cells', 'NK cells', 'Interstitial macrophages', 'CD4 T cells', 'Classical monocyte', 'Alveolar macrophages MT', 'DC2', 'Plasmacytoid DCs', 
  'T cell proliferating', 'Non_classical monocytes', 'Alveolar macrophages CCL3', 'Mast cells', 'Alveolar Mph proliferating', 'B cells', 'Plasma cells',  # Immune
  'SMC', 'Adventitial fibroblasts', 'Mesothelium', 'Myofibroblast', 'Alveolar fibroblasts'  # Stromal
)
ct_levels <- gsub(" ", "_", levels(lr.isoform.sub$Celltype))
df_DEI$Celltype <- factor(df_DEI$Celltype, levels = ct_levels)
fc <- colorRampPalette(c("#FA877533", "#610C04"))
fc_yellow <- colorRampPalette(c("#8B4000", "#FFD70033"))
fc_purple <- colorRampPalette(c("lightpink","#4B0082"))
fc_blue <- colorRampPalette(c("lightblue", "darkblue"))
fc_yellow(9)
cellcolors <- c(fc_yellow(9)[c(3,2,1,4,8,9,7,5,6)], fc(17)[c(1,17,15,16,14,2,13,4,11,
                                                             3,12,5,
                                                             6,7,9,10,8)],
                fc_purple(6)[c(4,2,6,3,5,1)], fc_blue(5)[c(5,1,2,4,3)])




cols <- c('darkseagreen1' , 'green4', 'green1', 'seagreen4','olivedrab', 'palegreen2', 'limegreen', 'darkgreen', 'springgreen4',  # Epithelial (9 colors)
          'royalblue1', 'dodgerblue4', 'steelblue4', 'blue', 'aquamarine4', 'lightblue',  # Endothelial (6 colors)
          'purple1', 'orchid1','mediumpurple3','mediumorchid2','darkviolet','magenta','orchid4', 'lavenderblush3','darkorchid4', 'slateblue', 'blueviolet', 
          'deeppink3', 'mediumvioletred', 'hotpink3', 'mediumorchid1', 'violetred1', 'plum2', 'plum3',  # Immune (17 colors)
          'tan', 'lightsalmon', 'sienna', 'gold', 'goldenrod3'  # Stromal (5 colors)
)
ggplot(df_DEI, aes(x = Celltype, y = Number, fill = Celltype)) + 
  geom_col(position = 'dodge') +
  geom_text(aes(label = Number), position = position_dodge(width = 0.9), vjust = -0.25, size = 5) + 
  coord_cartesian(ylim = c(0, max(df_DEI$Number) * 1.1)) + 
  ylab('Number of DEIs') + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10),
        legend.position="top",legend.justification="right") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = cellcolors) + NoLegend()
tmp.list <- lapply(DEI_sig.list, function(x) {x$X})
length(unique(Reduce(c, tmp.list)))



df_DEI_unique <- data.frame(Celltype = celltypes,
                            Number = sapply(DEI_unique.list, length))
df_DEI_unique$Celltype <- factor(df_DEI_unique$Celltype, levels = ct_levels)
ggplot(df_DEI_unique, aes(x = Celltype, y = Number, fill = Celltype)) + 
  geom_col(position = 'dodge') +
  geom_text(aes(label = Number), position = position_dodge(width = 0.9), vjust = -0.25, size = 5) + 
  coord_cartesian(ylim = c(0, max(df_DEI_unique$Number) * 1.1)) + 
  ylab('Number of Unique Isoforms') + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color= cols),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10),
        legend.position="top",legend.justification="right") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = cols) + NoLegend()

length(unique(Reduce(c, DEI_unique.list)))
sum(df_DEI_unique$Number)

# For visualization by logFC
Top_DEI <- lapply(DEI_sig.list, function(x){
  x <- x[order(x$log2FoldChange, decreasing = TRUE),]
  x$X[c(1:3)]
})
Top_DEI <- Top_DEI[gsub(" ", "_", ct_levels)]
length(unique(Reduce(c,Top_DEI)))
# Top DEI For visualization by padj
Top_DEI <- lapply(DEI_sig.list, function(x){
  x <- x[order(x$padj),]
  x$X[c(1:3)]
})
Top_DEI <- Top_DEI[gsub(" ", "_", ct_levels)]
length(unique(Reduce(c,Top_DEI)))


# matching with eIsoforms
DEI_visual <- lapply(DEI_unique.list, function(x){
  if (length(x) < 3) {
    return(x)
  }else{return(x[c(1:3)])}
})
DEI_visual <- DEI_visual[ct_levels]
lr <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final_v2_wRNAexpr.RDS")
lr.isoform <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr.sc.data.Isoform.22samples.combined.rds")

DefaultAssay(lr.isoform) <- "Isoform"
# Subset the original isoform expression profiles to match your Seurat object
lr.isoform.sub <- subset(lr.isoform, cells = colnames(lr))
lr.isoform.sub$Sample <- lr$Sample
lr.isoform.sub$Celltype <- lr$Celltype
# Normalize the Isoform assay data within the Seurat object
lr.isoform.sub <- NormalizeData(lr.isoform.sub, assay = "Isoform", normalization.method = "LogNormalize", scale.factor = 10000)
saveRDS(lr.isoform.sub, file = "/data/Choi_lung/scLongreads/Seurat/final/lr_final_w_isoform_expr.RDS")
lr.isoform.sub <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final_w_isoform_expr.RDS")
# Check that the normalization has been applied
print(lr.isoform.sub@assays$Isoform)
lr.isoform.sub$Celltype <- factor(lr.isoform.sub$Celltype, levels = ct_levels)
# Step 5: Generate the DotPlot with the reordered markers
pdf("/data/Choi_lung/scLongreads/DEI/Dotplot_unique_DEI.pdf", width = 30)
DotPlot(lr.isoform.sub, assay = "Isoform", features = DEI_visual, 
        cols = c('lightgrey', 'blue'), scale = TRUE, scale.by = 'radius', 
        group.by = 'Celltype') + RotatedAxis()
dev.off()
pdf("/data/Choi_lung/scLongreads/DEI/Dotplot_top_DEI_bypvalue.pdf", width = 30)
DotPlot(lr.isoform.sub, assay = "Isoform", features = unique(Reduce(c,Top_DEI)), 
        cols = c('lightgrey', 'blue'), scale = TRUE, scale.by = 'radius', 
        group.by = 'Celltype') + RotatedAxis()
dev.off()


files <- list.files("/Users/lib14/Downloads/For_Bolun/DEG/Unfiltered/", pattern = ".csv")

files <- files[-c(1,2)]
celltypes <- str_split_fixed(files, "_", 4)[,4]
celltypes <- str_split_fixed(celltypes, "\\.", 2)[,1]
celltypes <- gsub(" ", "_", celltypes)

DEG.list <- list()
Sig_DEG <- NULL
for (i in 1:37) {
  DEI_rst <- read.csv(paste0("/Users/lib14/Downloads/For_Bolun/DEG/Unfiltered/",files[i]))
  
  DEG.list[[celltypes[i]]] <- DEI_rst
  Sig_DEG <- c(length(which(DEI_rst$padj < 0.05 & DEI_rst$log2FoldChange > 1)), Sig_DEG)
}
names(Sig_DEG) <- celltypes

DEG_sig.list <- list()
DEG_sigpadj.list <- list()
DEG_unique.list <- list()
for (i in 1:37) {
  DEI_df <- DEG.list[[i]]
  DEI_sub <- DEI_df[which(DEI_df$padj < 0.05 & DEI_df$log2FoldChange > 0),]
  DEI_sub$Celltype <- celltypes[i]
  DEG_sigpadj.list[[i]] <- DEI_sub
  DEI_sub <- DEI_df[which(DEI_df$padj < 0.05 & DEI_df$log2FoldChange > 1),]
  DEG_sig.list[[i]] <- DEI_sub
  tmp <- do.call(rbind, DEG.list[-i])
  DEI_other <- unique(tmp$X[which(tmp$padj < 0.05 & tmp$log2FoldChange > 0)])
  DEG_unique <- setdiff(DEI_sub$X, DEI_other)
  DEG_unique.list[[i]] <- DEG_unique
}
names(DEG_sig.list) <- celltypes
names(DEG_sigpadj.list) <- celltypes
names(DEG_unique.list) <- celltypes
`%nin%` <- Negate(`%in%`)
DEI_sig_onlyiso.list <- list()
for (celltype in celltypes) {
  DEI_sig <- DEI_sig.list[[celltype]]
  DEI_sig$gene_name <- TALON_afterqc_orf_secondpass2[DEI_sig$X,]$annot_gene_name
  DEG_sigpadj <- DEG_sigpadj.list[[celltype]]
  DEI_sig_onlyiso <- subset(DEI_sig, gene_name %nin% DEG_sigpadj$X)
  DEI_sig_onlyiso$Celltype <- celltype
  DEI_sig_onlyiso.list[[celltype]] <- DEI_sig_onlyiso
}
DEI_sig_onlyiso <- do.call(rbind, DEI_sig_onlyiso.list)
length(unique(DEI_sig_onlyiso$X))
S100A4_isoforms <- TALON_afterqc_orf_secondpass2$transcript_name_unique[grep("S100A4", TALON_afterqc_orf_secondpass2$transcript_name_unique)]
DotPlot(lr.isoform.sub, assay = "Isoform", features = "S100A4-201", 
        cols = c('lightgrey', 'blue'), scale = TRUE, scale.by = 'radius', 
        group.by = 'Celltype') + RotatedAxis()
DotPlot(lr.isoform.sub, assay = "Isoform", features = S100A4_isoforms, 
        cols = c('lightgrey', 'blue'), scale = TRUE, scale.by = 'radius', 
        group.by = 'Celltype') + RotatedAxis()
DotPlot(lr, assay = "RNA", features = "S100A4", 
        cols = c('lightgrey', 'blue'), scale = TRUE, scale.by = 'radius', 
        group.by = 'Celltype') + RotatedAxis()
RPL18_isoforms <- TALON_afterqc_orf_secondpass2$transcript_name_unique[grep("HSD11B1L-", TALON_afterqc_orf_secondpass2$transcript_name_unique)]
DotPlot(lr.isoform.sub, assay = "Isoform", features = "RPL18-208", 
        cols = c('lightgrey', 'blue'), scale = TRUE, scale.by = 'radius', 
        group.by = 'Celltype') + RotatedAxis()
DotPlot(lr.isoform.sub, assay = "Isoform", features = RPL18_isoforms, 
        cols = c('lightgrey', 'blue'), scale = TRUE, scale.by = 'radius', 
        group.by = 'Celltype') + RotatedAxis()

DotPlot(lr, assay = "RNA", features = "RPL18", 
        cols = c('lightgrey', 'blue'), scale = TRUE, scale.by = 'radius', 
        group.by = 'Celltype') + RotatedAxis()
SRGN_isoforms <- TALON_afterqc_orf_secondpass2$transcript_name_unique[grep("SRGN-", TALON_afterqc_orf_secondpass2$transcript_name_unique)]
DotPlot(lr.isoform.sub, assay = "Isoform", features = SRGN_isoforms, 
        cols = c('lightgrey', 'blue'), scale = TRUE, scale.by = 'radius', 
        group.by = 'Celltype') + RotatedAxis()
df_DEI_onlyiso <- data.frame(Celltype = celltypes,
                             Number = sapply(DEI_sig_onlyiso.list, nrow))
ct_levels <- c(
  'Secretory transitional cells', 'Club', 'AT2', 'Multiciliated', 'AT1', 'Alveolar transitional cells', 'Goblet', 'Basal', 'Neuroendocrine',  # Epithelial
  'Lymphatic EC', 'EC aerocyte capillary', 'EC arterial', 'EC venous systemic', 'EC venous pulmonary', 'EC general capillary',  # Endothelial
  'CD8 T cells', 'Alveolar macrophages', 'NK T cells', 'NK cells', 'Interstitial macrophages', 'CD4 T cells', 'Classical monocyte', 'Alveolar macrophages MT', 'DC2', 'Plasmacytoid DCs', 
  'T cell proliferating', 'Non_classical monocytes', 'Alveolar macrophages CCL3', 'Mast cells', 'Alveolar Mph proliferating', 'B cells', 'Plasma cells',  # Immune
  'SMC', 'Adventitial fibroblasts', 'Mesothelium', 'Myofibroblast', 'Alveolar fibroblasts'  # Stromal
)
ct_levels <- gsub(" ", "_", ct_levels)
df_DEI_onlyiso$Celltype <- factor(df_DEI$Celltype, levels = ct_levels)
cols <- c('darkseagreen1' , 'green4', 'green1', 'seagreen4','olivedrab', 'palegreen2', 'limegreen', 'darkgreen', 'springgreen4',  # Epithelial (9 colors)
          'royalblue1', 'dodgerblue4', 'steelblue4', 'blue', 'aquamarine4', 'lightblue',  # Endothelial (6 colors)
          'purple1', 'orchid1','mediumpurple3','mediumorchid2','darkviolet','magenta','orchid4', 'lavenderblush3','darkorchid4', 'slateblue', 'blueviolet', 
          'deeppink3', 'mediumvioletred', 'hotpink3', 'mediumorchid1', 'violetred1', 'plum2', 'plum3',  # Immune (17 colors)
          'tan', 'lightsalmon', 'sienna', 'gold', 'goldenrod3'  # Stromal (5 colors)
)
ggplot(df_DEI_onlyiso, aes(x = Celltype, y = Number, fill = Celltype)) + 
  geom_col(position = 'dodge') +
  geom_text(aes(label = Number), position = position_dodge(width = 0.9), vjust = -0.25, size = 5) + 
  coord_cartesian(ylim = c(0, max(df_DEI_onlyiso$Number) * 1.1)) + 
  ylab('Number of Unique Isoforms') + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color= cols),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10),
        legend.position="top",legend.justification="right") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = cols) + NoLegend()



for (i in 1:37) {
  DEI_df <- DEI_sig.list[[i]]
  DEI_df$Celltype <- names(DEI_sig.list)[i]
  DEI_sig.list[[i]] <- DEI_df
}
DEI_sig <- Reduce(rbind,DEI_sig.list)
colnames(TALON_afterqc_orf_secondpass)
TALON_afterqc_orf_secondpass[,c("annot_gene_id", "annot_transcript_id", "annot_gene_name", "transcript_name_unique")]
colnames(DEI_sig)[1] <- "transcript_name_unique"
tmp <- left_join(DEI_sig, TALON_afterqc_orf_secondpass[,c("annot_gene_id", "annot_transcript_id", "annot_gene_name", "transcript_name_unique")], 
                 by = "transcript_name_unique")

files <- list.files("/data/Choi_lung/Elelta/eQTLProjects/BolunTask/PostIsoformFiltering/Subgroup/results5/unfilteredresults/", pattern = ".csv")
celltypes <- str_split_fixed(files, "_", 4)[,4]
celltypes <- str_split_fixed(celltypes, "\\.", 2)[,1]
celltypes <- gsub(" ", "_", celltypes)

DEI.subgroup.list <- list()

for (i in 1:37) {
  DEI_rst <- read.csv(paste0("/data/Choi_lung/Elelta/eQTLProjects/BolunTask/PostIsoformFiltering/Subgroup/results5/unfilteredresults/",files[i]))
  DEI.subgroup.list[[celltypes[i]]] <- DEI_rst
}


DEI.subgroup_sig.list <- list()
for (i in 1:37) {
  DEI_df <- DEI.subgroup.list[[i]]
  DEI_sub <- DEI_df[which(DEI_df$padj < 0.05 & DEI_df$log2FoldChange > 1),]
  DEI.subgroup_sig.list[[i]] <- DEI_sub
}
names(DEI.subgroup_sig.list) <- celltypes
DEI.combined_sig.list <- list()
for (i in 1:37) {
  DEI_whole <- DEI_sig.list[[i]]
  DEI_subgroup<- DEI.subgroup_sig.list[[i]]
  DEI_combined <- left_join(DEI_whole, DEI_subgroup, by = "X", suffix = c("whole","sub"))
  DEI_combined$Celltype <- celltypes[i]
  DEI_combined <- DEI_combined[!is.na(DEI_combined$padjsub),]
  DEI.combined_sig.list[[i]] <- DEI_combined
}
names(DEI.combined_sig.list) <- celltypes



df_DEI <- data.frame(Celltype = celltypes,
                     Number = sapply(DEI.subgroup_sig.list, nrow))
ct_levels <- c(
  'Secretory transitional cells', 'Club', 'AT2', 'Multiciliated', 'AT1', 'Alveolar transitional cells', 'Goblet', 'Basal', 'Neuroendocrine',  # Epithelial
  'Lymphatic EC', 'EC aerocyte capillary', 'EC arterial', 'EC venous systemic', 'EC venous pulmonary', 'EC general capillary',  # Endothelial
  'CD8 T cells', 'Alveolar macrophages', 'NK T cells', 'NK cells', 'Interstitial macrophages', 'CD4 T cells', 'Classical monocyte', 'Alveolar macrophages MT', 'DC2', 'Plasmacytoid DCs', 
  'T cell proliferating', 'Non_classical monocytes', 'Alveolar macrophages CCL3', 'Mast cells', 'Alveolar Mph proliferating', 'B cells', 'Plasma cells',  # Immune
  'SMC', 'Adventitial fibroblasts', 'Mesothelium', 'Myofibroblast', 'Alveolar fibroblasts'  # Stromal
)
ct_levels <- gsub(" ", "_", ct_levels)
df_DEI$Celltype <- factor(df_DEI$Celltype, levels = ct_levels)
cols <- c('darkseagreen1' , 'green4', 'green1', 'seagreen4','olivedrab', 'palegreen2', 'limegreen', 'darkgreen', 'springgreen4',  # Epithelial (9 colors)
          'royalblue1', 'dodgerblue4', 'steelblue4', 'blue', 'aquamarine4', 'lightblue',  # Endothelial (6 colors)
          'purple1', 'orchid1','mediumpurple3','mediumorchid2','darkviolet','magenta','orchid4', 'lavenderblush3','darkorchid4', 'slateblue', 'blueviolet', 
          'deeppink3', 'mediumvioletred', 'hotpink3', 'mediumorchid1', 'violetred1', 'plum2', 'plum3',  # Immune (17 colors)
          'tan', 'lightsalmon', 'sienna', 'gold', 'goldenrod3'  # Stromal (5 colors)
)
ggplot(df_DEI, aes(x = Celltype, y = Number, fill = Celltype)) + 
  geom_col(position = 'dodge') +
  geom_text(aes(label = Number), position = position_dodge(width = 0.9), vjust = -0.25, size = 5) + 
  coord_cartesian(ylim = c(0, max(df_DEI$Number) * 1.1)) + 
  ylab('Number of DEIs') + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10),
        legend.position="top",legend.justification="right") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = cols) + NoLegend()


# both in whole group and subgroup
df_DEI <- data.frame(Celltype = celltypes,
                     Number = sapply(DEI.combined_sig.list, nrow))
ct_levels <- c(
  'Secretory transitional cells', 'Club', 'AT2', 'Multiciliated', 'AT1', 'Alveolar transitional cells', 'Goblet', 'Basal', 'Neuroendocrine',  # Epithelial
  'Lymphatic EC', 'EC aerocyte capillary', 'EC arterial', 'EC venous systemic', 'EC venous pulmonary', 'EC general capillary',  # Endothelial
  'CD8 T cells', 'Alveolar macrophages', 'NK T cells', 'NK cells', 'Interstitial macrophages', 'CD4 T cells', 'Classical monocyte', 'Alveolar macrophages MT', 'DC2', 'Plasmacytoid DCs', 
  'T cell proliferating', 'Non_classical monocytes', 'Alveolar macrophages CCL3', 'Mast cells', 'Alveolar Mph proliferating', 'B cells', 'Plasma cells',  # Immune
  'SMC', 'Adventitial fibroblasts', 'Mesothelium', 'Myofibroblast', 'Alveolar fibroblasts'  # Stromal
)
ct_levels <- gsub(" ", "_", ct_levels)
df_DEI$Celltype <- factor(df_DEI$Celltype, levels = ct_levels)
cols <- c('darkseagreen1' , 'green4', 'green1', 'seagreen4','olivedrab', 'palegreen2', 'limegreen', 'darkgreen', 'springgreen4',  # Epithelial (9 colors)
          'royalblue1', 'dodgerblue4', 'steelblue4', 'blue', 'aquamarine4', 'lightblue',  # Endothelial (6 colors)
          'purple1', 'orchid1','mediumpurple3','mediumorchid2','darkviolet','magenta','orchid4', 'lavenderblush3','darkorchid4', 'slateblue', 'blueviolet', 
          'deeppink3', 'mediumvioletred', 'hotpink3', 'mediumorchid1', 'violetred1', 'plum2', 'plum3',  # Immune (17 colors)
          'tan', 'lightsalmon', 'sienna', 'gold', 'goldenrod3'  # Stromal (5 colors)
)
ggplot(df_DEI, aes(x = Celltype, y = Number, fill = Celltype)) + 
  geom_col(position = 'dodge') +
  geom_text(aes(label = Number), position = position_dodge(width = 0.9), vjust = -0.25, size = 5) + 
  coord_cartesian(ylim = c(0, max(df_DEI$Number) * 1.1)) + 
  ylab('Number of DEIs') + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10),
        legend.position="top",legend.justification="right") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = cols) + NoLegend()


DEI.combined_sig <- do.call(rbind, DEI.combined_sig.list)
length(unique(DEI.combined_sig$X))
length(unique(DEI.combined_sig$gene))
length(unique(DEG_sig$X))
idx1 <- which(TALON_afterqc_orf_secondpass$structural_category == "fusion")
TALON_afterqc_orf_secondpass$transcript_catalog <- TALON_afterqc_orf_secondpass$transcript_novelty
TALON_afterqc_orf_secondpass$transcript_catalog[idx1] <- "fusion"
TALON_afterqc_orf_secondpass$annot_gene_name[idx1] <- TALON_afterqc_orf_secondpass$gene[idx1]
DEI.combined_sig$gene <- TALON_afterqc_orf_secondpass[DEI.combined_sig$X,]$annot_gene_name

df_DEI_info <- DEI.combined_sig %>% select(X, gene) %>% distinct()
df <- as.data.frame(table(df_DEI_info$gene))
plot(df$Freq)
genes_in_gtf <- gtf %>% filter(type == "gene") %>% select(gene_name)
coding_genes_in_gtf <- gtf %>% filter(type == "gene" & gene_type == "protein_coding") %>% select(gene_name)
DEG_sig$X[which(!(DEG_sig$X %in% genes_in_gtf$gene_name))] <- str_split_fixed(DEG_sig$X[which(!(DEG_sig$X %in% genes_in_gtf$gene_name))], "\\.", 2)[,1]
DEG_sigpadj$X[which(!(DEG_sigpadj$X %in% genes_in_gtf$gene_name))] <- str_split_fixed(DEG_sigpadj$X[which(!(DEG_sigpadj$X %in% genes_in_gtf$gene_name))], "\\.", 2)[,1]
setdiff(unique(DEI.combined_sig$gene), unique(DEG_sig$X))
setdiff(unique(DEG_sig$X), unique(DEI.combined_sig$gene))
stringent_isoDEG <- setdiff(unique(DEI.combined_sig$gene), unique(DEG_sigpadj$X))
stringent_isoDEG_coding <- stringent_isoDEG[which(stringent_isoDEG %in% coding_genes_in_gtf$gene_name)] # 49
stringent_isoDEG[!grepl("^A[CLFP]",stringent_isoDEG)]
View(DEI.combined_sig[which(DEI.combined_sig$gene %in% stringent_isoDEG),])
DEI_sig <- do.call(rbind, DEI_sig.list)
DEG_sig <- do.call(rbind, DEG_sig.list)
DEG_sigpadj <- do.call(rbind, DEG_sigpadj.list)
DEG_sigpadj$gene_ct <- paste(DEG_sigpadj$X, DEG_sigpadj$Celltype, sep = ":")
DEI_sig$gene_id <- TALON_afterqc_orf_secondpass[DEI_sig$X,]$annot_gene_name
DEI_sig$gene_ct <- paste(DEI_sig$gene_id, DEI_sig$Celltype, sep = ":")
tmp <- left_join(DEI_sig, DEG_sigpadj, by = "gene_ct")
tmp <- tmp[which(tmp$log2FoldChange.x > tmp$log2FoldChange.y),]
tmp$delta_FC <- 2^(tmp$log2FoldChange.x) - 2^(tmp$log2FoldChange.y)

df_DEI <- data.frame(Celltype = names(table(tmp$Celltype.x)),
                     Number = table(tmp$Celltype.x))
df_DEI$Celltype <- factor(df_DEI$Celltype, levels = ct_levels)
ggplot(df_DEI, aes(x = Celltype, y = Number.Freq, fill = Celltype)) + 
  geom_col(position = 'dodge') +
  geom_text(aes(label = Number.Freq), position = position_dodge(width = 0.9), vjust = -0.25, size = 5) + 
  coord_cartesian(ylim = c(0, max(df_DEI$Number.Freq) * 1.1)) + 
  ylab('Number of DEIs with larger FC than gene level') + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10),
        legend.position="top",legend.justification="right") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = cols) + NoLegend()
library(DESeq2)
colData <- data.frame(indv = str_split_fixed(colnames(count_mtx_sum), "-", 2)[,2],
                      celltype = str_split_fixed(colnames(count_mtx_sum), "-", 2)[,1])
dds <- DESeqDataSetFromMatrix(countData = count_mtx_sum, colData = colData, design = ~ celltype)
# Estimate size factors
dds <- estimateSizeFactors(dds, type = "poscounts")
normalized_counts <- counts(dds, normalized=TRUE)
library_size <- colSums(count_mtx_sum)
gene_list <- c("ENSG00000105974", "ENSG00000119888", "ENSG00000149021","ENSG00000196188")

idx <- which(TALON_afterqc_orf_secondpass$annot_gene_id %in% gene_list)
NB.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_list.rds")
isoforms.list <- lapply(NB.list, function(x)x$phenotype_id)
isoform_tested <- unique(Reduce(c, isoforms.list))
isoform_tested_names <- TALON_afterqc_orf_secondpass$transcript_name_unique[which(TALON_afterqc_orf_secondpass2$annot_transcript_id %in% isoform_tested)]
isoform_list <- TALON_afterqc_orf_secondpass$transcript_name_unique[idx]
isoform_list <- intersect(isoform_list, isoform_tested_names)
count_iso_raw <- as.matrix(count_mtx_sum[isoform_list,])
count_iso_raw <- rbind(count_iso_raw, library_size)
count_iso_norm <- normalized_counts[isoform_list,]

write.table(count_iso_raw, "/data/Choi_lung/scLongreads/DEI/Count_mtx_isoforms_4genes.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(count_iso_norm, "/data/Choi_lung/scLongreads/DEI/DEseq2_norm_count_mtx_isoforms_4genes.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

write.table(DEI.list[["Alveolar_transitional_cells"]], "/data/Choi_lung/scLongreads/DEI/DEseq2_Alveolar_transitional_cells_rst.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


# heat map for suppl
DEI_top_test <- DEI.combined_sig %>% group_by(Celltype) %>% mutate(rank = rank(-log2FoldChangewhole))

DEI_top_test <- subset(DEI_top_test, rank < 6)
DEI_top_test$Celltype <- factor(DEI_top_test$Celltype, levels = gsub(" ", "_", levels(lr.isoform.sub$Celltype)))
DEI_top_test <- DEI_top_test[order(DEI_top_test$Celltype),]
DEIs <- unique(DEI_top_test$X)

norm_mtx_sum_hm <- norm_mtx_sum[DEIs,]
indv_group <- str_split_fixed(colnames(norm_mtx_sum_hm), "-", 2)[,1]

indv_group <- factor(indv_group, levels = gsub(" ", "_", levels(lr.isoform.sub$Celltype)))
idxs <- order(indv_group)
annot_col <- data.frame(row.names = colnames(norm_mtx_sum_hm), Celltype = indv_group)
annot_col$Celltype <- factor(annot_col$Celltype, levels =  gsub(" ", "_", levels(lr.isoform.sub$Celltype)))
df <- as.data.frame(table(annot_col$Celltype))

names(cellcolors) <- gsub(" ", "_", levels(lr.isoform.sub$Celltype))
ann_colors = list(
  Celltype = cellcolors)
norm_mtx_sum_hm <- norm_mtx_sum_hm[,idxs]
library(pheatmap)

p <- pheatmap(as.matrix(norm_mtx_sum_hm), cluster_rows = FALSE, scale = "row",
              cluster_cols = FALSE,show_colnames = FALSE, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                    "RdBu")))(100),
              annotation_col = annot_col, annotation_colors = ann_colors,  gaps_col = cumsum(df$Freq))
p <- pheatmap(as.matrix(norm_mtx_sum_hm), cluster_rows = FALSE, scale = "row",
              cluster_cols = FALSE,show_colnames = FALSE, color = colorRampPalette(c("#4575B4", "white", "#800000"))(100),
              annotation_col = annot_col, annotation_colors = ann_colors,  gaps_col = cumsum(df$Freq))

pdf("/data/lib14/project/scLongread/FigureSA_subgroup.pdf",height = 25, width = 20)
print(p)
dev.off()
DEIs_subgroup <- DEIs
DEI_logFC_larger_num.list <- list()
DEI_logFC_larger_num <- NULL
DEI_gene_num <- NULL
DEI_gene_num.list <- list()
for (celltype in celltypes) {
  AT1_DEI <- DEI_sig.list[[celltype]]
  AT1_DEI_unfilter <- DEG_sigpadj.list[[celltype]]
  AT1_DEI$gene <- TALON_afterqc_orf_secondpass2[AT1_DEI$X,"annot_gene_name"]
  AT1_DEI$Celltype <- celltype
  AT1_DEI_unfilter$Celltype <- celltype
  rownames(AT1_DEI_unfilter) <- AT1_DEI_unfilter$X
  tmp <- AT1_DEI %>% group_by(gene) %>% summarise(max_logFC = max(log2FoldChange))
  tmp <- tmp[!is.na(tmp$gene),]
  rownames(tmp) <- tmp$gene
  DEI_gene_num.list[[celltype]] <- tmp$gene
  DEI_gene_num <- c(DEI_gene_num, nrow(tmp))
  gl <- intersect(tmp$gene, AT1_DEI_unfilter$X)
  N <- length(which(tmp[gl,"max_logFC"] > AT1_DEI_unfilter[gl, "log2FoldChange"]))
  DEI_logFC_larger_num <- c(DEI_logFC_larger_num, N)
  DEI_logFC_larger_num.list[[celltype]] <- gl[which(tmp[gl,"max_logFC"] > AT1_DEI_unfilter[gl, "log2FoldChange"])]
  AT1_DEI -> DEI_sig.list[[celltype]]
  AT1_DEI_unfilter -> DEG_sigpadj.list[[celltype]]
}
DEI_logFC_larger_num.list <- do.call(c, DEI_logFC_larger_num.list)
DEI_gene_num.list <- do.call(c, DEI_gene_num.list)
length(unique(DEI_gene_num.list))
length(unique(DEI_logFC_larger_num.list))
df <- data.frame(DEI_gene_num = DEI_gene_num,
                 DEI_largerthan_DEG = DEI_logFC_larger_num)
df$pct <- round((df$DEI_largerthan_DEG/df$DEI_gene_num)*100,2)
AT1_DEI <- DEI_sig.list$AT1
AT1_DEI$gene <- TALON_afterqc_orf_secondpass2[AT1_DEI$X,"annot_gene_name"]
AT1_DEI_unfilter <- DEG_sigpadj.list$AT1
rownames(AT1_DEI_unfilter) <- AT1_DEI_unfilter$X
tmp <- AT1_DEI %>% group_by(gene) %>% summarise(max_logFC = max(log2FoldChange))
tmp <- tmp[!is.na(tmp$gene),]
rownames(tmp) <- tmp$gene
gl <- intersect(tmp$gene, AT1_DEI_unfilter$X)
gl[which(tmp[gl,"max_logFC"] > AT1_DEI_unfilter[gl, "log2FoldChange"])]

for (i in 1:37) {
  DEI_whole <- DEI.list[[i]]
  DEG <- DEG.list[[i]]
  DEI_whole$Celltype <- celltypes[i]
  DEG$Celltype <- celltypes[i]
  DEI.list[[i]] <- DEI_whole
  DEG.list[[i]] <- DEG
}
DEG_sigpadj <- do.call(rbind, DEG.list)
DEI_sig <- do.call(rbind, DEI.list)

DEG_sigpadj <- do.call(rbind, DEG_sigpadj.list)
DEI_sig <- do.call(rbind, DEI_sig.list)
DEI_sig$gene <- TALON_afterqc_orf_secondpass2[DEI_sig$X,"annot_gene_name"]
DEI_sig$gene_ct_id <- paste(DEI_sig$gene, DEI_sig$Celltype, sep = "-")
DEG_sigpadj$gene_ct_id <- paste(DEG_sigpadj$X, DEG_sigpadj$Celltype, sep = "-")

DEIs_sig <- left_join(DEI_sig, DEG_sigpadj, by = "gene_ct_id", suffix = c(".DEI",".DEG"))
length(unique(DEIs_sig$X.DEI[which(DEIs_sig$log2FoldChange.DEI > DEIs_sig$log2FoldChange.DEG)]))
length(unique(DEIs_sig$X.DEI)
)
DEIs_sig_sub <- subset(DEIs_sig, gene %in% HVG_3k)
length(unique(DEIs_sig_sub$gene))
DEIs_sig_sub.sig <- subset(DEIs_sig_sub, padj.DEI < 0.05 & log2FoldChange.DEI > 0)

length(unique(DEIs_sig_sub.sig$gene))
length(unique(DEIs_sig_sub.sig$X.DEI))
DEIs_sig_sub.sig.final <- DEIs_sig_sub.sig[which(DEIs_sig_sub.sig$log2FoldChange.DEI > DEIs_sig_sub.sig$log2FoldChange.DEG),]
length(unique(DEIs_sig_sub.sig.final$gene))
length(unique(DEIs_sig_sub.sig$X.DEI))
most_diff_DEI <- DEIs_sig_sub.sig.final %>% group_by(gene) %>% summarise(most_diff_isoform = X.DEI[which.max(log2FoldChange.DEI)])

saveRDS(most_diff_DEI, file = "/data/Choi_lung/scLongreads/Seurat/final/most_diff_DEI.rds")
most_diff_DEI <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/most_diff_DEI.rds")
######################
#### permutation #####
######################


files <- list.files("/data/Choi_lung/scLongreads/DEI/WholegroupPhaseI/output_permutation/UnfilteredResults/SMC/", pattern = ".csv")


celltypes <- str_split_fixed(files, "_", 4)[,4]
celltypes <- str_split_fixed(celltypes, "\\.", 2)[,1]
celltypes <- gsub(" ", "_", celltypes)

DEI_per.list <- list()
Sig_DEI <- NULL
for (i in 1:1000) {
  DEI_rst <- read.csv(paste0("/data/Choi_lung/scLongreads/DEI/WholegroupPhaseI/output_permutation/UnfilteredResults/SMC/",files[i]))
  DEI_per.list[[i]] <- DEI_rst
  Sig_DEI <- c(length(which(DEI_rst$padj < 0.05)), Sig_DEI)
}
# Sig_Num_list <- list()
Sig_Num_list[["Mast_cells"]] <- Sig_DEI
Sig_Num_list[["SMC"]] <- Sig_DEI
Sig_Num_list[["Myofibroblast"]] <- Sig_DEI
hist(Sig_DEI)
tmp.list <- lapply(DEI_per.list, function(x){
  x <- as.data.frame(x[,c(7)])
  return(x)
})
# names(cellcolors) <- gsub(" ", "_",levels(lr$Celltype))
permutation_100times <- do.call(cbind, tmp.list)
permutation_100times_sig <- rowSums((permutation_100times < 0.05), na.rm = TRUE)

hist(permutation_100times_sig, breaks = 100)
tmp <- data.frame(Sig_times = permutation_100times_sig,
                  Celltype = "Myofibroblast")

# df_per <- data.frame(Sig_times = permutation_100times_sig,
#                      Celltype = "AT2")
df_per <- rbind(df_per, tmp)
library(ggplot2)
library(dplyr)
library(patchwork)
ggplot(df_per, aes(x = Sig_times, fill = Celltype)) +
  geom_histogram(aes(y = ..density..),position="identity",binwidth = 1) + 
  geom_vline(xintercept = 50, linetype = 2) + 
  scale_fill_manual(values = cellcolors) + facet_wrap(~Celltype, nrow = 1) + 
  ylab("Fraction of tested isoforms") + xlab("# of permutated data where a isoform is falsely claimed as DEI") +
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")
ggsave("/data/lib14/project/scLongread/FigSMethodB.pdf", width = 20,height = 7)
table(df_per$Sig_times)
df_per$FDR_control <- TRUE
df_per$FDR_control[which(df_per$Sig_times < 50)] <- FALSE
table(df_per$Celltype, df_per$FDR_control)

rst <- DEI.list$Myofibroblast
rst <- DEI.list$Mast_cells
rst <- DEI.list$SMC
plot((rst$log2FoldChange), -log10(rst$padj))
rst$Sig_times <- df_per$Sig_times[which(df_per$Celltype == "SMC")]
rst$FP <- df_per$FDR_control[which(df_per$Celltype == "SMC")]
rst$gene_type <- "ns"
rst$gene_type[which(rst$padj < 0.05 & rst$log2FoldChange > 1)] <- "up"
rst$gene_type[which(rst$padj < 0.05 & rst$log2FoldChange < -1)] <- "down"
table(rst$gene_type, rst$FP)
rst$gene_type[which(rst$FP & rst$padj < 0.05 & abs(rst$log2FoldChange) > 1)] <- "FP"
# Add colour, size and alpha (transparency) to volcano plot --------------------
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey", "FP"="red") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1, "FP"=2) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5,"FP"=1)

rst %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  geom_point(shape = 21, colour ="black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-10, 10, 2)),       
                     limits = c(-10, 10)) +
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")+
  ggtitle("SMC")
ggsave("/data/lib14/project/scLongread/FigSMethodC3_FP_1kperm.pdf", width = 8,height = 7)


df <- data.frame(DEI = sapply(DEI_sig.list[gsub(" ","_",levels(lr$Celltype))], nrow), 
                 Total = sapply(DEI.list[gsub(" ","_",levels(lr$Celltype))], nrow),
                 cellnum = table(lr$Celltype))
df$sig_pct <- df$DEI/df$Total
df$sig_pct1 <- round((df$DEI/df$Total)*100, 2)
cor(df$sig_pct1, df$cellnum.Freq)
cor(df$DEI, df$cellnum.Freq)
df$cellnum.Var1
ggplot(df, aes(x = cellnum.Var1, y = sig_pct1, fill = cellnum.Var1)) + 
  geom_bar(stat="identity",position = "identity") + scale_fill_manual(values = cellcolors) +
  coord_cartesian(ylim = c(0, max(df$sig_pct1) * 1.1)) +
  ylab('Percentage of significant') + xlab("Cell type")+
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        # legend.position="none",
        legend.justification="right") +
  scale_y_continuous(expand = c(0,0))

# Supplementary figure 5a
DEI_DEG_comp.list <- list()
celltypes <- gsub(" ", "_", celltypes)
for (celltype in celltypes) {
  AT1_DEI <- DEI_sig.list[[celltype]]
  AT1_DEI$gene <- TALON_afterqc_orf_secondpass[AT1_DEI$X, "annot_gene_name"]
  # AT1_DEI_sig_top
  AT1_DEI_sig_top <- AT1_DEI %>% group_by(gene) %>% mutate(log2FC_DEI_max = max(log2FoldChange))
  AT1_DEI_sig_top <- AT1_DEI_sig_top[which(AT1_DEI_sig_top$log2FoldChange == AT1_DEI_sig_top$log2FC_DEI_max),]
  AT1_DEG <- DEG.list[[celltype]]
  AT1_DEG <- subset(AT1_DEG, log2FoldChange > 0)
  colnames(AT1_DEG)[1] <- "gene"
  AT1_DEI_sig_top_sum <- left_join(AT1_DEI_sig_top, AT1_DEG, by = "gene")
  AT1_DEI_sig_top_sum <- AT1_DEI_sig_top_sum[!is.na(AT1_DEI_sig_top_sum$padj.y),]
  AT1_DEI_sig_top_sum <- subset(AT1_DEI_sig_top_sum, log2FoldChange.y > 0 & log2FoldChange.x > 0)
  plot(AT1_DEI_sig_top_sum$log2FoldChange.x, AT1_DEI_sig_top_sum$log2FoldChange.y)
  AT1_DEI_sig_top_sum$DEI <- FALSE
  AT1_DEI_sig_top_sum$DEI <- (AT1_DEI_sig_top_sum$log2FoldChange.x > 1 & AT1_DEI_sig_top_sum$padj.x < 0.05)
  AT1_DEI_sig_top_sum <- subset(AT1_DEI_sig_top_sum, padj.x < 0.05 & padj.y < 0.05)
  AT1_DEI_sig_top_sum$Celltype <- celltype
  DEI_DEG_comp.list[[celltype]] <- AT1_DEI_sig_top_sum
}
# AT1_DEI <- DEI_sig.list$AT1
# AT1_DEI$gene <- TALON_afterqc_orf_secondpass[AT1_DEI$X, "annot_gene_name"]
# # AT1_DEI_sig_top
# AT1_DEI_sig_top <- AT1_DEI %>% group_by(gene) %>% mutate(log2FC_DEI_max = max(log2FoldChange))
# AT1_DEI_sig_top <- AT1_DEI_sig_top[which(AT1_DEI_sig_top$log2FoldChange == AT1_DEI_sig_top$log2FC_DEI_max),]
# AT1_DEG <- DEG.list
# AT1_DEG <- subset(AT1_DEG, log2FoldChange > 0)
# colnames(AT1_DEG)[1] <- "gene"
# AT1_DEI_sig_top_sum <- left_join(AT1_DEI_sig_top, AT1_DEG, by = "gene")
# AT1_DEI_sig_top_sum <- AT1_DEI_sig_top_sum[!is.na(AT1_DEI_sig_top_sum$padj.y),]
# AT1_DEI_sig_top_sum <- subset(AT1_DEI_sig_top_sum, log2FoldChange.y > 0 & log2FoldChange.x > 0)
# plot(AT1_DEI_sig_top_sum$log2FoldChange.x, AT1_DEI_sig_top_sum$log2FoldChange.y)
# AT1_DEI_sig_top_sum$DEI <- FALSE
# AT1_DEI_sig_top_sum$DEI <- (AT1_DEI_sig_top_sum$log2FoldChange.x > 1 & AT1_DEI_sig_top_sum$padj.x < 0.05)
# AT1_DEI_sig_top_sum <- subset(AT1_DEI_sig_top_sum, padj.x < 0.05 & padj.y < 0.05)
# cellcolors
# AT1_DEI_sig_top_sum$delta <- AT1_DEI_sig_top_sum$log2FoldChange.x - AT1_DEI_sig_top_sum$log2FoldChange.y
# AT1_DEI_sig_top_sum$delta_rank <- rank(-AT1_DEI_sig_top_sum$delta)
# AT1_DEI_sig_top_sum$label <- NA
# AT1_DEI_sig_top_sum$label[which(AT1_DEI_sig_top_sum$delta_rank <= 21)] <- AT1_DEI_sig_top_sum$gene[which(AT1_DEI_sig_top_sum$delta_rank <= 21)]
DEI_DEG_comp <- do.call(rbind, DEI_DEG_comp.list)
DEI_DEG_comp$Celltype <- factor(DEI_DEG_comp$Celltype, levels = gsub(" ","_",levels(metadata$Celltype)))
ggplot(DEI_DEG_comp, aes(x = log2FoldChange.y, y = log2FoldChange.x)) + 
  geom_point(aes(color = Celltype), size = .5) + scale_color_manual(values = c(cellcolors))+
  xlab('log2 fold change at gene level') + ylab("log2 fold change of top DEI within each gene")+
  geom_abline(intercept = 0, slope = 1, color = "black")+ facet_wrap(~Celltype, ncol = 8)+
  coord_cartesian(ylim = c(0, max(DEI_DEG_comp$log2FoldChange.x))) + theme_bw() +
  theme(strip.background = element_rect(fill='transparent'),legend.position="none")
ggsave("figures/Fig.S5a.pdf", width = 16, height = 10)
library(ggrepel)
ggplot(AT1_DEI_sig_top_sum[AT1_DEI_sig_top_sum$DEI,], aes(x = log2FoldChange.y, y = log2FoldChange.x)) + 
  geom_point(aes(color = DEI)) + scale_color_manual(values = c(cellcolors[1]))+
  geom_text_repel(aes(label = label)) + 
  xlab('log2 fold change at gene level') + ylab("log2 fold change of top DEI within each gene")+
  geom_abline(intercept = 0, slope = 1, color = "black")+
  coord_cartesian(ylim = c(0, max(AT1_DEI_sig_top_sum$log2FoldChange.x))) +
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12,face = "plain"),
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")
ggsave("DEI_DEG_compare.pdf", width = 7, height = 7)

