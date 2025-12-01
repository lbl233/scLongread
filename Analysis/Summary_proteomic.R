library("seqinr")
library(stringi)
library(stringr)
library(dplyr)
fastafile<- read.fasta(file = "Need_clean/NCI_lung_eIsoforms.fasta", 
                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
cds_fastafile<- read.fasta(file = "/data/Choi_lung/scLongreads/TALON_workspace/test10s/transcripts.fasta.transdecoder.cds", 
                       seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
NB.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")
eIsoforms <- unique(do.call(rbind, isoQTL_sig_list)$phenotype_id)
pro_tested <- unique(str_split_fixed(names(fastafile), "\\.", 2)[,1])
sum(grepl("TALONT", pro_tested))
sum(grepl("TALONT", subsetlist))
18697/length(subsetlist)
NCI134 <- read.table("Need_clean/NCI134.txt", header = TRUE, sep = "\t", comment.char = "")
NCI134 <- subset(NCI134, Protein.FDR.Confidence..Combined == "High" )
NCI135 <- read.table("Need_clean/NCI135.txt", header = TRUE, sep = "\t", comment.char = "")
NCI135 <- subset(NCI135, Protein.FDR.Confidence..Combined == "High" )
NCI138 <- read.table("Need_clean/NCI138.txt", header = TRUE, sep = "\t", comment.char = "")
NCI138 <- subset(NCI138, Protein.FDR.Confidence..Combined == "High" )
NCI139 <- read.table("Need_clean/NCI139.txt", header = TRUE, sep = "\t", comment.char = "")
NCI139 <- subset(NCI139, Protein.FDR.Confidence..Combined == "High" )
all_val_pro <- unique(Reduce(union, list(NCI134$Accession, NCI135$Accession, NCI138$Accession, NCI139$Accession)))

pro_detection <- data.frame(Accession = all_val_pro,
                            NCI134 = (all_val_pro %in% NCI134$Accession),
                            NCI135 = (all_val_pro %in% NCI135$Accession),
                            NCI138 = (all_val_pro %in% NCI138$Accession),
                            NCI139 = (all_val_pro %in% NCI139$Accession))
idx <- which(rowSums(pro_detection[,c(2:5)]) >= 2)
pro_detection$MoreThanOne <- rowSums(pro_detection[,c(2:5)]) >= 2
pro_detection$MoreThanTwo <- rowSums(pro_detection[,c(2:5)]) >= 3
pro_detection$All <- rowSums(pro_detection[,c(2:5)]) == 4
pro_detection$transcript_id <- str_split_fixed(pro_detection$Accession, pattern = "\\.", 2)[,1]
pro_detection <- pro_detection %>% group_by(transcript_id) %>% mutate(transcript_name = Search_transcript_name2(transcript_id))
sum(pro_detection$MoreThanOne)/length(subsetlist)

pro_detection <- left_join(pro_detection, NCI134[,c("Accession","X..Unique.Peptides")], by = "Accession")
colnames(pro_detection)[11] <- "Num_unique_peptides_NCI134"
pro_detection <- left_join(pro_detection, NCI135[,c("Accession","X..Unique.Peptides")], by = "Accession")
colnames(pro_detection)[12] <- "Num_unique_peptides_NCI135"
pro_detection <- left_join(pro_detection, NCI138[,c("Accession","X..Unique.Peptides")], by = "Accession")
colnames(pro_detection)[13] <- "Num_unique_peptides_NCI138"
pro_detection <- left_join(pro_detection, NCI139[,c("Accession","X..Unique.Peptides")], by = "Accession")
colnames(pro_detection)[14] <- "Num_unique_peptides_NCI139"
pro_detection <- pro_detection %>% group_by(Accession) %>% mutate(Max_Num_unique_peptides = max(c(Num_unique_peptides_NCI134,
                                                                                 Num_unique_peptides_NCI135,
                                                                                 Num_unique_peptides_NCI138,
                                                                                 Num_unique_peptides_NCI139), na.rm = TRUE))
pro_detection$Novelty <- "known"
pro_detection$Novelty[grepl("TALONT", pro_detection$Accession)] <- "Novel"
table(pro_detection$Novelty)
pro_detection_only_one_unique_peptide <- subset(pro_detection, Max_Num_unique_peptides < 2)
pro_detection_hc <- subset(pro_detection, Max_Num_unique_peptides >= 2 | MoreThanOne == TRUE)
eIsoform_pro_val <- unique(eIsoforms[which(eIsoforms %in% pro_detection_hc$transcript_id)])
sum(grepl("TALONT", eIsoforms))
sum(grepl("TALONT", eIsoform_pro_val))
table(pro_detection_only_one_unique_peptide$Novelty)
cds_fastafile <- cds_fastafile[pro_detection$Accession]
# Novel isoform
sum(grepl("TALONT", pro_detection$Accession))/18697 # 9.7%
sum(grepl("TALONT", pro_detection$Accession[pro_detection$MoreThanOne]))/18697 # 8.0%
sum(grepl("TALONT", pro_detection$Accession[pro_detection$All]))/18697 # 5.6%

# Known isoform
sum(!grepl("TALONT", pro_detection$Accession))/(length(subsetlist) - 18697)  # 24.7%
sum(!grepl("TALONT", pro_detection$Accession[pro_detection$MoreThanOne]))/(length(subsetlist) - 18697)  # 21.5%
sum(!grepl("TALONT", pro_detection$Accession[pro_detection$All]))/(length(subsetlist) - 18697)  # 16.6%

# Novel isoform
sum(grepl("TALONT", pro_detection_hc$Accession))/18697 # 8.7%
sum(grepl("TALONT", pro_detection_hc$Accession[pro_detection_hc$MoreThanOne]))/18697 # 8.0%
sum(grepl("TALONT", pro_detection_hc$Accession[pro_detection_hc$All]))/18697 # 5.6%

# Known isoform
sum(!grepl("TALONT", pro_detection_hc$Accession))/(length(subsetlist) - 18697)  # 22.8%
sum(!grepl("TALONT", pro_detection_hc$Accession[pro_detection_hc$MoreThanOne]))/(length(subsetlist) - 18697)  # 21.5%
sum(!grepl("TALONT", pro_detection_hc$Accession[pro_detection_hc$All]))/(length(subsetlist) - 18697)  # 16.6%



df_overall <- data.frame(pct = c(8.7, 8.0, 5.6, 22.8, 21.5, 16.6), 
                         Novelty = rep(c("Novel", "Known"), each = 3),
                         detection = rep(c("Any", "More than one", "All"), 2))
df_overall$detection <- factor(df_overall$detection, levels = c("Any", "More than one", "All"))
ggplot(subset(df_overall, detection != "More than one"), aes(x = detection, y = pct, fill = Novelty)) +
  geom_col(position = "dodge") +
  geom_text(
    aes(label = pct),
    colour = "black", size = 4.5,
    vjust = 1.5, position = position_dodge(.9)
  )+ylab("Percentage of validated isoforms")+ ggtitle("Validation of 42,824 isoforms")+xlab("")+
  theme(axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        
        axis.text.x = element_text(color = "black", size = 12),
        strip.background = element_blank(), strip.text = element_text(size = 14))

ggsave("/data/lib14/project/scLongread/Fig2G.pdf", width = 7,height = 7)

lr.isoform.sub <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final_w_isoform_expr.RDS")
lr <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final_v2_wRNAexpr.RDS")
lr.isoform.sub$Sample_NCI <- lr$Sample_NCI
count_mtx <- lr.isoform.sub@assays$Isoform@counts
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass2$transcript_name_unique,]
group <- as.character(lr.isoform.sub$Sample_NCI) %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx)
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
sum(count_mtx_sum)
count_mtx_sum <- Matrix::t(count_mtx_sum)
dim(count_mtx_sum)
indv_expr_tf <- (count_mtx_sum > 0)

colnames(indv_expr_tf)
head(rownames(indv_expr_tf))
iso_detection <- as.data.frame(indv_expr_tf[pro_detection$transcript_name,c("NCI_134", "NCI_135", "NCI_138", "NCI_139")])
iso_detection_all <- as.data.frame(indv_expr_tf[pro_detection$transcript_name,])

pro_detection$NCI134_trans <- TRUE
pro_detection$NCI135_trans <- TRUE
pro_detection$NCI138_trans <- TRUE
pro_detection$NCI139_trans <- TRUE
pro_detection$NCI134_trans[which(!iso_detection_all$NCI_134)] <- FALSE
pro_detection$NCI135_trans[which(!iso_detection_all$NCI_135)] <- FALSE
pro_detection$NCI138_trans[which(!iso_detection_all$NCI_138)] <- FALSE
pro_detection$NCI139_trans[which(!iso_detection_all$NCI_139)] <- FALSE
NCI134_tranNO_proYes <- unique(pro_detection$transcript_id[which(pro_detection$NCI134 & !(pro_detection$NCI134_trans))])
NCI135_tranNO_proYes <- unique(pro_detection$transcript_id[which(pro_detection$NCI135 & !(pro_detection$NCI135_trans))])
NCI138_tranNO_proYes <- unique(pro_detection$transcript_id[which(pro_detection$NCI138 & !(pro_detection$NCI138_trans))])
NCI139_tranNO_proYes <- unique(pro_detection$transcript_id[which(pro_detection$NCI139 & !(pro_detection$NCI139_trans))])

sum(grepl("TALONT",NCI134_tranNO_proYes))
sum(grepl("TALONT",NCI135_tranNO_proYes))
sum(grepl("TALONT",NCI138_tranNO_proYes))
sum(grepl("TALONT",NCI139_tranNO_proYes))


count_mtx <- lr.isoform.sub@assays$Isoform@counts
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass2$transcript_name_unique,]
lr.isoform.sub$orig.ident <- as.character(lr.isoform.sub$orig.ident)
group <- lr.isoform.sub$orig.ident %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx)
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
sum(count_mtx_sum)
count_mtx_sum <- Matrix::t(count_mtx_sum)
dim(count_mtx_sum)
iso_detection_batch <- as.data.frame((count_mtx_sum > 0)[pro_detection$transcript_name,])
pro_detection <- cbind(pro_detection, iso_detection_batch[,c("20_NCI_130_135","21_NCI_136_140")])

Batch20_tranNO_proYes <- unique(pro_detection$transcript_id[which(pro_detection$NCI134 & !(pro_detection$`20_NCI_130_135`))])
Batch20_1_tranNO_proYes <- unique(pro_detection$transcript_id[which(pro_detection$NCI135 & !(pro_detection$`20_NCI_130_135`))])
Batch21_tranNO_proYes <- unique(pro_detection$transcript_id[which(pro_detection$NCI138 & !(pro_detection$`21_NCI_136_140`))])
Batch21_1_tranNO_proYes <- unique(pro_detection$transcript_id[which(pro_detection$NCI139 & !(pro_detection$`21_NCI_136_140`))])

tranNO_proYes_df <- data.frame(Number = c(length(NCI134_tranNO_proYes), length(NCI135_tranNO_proYes), length(NCI138_tranNO_proYes), length(NCI139_tranNO_proYes),
                                          length(Batch20_tranNO_proYes), length(Batch20_1_tranNO_proYes), length(Batch21_tranNO_proYes), length(Batch21_1_tranNO_proYes)),
                               Sample = rep(c("NCI134", "NCI135", "NCI138", "NCI139"), 2), 
                               Batch = rep(c("Batch20", "Batch20", "Batch21", "Batch21"), 2),
                               Level = rep(c("Sample", "Batch"), each = 4))
ggplot(tranNO_proYes_df, aes(x = Sample, y = Number, fill = Level)) +
  geom_col(position = "dodge") +
  geom_text(
    aes(label = Number),
    colour = "black", size = 4.5,
    vjust = 1.5, position = position_dodge(.9)
  )+ylab("Number of isoforms")+ ggtitle("Isoforms detected at protein level but not transcript level")+xlab("") + 
  theme(axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        
        axis.text.x = element_text(color = "black", size = 12),
        strip.background = element_blank(), strip.text = element_text(size = 14))
ggsave("/data/lib14/project/scLongread/Fig2G-2.pdf", width = 8,height = 7)
sum(pro_detection_hc$NCI134)
sum(pro_detection_hc$NCI135)
sum(pro_detection_hc$NCI138)
sum(pro_detection_hc$NCI139)
"MLLFLLSALVLLTQPLGYLEAEMKTYSHRTMPSACTLVMCSSVESGLPGRDGRDGREGPRGEKGDPGLPGAAGQAGMPGQAGPVGPKGDNGSVGEPGPKGDTGPSGPPGPPGVPGPAGREGPLGKQGNIGPQGKPGPKGEAGPKGEVGAPGMQGSAGARGLAGPKGERGVPGERGVPGNTGAAGSAGAMGPQGSPGARGPPGLKGDKGIPGDKGAKGESGLPDVASLRQQVEALQGQVQHLQAAFSQYKKVELFPNGQSVGEKIFKTAGFVKPFTEAQLLCTQAGGQLASPRSAAENAALQQLVVAKNEAAFLSMTDSKTEGKFTYPTGESLVYSNWAPGEPNDDGGSEDCVEIFTNGKWNDRACGEKRLVVCEF"
"MILNKALLLGALALTTVMSPCGGEDIVADHVASCGVNLYQFYGPSGQYTHEFDGDEQFYVDLERKETAWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRYNSTAATNEVPEVTVFSKFPVTLGQPNTLICLVDNIFPPVVNITWLSNGHSVTEGVSETSFLSKSDHSFFKISYLTFLPSADEIYDCKVEHWGLDEPLLKHWEPEIPAPMSELTETLVCALGLSVGLMGIVVGTVFIIQGLRSVGASRHQGLL"
"MILNKALLLGALALTTVMSPCGGEDIVADHVASCGVNLYQFYGPSGQYTHEFDGDEQFYVDLERKETAWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRYNSTAATNEVPEVTVFSKSPVTLGQPNTLICLVDNIFPPVVNITWLSNGQSVTEGVSETSFLSKSDHSFFKISYLTFLPSADEIYDCKVEHWGLDQPLLKHWEPEIPAPMSELTETVVCALGLSVGLMGIVVGTVFIIQGLRSVGASRHQGPL"
tmp <- unique(DEIs_sig$X.DEI[which(DEIs_sig$log2FoldChange.DEI > DEIs_sig$log2FoldChange.DEG)])[which(unique(DEIs_sig$X.DEI[which(DEIs_sig$log2FoldChange.DEI > DEIs_sig$log2FoldChange.DEG)]) %in% pro_detection$transcript_name)]
df <- data.frame(transcrip_name = tmp,
                 gene_name = TALON_afterqc_orf_secondpass2[tmp,"annot_gene_name"])
View(df[which(df$gene_name %in% c(Epi_sig_list, Imm_sig_list, Endo_sig_list,Stroma_sig_list)),])
"MSIRVTQKSYKVSTSGPRAFSSRSYTSGPGSRISSSSFSRVGSSNFRGGLGGGYGGASGMGGITAVTVNQSLLSPLVLEVDPNIQAVRTQEKEQIKTLNNKFASFIDKVRFLEQQNKMLETKWSLLQQQKTARSNMDNMFESYINNLRRQLETLGQEKLKLEAELGNMQGLVEDFKNKYEDEINKRTEMENEFVLIKKDVDEAYMNKVELESRLEGLTDEINFLRQLYEEEIRELQSQISDTSVVLSMDNSRSLDMDSIIAEVKAQYEDIANRSRAEAESMYQIKYEELQSLAGKHGDDLRRTKTEISEMNRNISRLQAEIEGLKGQRASLEAAIADAEQRGELAIKDANAKLSELEAALQRAKQDMARQLREYQELMNVKLALDIEIATYRKLLEGEESRLESGMQNMSIHTKTTSGYAGGLSSAYGGLTSPGLSYSLGSSFGSGAGSSSFSRTSSSRAVVVKKIETRDGKLVSESSDVLPK"


# Glinos et al 2022
GTEx_pro_val <- read.table("/data/lib14/R/Rscripts/Glinos_TS3.txt", sep = "\t", header = TRUE)
GTEx_pro_val_sub <- subset(GTEx_pro_val, Mass_Spec == "Yes")
gff_combined_sub <- subset(gff_combined, class_code == "=")
length(which(gff_combined_sub$cmp_ref %in% GTEx_pro_val_sub$Transcript))
validated_overlaps <- unique(GTEx_pro_val_sub$Transcript[which(GTEx_pro_val_sub$Transcript %in% gff_combined_sub$cmp_ref)])
View(GTEx_pro_val_sub[which(GTEx_pro_val_sub$Transcript %in% validated_overlaps),])
length(validated_overlaps) # 756
gff_combined_sub_proval <- subset(gff_combined_sub, cmp_ref %in% GTEx_pro_val_sub$Transcript)
length(which(gff_combined_sub_proval$transcript_id %in% pro_detection$transcript_id))

novel_isoform_val_bytwo <- unique(pro_detection_hc$transcript_id)[which(unique(pro_detection_hc$transcript_id) %in% gff_combined_sub_proval$transcript_id)]

lungfun_coloc_proval <- c("ENST00000344691","TALONT001781255","TALONT003093589","ENST00000296930","ENST00000616870","ENST00000338631", "ENST00000362034","ENST00000367063","TALONT001769911")
sapply(lungfun_coloc_proval, function(x) Search_transcript_name2(x))
lungcancer_TWAS_proval <- unique(TWAS_tb$ID[which(TWAS_tb$ID %in% pro_detection_hc$transcript_id)])
sapply(lungcancer_TWAS_proval, function(x) Search_transcript_name2(x))
