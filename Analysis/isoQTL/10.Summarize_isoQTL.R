# Summarize the isoQTL results in epithelial populations

# Bolun Li
# Jul 1 2024


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
source("/data/lib14/R/Rscripts/utilities.R")
# compare expr and ratio
file_paths <- list.files("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/output_iso/", 
                         pattern = "cis_qtl.txt.gz")
file_paths <- file_paths[c(14:25)]
file_paths <- file_paths[c(14:37,50:73)]
file_paths <- file_paths[c(26:37,50:73)] # cutoff
file_paths <- file_paths[c(2:25,86:97)]
rst.sig.iso.list <- list()
rst.iso.list <- list()
for (i in c(1:length(file_paths))) {
  file_path <- paste0("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/output_iso/", file_paths[i])
  covs_used <- str_split_fixed(str_split_fixed(file_paths[i], "_", n = 5)[,5], "\\.", n = 2)[,1]
  rst <- read.table(gzfile(file_path),sep="\t", header = TRUE)
  rst.iso.list[[covs_used]] <- rst
  rst <- subset(rst, qval < 0.05)
  rst.sig.iso.list[[covs_used]] <- rst
}

file_paths <- list.files("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/output_iso20/", 
                         pattern = "cis_qtl.txt.gz")
file_paths <- file_paths[1:12]
file_paths <- file_paths[17:19]
# rst.sig.iso.list <- list()
# rst.iso.list <- list()
for (i in c(1:length(file_paths))) {
  file_path <- paste0("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/output_iso20/", file_paths[i])
  covs_used <- str_split_fixed(str_split_fixed(file_paths[i], "_", n = 5)[,5], "\\.", n = 2)[,1]
  rst <- read.table(gzfile(file_path),sep="\t", header = TRUE)
  rst.iso.list[[covs_used]] <- rst
  rst <- subset(rst, qval < 0.05)
  rst.sig.iso.list[[covs_used]] <- rst
}

# isoform expressed by > 30% indvs
iso30_ratio <- rst.sig.iso.list$iso30_45PF
iso30_expr <- rst.sig.iso.list$iso30_expr_35PF

length(intersect(iso30_ratio$phenotype_id, iso30_expr$phenotype_id))

# Q-Q plot
exp = -(log10(c(1:length(rst.iso.list[["iso30_expr_35PF"]]$pval_nominal))/(length(rst.iso.list[["iso30_expr_35PF"]]$pval_nominal)+ 1)))
obs = -log10(sort(rst.iso.list[["iso30_expr_35PF"]]$pval_nominal))
plot(exp, obs, ylab="Observed(−logP)",xlab="Expected(−logP)",ylim=c(0,50),xlim=c(0,4),main = "Nominal P value") 
lines(c(0,4),c(0,4),col=1,lwd=2)

df_isoQTL <- data.frame(covs = names(rst.sig.iso.list),
                        Sig = sapply(rst.sig.iso.list, nrow))
df_isoQTL$cutoff <- str_split_fixed(df_isoQTL$covs, "_", 2)[,1]
df_isoQTL$cutoff <- gsub("iso", "", df_isoQTL$cutoff)
df_isoQTL$cutoff <- paste0(df_isoQTL$cutoff, "%")
df_isoQTL$cutoff  <- factor(df_isoQTL$cutoff, levels = paste0(c(20,30,40,50), "%"))
df_isoQTL$covs_pf <- str_split_fixed(df_isoQTL$covs, "_", 2)[,2]
df_isoQTL$value_type <- "ratio"
df_isoQTL$value_type[grep("expr", df_isoQTL$covs_pf )] <- "expr"
df_isoQTL$covs_pf[grep("expr", df_isoQTL$covs_pf )] <- str_split_fixed(df_isoQTL$covs_pf[grep("expr", df_isoQTL$covs_pf )], "_", 2)[,2]

df_isoQTL$covs_pf <- factor(df_isoQTL$covs_pf, levels = paste0(seq(5,60,5), "PF"))
ggplot(df_isoQTL[c(8,16,29,40),], aes(x=cutoff, y=Sig)) +
  geom_bar(stat="identity", position=position_dodge(),fill = "#7CAE00") + xlab("Cutoff")+ylab("Number of significant isoforms") +
  geom_text(aes(label = Sig), , vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  theme_classic()
ggplot(df_isoQTL, aes(x=covs_pf, y=Sig, fill = value_type)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("covs")+ylab("Number of significant isoforms") +
  geom_text(aes(label = Sig), , vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  theme_classic()+facet_wrap(~cutoff)
ggplot(df_isoQTL, aes(x=covs_pf, y=Sig, fill = Celltype)) +
  geom_bar(stat="identity") + xlab("covs")+ylab("Number of significant isoforms") +
  geom_text(aes(label = Sig), , vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  theme_classic()+facet_wrap(~cutoff)
df_isoQTL$Celltype <- "AT2"

# compare across cell type
celltypes <- read.table(file = "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Celltypes.txt", sep = "\t")
celltypes <- read.table(file = "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/celltypes.txt", sep = "\t")
celltypes <- celltypes$V1
rst.sig.iso.ct.list <- list()
rst.iso.ct.list <- list()
for (celltype in celltypes) {
  celltype_for_path <- gsub(" ", "_", celltype)
  output_path <- paste0("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/", celltype_for_path, "/output/")
  file_paths <- list.files(output_path, 
                           pattern = "cis_qtl.txt.gz")
  file_paths <- file_paths[grepl("Standardized_quantiled_expr", file_paths)]
  for (i in c(1:length(file_paths))) {
    file_path <- paste0(output_path, file_paths[i])
    covs_used <- str_split_fixed(str_split_fixed(file_paths[i], "_", n = 5)[,5], "\\.", n = 2)[,1]
    covs_used <- paste(celltype_for_path, covs_used, sep = "-")
    rst <- read.table(gzfile(file_path),sep="\t", header = TRUE)
    rst.iso.ct.list[[covs_used]] <- rst
    rst <- subset(rst, qval < 0.05)
    rst.sig.iso.ct.list[[covs_used]] <- rst
  }
}

############
## permutation "significant" hits
permutation_sig_list <- rst.sig.iso.ct.list[c(3,19,29,41,52,71,74,94)]
names(permutation_sig_list) <- celltypes
saveRDS(permutation_sig_list, file = "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/permutation_sig_list.rds")
####

df_isoQTL_ct <- data.frame(covs = names(rst.sig.iso.ct.list),
                           Sig = sapply(rst.sig.iso.ct.list, nrow))
df_isoQTL_ct$Celltype <- str_split_fixed(df_isoQTL_ct$covs, "-", 2)[,1]
df_isoQTL_ct$covs_pf <- str_split_fixed(df_isoQTL_ct$covs, "-", 2)[,2]
df_isoQTL_ct$cutoff <- str_split_fixed(df_isoQTL_ct$covs_pf, "_", 2)[,1]
df_isoQTL_ct$covs_pf <- str_split_fixed(df_isoQTL_ct$covs_pf, "_", 2)[,2]
df_isoQTL_ct$value_type <- "ratio"
idxs <- which(df_isoQTL_ct$cutoff == "batch")
df_isoQTL_ct$value_type[idxs] <- "expr"
df_isoQTL_ct$cutoff[idxs] <- str_split_fixed(df_isoQTL_ct$covs_pf[idxs], "_", 2)[,1]
df_isoQTL_ct$covs_pf[idxs] <- str_split_fixed(df_isoQTL_ct$covs_pf[idxs], "_", 2)[,2]

df_isoQTL_ct$covs_pf <- factor(df_isoQTL_ct$covs_pf, levels = paste0(seq(5,60,5), "PF"))
df_isoQTL_ct <- df_isoQTL_ct[,colnames(df_isoQTL)]
df_isoQTL_ct_ept <- rbind(subset(df_isoQTL, (cutoff == "iso20") ), df_isoQTL_ct)
df_isoQTL_ct_ept$value_type <- as.character(df_isoQTL_ct_ept$value_type)
df_isoQTL_ct_ept$value_type[which(df_isoQTL_ct_ept$value_type == "expr")] <- "expression-sum"
df_isoQTL_ct_ept$value_type <- factor(df_isoQTL_ct_ept$value_type, levels = c("ratio", "expression-sum"))

ggplot(df_isoQTL_ct_ept, aes(x=covs_pf, y=Sig, fill = value_type)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("# of PEER factors")+ylab("Number of significant isoforms") +
  geom_text(aes(label = Sig), , vjust=-.3, color="black",
            position = position_dodge(0.9), size=3.5)+theme_bw()+ 
  theme(strip.text.x.top = element_text(size = 14, face = "bold", colour = "black"),
        strip.background = element_blank(),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 10, face = "bold",colour = "black"))+
  facet_wrap(~Celltype, nrow = 2)

AT2 <- subset(df_isoQTL_ct_ept, Celltype == "AT2")
AT2 <- subset(AT2, value_type == "expr")
AT2_base <- data.frame(Sig = c(162,165,164),
                       covs_pf = c("No cov", "Age", "Age+batch"))
AT2 <- rbind(AT2_base, AT2[,c(2,4)])
AT2$covs_pf <- factor(AT2$covs_pf, levels = c("No cov", "Age", "Age+batch", paste0(seq(5,60,5), "PF")))
ggplot(AT2, aes(x=covs_pf, y=Sig)) +
  geom_bar(stat="identity", position=position_dodge(),fill = "#00BFC4") + xlab("# of PEER factors")+ylab("Number of significant isoforms") +
  geom_text(aes(label = Sig), , vjust=-.3, color="black",
            position = position_dodge(0.9), size=3.5)+theme_bw()+ 
  theme(strip.text.x.top = element_text(size = 14, face = "bold", colour = "black"),
        strip.background = element_blank(),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 10, face = "bold",colour = "black"))

ATC <- rst.sig.iso.ct.list$`Alveolar_transitional_cells-batch_iso20_20PF`
STC <- rst.sig.iso.ct.list$`Secretory_transitional_cells-batch_iso20_5PF`

intersect(ATC$phenotype_id, STC$phenotype_id)
setdiff(ATC$phenotype_id, STC$phenotype_id)
setdiff(STC$phenotype_id, ATC$phenotype_id)

Search_transcript_name(intersect(ATC$phenotype_id, STC$phenotype_id))
rst.list <- c((rst.iso.ct.list),(rst.iso.list))
shared_isos <- Reduce(intersect, lapply(rst.list, function(x)x$phenotype_id))

library(ACAT)
library(qvalue)
pos_sub <- subset(pos, X %in% rst.iso.list[["5PF"]]$phenotype_id)
phenotype_groups <- pos_sub[,c(1,5)]
colnames(phenotype_groups) <- c("phenotype_id", "gene_id")
QTL_rst_mtx <- rst.iso.list[["25PF"]]
nominate_sGene_by_ACAT <- function(QTL_rst_mtx,     # should contain pval_beta
                                   phenotype_groups # should contain phenotype id and group (gene) id
){
  work_mtx <- left_join(QTL_rst_mtx, phenotype_groups, by = "phenotype_id")
  rst<- work_mtx %>% group_by(gene_id) %>% mutate(pval_aggr_ACAT = ACAT(pval_beta))
  rst <- distinct(rst, gene_id, .keep_all = TRUE)
  rst$FDR <- p.adjust(rst$pval_aggr_ACAT, method = "BH")
  qobj <- qvalue(rst$pval_aggr_ACAT)
  rst$qval_ACAT <- qobj$qvalues
  return(rst[,c("gene_id","pval_aggr_ACAT", "FDR", "qval_ACAT")])
}

test <- lapply(rst.iso.list, function(x) nominate_sGene_by_ACAT(x,phenotype_groups))
sig_sGene <- sapply(test, function(x){
  rst <- subset(x, qval_ACAT < 0.05)
  return(nrow(rst))
})

sGene.list <- sapply(test, function(x){
  rst <- subset(x, qval_ACAT < 0.05)
  return(unique(rst$gene_id))
})
df_sGene <- data.frame(covs = names(sig_sGene),
                       sGenes = sig_sGene)

df_sGene$covs <- factor(df_sGene$covs, levels = paste0(seq(5,60,5), "PF"))
ggplot(df_sGene, aes(x=covs, y=sGenes, fill = covs)) +
  geom_bar(stat="identity") + xlab("covs")+ylab("Number of sGenes") +
  geom_text(aes(label = sGenes), vjust=-1)+
  theme_classic()
df_sum <- data.frame(covs = names(rst.sig.iso.list),
                     Isoforms = sapply(rst.sig.iso.list, nrow))

df_sum$covs <- factor(df_sum$covs, levels = paste0(seq(5,60,5), "PF"))
ggplot(df_sum, aes(x=covs, y=Isoforms, fill = covs)) +
  geom_bar(stat="identity") + xlab("covs")+ylab("Number of significant Isoforms") +
  geom_text(aes(label = Isoforms), vjust=-1)+scale_fill_manual(values = alpha(scale_fill_hue()$palette(12),.7))+
  theme_classic()

################################
## Check p value in Q-Q plot
################################

all_qtl <- read.table(file = "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/output/AT2/qtl_results_all.txt", 
                      sep = "\t", header = TRUE)
permutated_qtl <- read.table(file = "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/output/AT2/permutation_qtl_results_all.txt", 
                             sep = "\t", header = TRUE)

qqnorm(-log10(all_qtl$pval_nominal), pch = 1, frame = FALSE)
qqline(-log10(rst.iso.list[["40PF"]]$pval_nominal), col = "steelblue", lwd = 2)
exp = -(log10(c(1:length(all_qtl$pval_nominal))/(length(all_qtl$pval_nominal)+ 1)))
exp = -(log10(c(1:10000)/(10000 + 1)))
obs = -log10(sort(all_qtl$pval_nominal))

set.seed(12345)
idx <- sample(1:length(all_qtl$pval_nominal), 20000000)
pvals <- all_qtl$pval_nominal[idx]
pvals_permutated <- permutated_qtl$pval_nominal[idx]
idx <- sample(1:nrow(mashr_1M_qtl@assays@data$lfsrs), 20000000)
pvals <- mashr_1M_qtl@assays@data$lfsrs[idx,"AT2"]
idx <- sample(1:nrow(mashr_p02_1M_qtl@assays@data$lfsrs), 20000000)
pvals <- mashr_p02_1M_qtl@assays@data$lfsrs[idx,"AT2"]

min(pvals[pvals != 0])
pvals[pvals == 0] = 1e-164
pvals_permutated <- mashr_permutation_qtl@assays@data$lfsrs[idx, "AT2"]
min(pvals_permutated[pvals_permutated != 0])
pvals_permutated[pvals_permutated == 0] = 1e-7
obs = -log10(sort(pvals))
exp = -(log10(c(1:length(pvals))/(length(pvals)+ 1)))
library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")

qqplotdata <- function(logpvector){
  o = sort(logpvector,decreasing=T)
  e = -log10(ppoints(length(o)))       
  qqdata <- data.frame(o,e)
  qqdata$o <- round(qqdata$o,3)
  qqdata$e <- round(qqdata$e,3)
  keepU <- which(!duplicated(qqdata))
  qqdata <- qqdata[keepU,]
  
  N <- length(logpvector) ## number of p-values
  ## create the confidence intervals
  qqdata$c975 <- NA
  qqdata$c025 <- NA
  
  ## the jth order statistic from a
  ## uniform(0,1) sample
  ## has a beta(j,n-j+1) distribution
  ## (Casella & Berger, 2002,
  ## 2nd edition, pg 230, Duxbury)
  
  for(i in 1:length(keepU)){
    j <- keepU[i]
    qqdata$c975[i] <- -log10(qbeta(0.975,j,N-j+1))
    qqdata$c025[i] <- -log10(qbeta(0.025,j,N-j+1))
  }
  return(qqdata)
}

qqplotdata_permutation <- function(logpvector, logpermutatedp){
  o = sort(logpvector,decreasing=T)
  e = logpermutatedp[order(logpvector,decreasing=T)]
  qqdata <- data.frame(o,e)
  qqdata$o <- round(qqdata$o,3)
  qqdata$e <- round(qqdata$e,3)
  keepU <- which(!duplicated(qqdata))
  qqdata <- qqdata[keepU,]
  
  N <- length(logpvector) ## number of p-values
  ## create the confidence intervals
  qqdata$c975 <- NA
  qqdata$c025 <- NA
  
  ## the jth order statistic from a
  ## uniform(0,1) sample
  ## has a beta(j,n-j+1) distribution
  ## (Casella & Berger, 2002,
  ## 2nd edition, pg 230, Duxbury)
  
  for(i in 1:length(keepU)){
    j <- keepU[i]
    qqdata$c975[i] <- -log10(qbeta(0.975,j,N-j+1))
    qqdata$c025[i] <- -log10(qbeta(0.025,j,N-j+1))
  }
  return(qqdata)
}

# yLine <- c(-log10(5E-8))
# colLine <- c("red")
dat$log10P = -log10(dat$P)
gwas = as.data.frame(dat)
# # Determine frequency bins and create variable for binned QQ plot
# 
# minMAF <- min(gwas$MAF)
# 
# freqbins <- c(c(0.5,0.05,0.005,0.001,0)[which(c(0.5,0.05,0.005,0.001,0) > floor(minMAF*1000000)/1000000)],floor(minMAF*1000000)/1000000)
# gwas$freqbin <- cut(gwas$MAF, freqbins,include.lowest=T)
# freqtable <- table(gwas$freqbin)
# freqtable <- freqtable[order(-as.numeric(gsub("[\\[\\(](.+),.+","\\1",names(freqtable))))]
# freqtable <- freqtable[freqtable > 0]

## Generate QQ plot data by frequency bin
fbin <- character(0)
fN <- integer(0)
fx <- numeric(0)
fy <- numeric(0)
fcol <- character(0)
legendcol <- character(0)
conf <- list()
allcols <- brewer.pal(4,"Set1")
ycol <- "log10P"
for(f in 1:length(freqtable)){
  # fbin <- c(fbin,"AT2")
  # fsnps <- which(gwas$freqbin ==names(freqtable)[f])
  plotdata <- qqplotdata(-log10(pvals_permutated))
  plotdata <- qqplotdata_permutation(-log10(pvals), -log10(pvals_permutated))
  fN <- c(fN,freqtable[f])
  fx <- c(fx,plotdata$e)
  fy <- c(fy,plotdata$o)
  fcol <- c(fcol,rep(allcols[f],length(plotdata$o)))
  conf[[f]] <- data.frame('x'=c(plotdata$e,rev(plotdata$e)),
                          'y'=c(plotdata$c975,rev(plotdata$c025)))
  legendcol <- c(legendcol,allcols[f])
}
legendtext <- paste0("MAF=",fbin,"; N SNPs=",format(fN,big.mark=",",scientific=FALSE))
opt =  list(break.top = 10,
            top.size = 0.125)


png(filename = "/data/lib14/R/Rscripts/AT2_qqplot.png", width = 8, height = 8, units = "in",res=300)
xlim <- c(0,max(fx,na.rm=T))
ylim <- c(0,max(fy,na.rm=T))
maxY <- max(fy,na.rm=T)
print("okkkk2")
par(mar=c(5.1,5.1,4.1,1.1))
print("okkkk3")

lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
lab2 <- lab2[lab2 > max(lab1)]

# resulting range of top scale in bottom scale units
top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
top.data = max(lab2)-opt$break.top

# function to rescale the top part
rescale = function(y) { opt$break.top+(y-opt$break.top)/(top.data/top.range)}
rescaled.y = rescale(fy[fy>opt$break.top])
plot(0,0,
     ylim=c(min(fy),opt$break.top*(1+opt$top.size)),xlim=xlim,axes=FALSE,
     xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
     ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
     cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
     main=opt$maintitle,pch=19)

# Plot confidence intervals	
for(p in 1:length(conf)){
  polygon(conf[[p]]$'x',ifelse(conf[[p]]$'y'>opt$break.top,rescale(conf[[p]]$'y'),conf[[p]]$'y'),
          col=grDevices::rgb(t(grDevices::col2rgb(allcols[p])),alpha=50,maxColorValue=255),
          border = NA)
}

# add points
points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)

# identify line & add axis break
lines(xlim,xlim,col="black",lty = 2)
axis(1,cex.axis=1.5,cex.lab=1.5)
par(las=1)
axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
par(las=0)
box()
par(las=0)
points(fx[fy>opt$break.top],rescaled.y,cex=1,col=fcol[fy>opt$break.top],pch=19)
par(las=1)
axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
lines(range(fx),c(opt$break.top,opt$break.top),col = "grey",lty = 6)
# abline(h=ifelse(yLine<opt$break.top,
#                 yLine,
#                 rescale(yLine)),
#        col=colLine,lwd=1.5,lty=2)
# legend("topleft",legend=legendtext,col=legendcol,pch=15,bty="n")
text(5,1,expression(paste(lambda[1000]," = ")),cex = 1.5)

x = pvals
z = qnorm(x / 2)
lambda = round(median(z^2) / qchisq(0.5,1), 3)
N.effect = 126
lambda_1000 = round(1+1000*(lambda-1)/N.effect  ,3)
text(5.7,1,paste(lambda_1000),cex = 1.5)

title(paste0(trait_name[i2]," for ",eth_name[i1]))
dev.off()

# plot(1,1)
# text(1,0.8,expression(paste(lambda[1000]," = ",buquote(.(lambda_1000)))),cex = 1.5)

lambda_vec = c(lambda,lambda_1000)





obs = -log10(sort(rst.iso.list[["40PF"]]$pval_beta))
obs = -log10(sort(rst.iso.list[["40PF"]]$pval_perm))
set.seed(1)
idx <- sample(1:length(all_qtl$pval_nominal), 40000)
pdf('/data/lib14/R/Rscripts/qqplot.pdf')
plot(exp, obs, ylab="Observed(−logP)",xlab="Expected(−logP)",ylim=c(0,50),xlim=c(0,6),main = "Nominal P value") 
plot(exp, obs, ylab="Observed(−logP)",xlab="Expected(−logP)",ylim=c(0,50),xlim=c(0,4),main = "Permutation-adjusted P value (beta approximation)") 
plot(exp, obs, ylab="Observed(−logP)",xlab="Expected(−logP)",ylim=c(0,5),xlim=c(0,4),main = "Permutation-adjusted P value (direct)") 

lines(c(0,6),c(0,6),col=1,lwd=2)
dev.off()
length(intersect(unique(subset(work_mtx, qval < 0.05)$gene_id), sGene.list$`25PF`))

Search_transcript_name(setdiff(iso30_ratio$phenotype_id, iso30_expr$phenotype_id))



#################
library(snpStats)
fam <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.fam"
bim <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bim"
bed <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bed"
sample <- read.plink(bed, bim, fam)
AT1 <- read.table(gzfile("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT1/phenotype.bed.gz"), header = TRUE, comment.char = "", sep = "\t")
rownames(iso30_ratio) <- iso30_ratio$phenotype_id
rownames(iso30_expr) <- iso30_expr$phenotype_id

snp_ratio <- iso30_ratio[intersect(iso30_ratio$phenotype_id, iso30_expr$phenotype_id),]
snp_expr <- iso30_expr[intersect(iso30_ratio$phenotype_id, iso30_expr$phenotype_id),]
tmp <- left_join(snp_ratio, snp_expr, by = "phenotype_id")
which(tmp$variant_id.x == tmp$variant_id.y)
tmp <- tmp[which(tmp$variant_id.x == tmp$variant_id.y),]
which(tmp$slope.x*tmp$slope.y > 0)
# tmp$transcript_name <- Search_transcript_name(tmp$phenotype_id)

# rs2070687-SFTPC-TALONT000730493
# rs2070686-SFTPC-202-ENST00000437090
isoQTLs <- c("rs2070686", "rs2070687")
isoQTLs <- c("rs2710880")
isoQTLs <- c("rs9786995")
genotype_mtx <- sample$genotypes[,isoQTLs]
snps_mtx <- genotype_mtx@.Data
snps_mtx <- snps_mtx[colnames(count_mtx_sum),]
snps_mtx <- snps_mtx[colnames(AT1)[c(5:130)],]
isoQTL_sum <- data.frame(rs2070686 = snps_mtx[,1],
                         rs2070687 = snps_mtx[,2],
                         SFTPC_202_ENST00000437090_ratio = count_mtx_per["ENST00000437090",],
                         SFTPC_202_ENST00000437090_expr = count_mtx_sum["ENST00000437090",],
                         SFTPC_TALONT000730493_ratio = count_mtx_per["TALONT000730493",],
                         SFTPC_TALONT000730493_expr = count_mtx_sum["TALONT000730493",])
isoQTL_sum <- data.frame(rs2710880 = snps_mtx,
                         ENST00000304952 = as.numeric(AT1[which(AT1$trascript_id == "ENST00000304952"),c(5:130)]))
which(isoQTL_sum$rs2070686 == "01")
isoQTL_sum$rs2710880_AG <- "0|0"
isoQTL_sum$rs2710880_AG[which(isoQTL_sum$rs2710880 == "01")] <- "1|1"
isoQTL_sum$rs2710880_AG[which(isoQTL_sum$rs2710880 == "02")] <- "0|1"

isoQTL_sum$rs2070686_GT <- "0|0"
isoQTL_sum$rs2070686_GT[which(isoQTL_sum$rs2070686 == "01")] <- "1|1"
isoQTL_sum$rs2070686_GT[which(isoQTL_sum$rs2070686 == "02")] <- "0|1"
table(isoQTL_sum$rs2070686_GT)
isoQTL_sum$rs2070687_GT <- "0|0"
isoQTL_sum$rs2070687_GT[which(isoQTL_sum$rs2070687 == "01")] <- "1|1"
isoQTL_sum$rs2070687_GT[which(isoQTL_sum$rs2070687 == "02")] <- "0|1"
table(isoQTL_sum$rs2070687_GT)
ggplot(isoQTL_sum, aes(x=rs2070686_GT, y=SFTPC_202_ENST00000437090_ratio, fill=rs2070686_GT)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="rs2070686", y = "Ratio of SFTPC-202")+
  theme_classic()
ggplot(isoQTL_sum, aes(x=rs2070686_GT, y=SFTPC_202_ENST00000437090_expr, fill=rs2070686_GT)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="rs2070686", y = "Normalized expression of SFTPC-202")+
  theme_classic()
ggplot(isoQTL_sum, aes(x=rs2710880_AG, y=ENST00000304952, fill=rs2710880_AG)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="rs2710880", y = "Normalized expression of ENST00000304952")+
  theme_classic()
ggplot(isoQTL_sum, aes(x=rs2070687_GT, y=SFTPC_TALONT000730493_ratio, fill=rs2070687_GT)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="rs2070687", y = "Ratio of SFTPC-TALONT000730493")+
  theme_classic()
ggplot(isoQTL_sum, aes(x=rs2070687_GT, y=SFTPC_TALONT000730493_expr, fill=rs2070687_GT)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="rs2070687", y = "Normalized expression of SFTPC-TALONT000730493")+
  theme_classic()















# eIsoform list
rst.list <- c(rst.iso.ct.list, rst.iso.list)
rst.sig.list <- c(rst.sig.iso.ct.list, rst.sig.iso.list)
eIso.list <- rst.list[c(3,19,257,29,38,59,64,76,85,89,102,111,122,127,131,
                        140,141,150,160,176,179,190,196,207,208,222,226,
                        239,252)]
eIso.sig.list <- rst.sig.list[c(3,19,314,29,38,59,64,82,94,103,116,125,136,144,
                                156,165,166,175,184,194,210,213,228,233,254,264,
                                254,264,265, 279,283,295,309)]
celltypes <- c("Alveolar transitional cells", "AT1", "AT2", "Club", "Goblet", "Secretory transitional cells", "Basal", "Multiciliated",celltypes)
names(eIso.list) <- celltypes
names(eIso.sig.list) <- celltypes
df_sum <- data.frame(celltype = names(epi_sig_list),
                     eIsoform_mashr_300K = sapply(epi_p02_300K_mashr_sig_list[gsub(" ", "_", names(epi_sig_list))], function(x){length(unique(x$phenotype))}),
                     eIsoform_mashr_100K = sapply(epi_p02_mashr_sig_list[gsub(" ", "_", names(epi_sig_list))], function(x){length(unique(x$phenotype))}),
                     eIsoform_tensorqtl= sapply(epi_sig_list, nrow),
                     eIsoform_overlapped_300K_lfsr.1 = overlapped,
                     eIsoform_overlapped_300K = overlapped1,
                     eIsoform_overlapped_100K = overlapped2)
df_sum <- data.frame(celltype = names(eIso.sig.list),
                     eIsoform_tensorqtl= sapply(eIso.sig.list, nrow))
df_sum$celltype <- factor(df_sum$celltype, levels = celltypes)
ggplot(df_sum, aes(x=celltype, y=eIsoform_tensorqtl, fill = celltype)) +
  geom_bar(stat="identity") + xlab("celltype")+ylab("Number of eIsoforms") +
  geom_text(aes(label = eIsoform_tensorqtl), vjust=-1)+theme_classic()+RotatedAxis()
df_sum <- data.frame(celltype = names(epi_sig_list),
                     eIsoform_tensorqtl= sapply(epi_sig_list, nrow),
                     eIsoform_mashr_1M = sapply(epi_p02_1M_mashr_sig_list[gsub(" ", "_", names(epi_sig_list))], function(x){length(unique(x$phenotype))}),
                     eIsoform_mashr_1M_all = sapply(epi_all_1M_mashr_sig_list[gsub(" ", "_", names(epi_sig_list))], function(x){length(unique(x$phenotype))}),
                     eIsoform_mashr_300K = sapply(epi_p02_300K_mashr_sig_list[gsub(" ", "_", names(epi_sig_list))], function(x){length(unique(x$phenotype))}),
                     eIsoform_mashr_100K = sapply(epi_p02_mashr_sig_list[gsub(" ", "_", names(epi_sig_list))], function(x){length(unique(x$phenotype))}),
                     eIsoform_mashr_200K = sapply(epi_p02_200K_mashr_sig_list[gsub(" ", "_", names(epi_sig_list))], function(x){length(unique(x$phenotype))}),
                     eIsoform_mashr_10K = sapply(epi_p02_10K_mashr_sig_list[gsub(" ", "_", names(epi_sig_list))], function(x){length(unique(x$phenotype))}),
                     eIsoform_mashr_50K = sapply(epi_p02_50K_mashr_sig_list[gsub(" ", "_", names(epi_sig_list))], function(x){length(unique(x$phenotype))}))
df_sum <- data.frame(celltype = names(epi_sig_list),
                     eIsoform_tensorqtl= sapply(epi_sig_list, nrow),
                     eGene = c(283,403,494,387,253,268,176,746),
                     eIsoform_mashr = sapply(epi_p02_mashr_sig_list[gsub(" ", "_", names(epi_sig_list))], function(x){length(unique(x$phenotype))}))

c(283,403,494,387,253,268,176,746)


df_sum$eIsoform_300K_ratio = round((df_sum$eIsoform_overlapped_300K/df_sum$eIsoform_tensorqtl)*100,2)
df_sum$eIsoform_100K_ratio = round((df_sum$eIsoform_overlapped_100K/df_sum$eIsoform_tensorqtl)*100,2)
df_sum$eIsoform_300K_lfsr.1_ratio = round((df_sum$eIsoform_overlapped_300K_lfsr.1/df_sum$eIsoform_tensorqtl)*100,2)
pdf('R/Rscripts/Epi_iso_sum.pdf', width = 10)
ggplot(df_sum, aes(x=celltype, y=num_isoform, fill = celltype)) +
  geom_bar(stat="identity") + xlab("celltype")+ylab("Number of Isoform (>20%)") +
  geom_text(aes(label = num_isoform), vjust=-1)+theme_classic()
ggplot(df_sum, aes(x=celltype, y=eIsoform_tensorqtl, fill = celltype)) +
  geom_bar(stat="identity") + xlab("celltype")+ylab("Number of eIsoform in long-read") +
  geom_text(aes(label = eIsoform_tensorqtl), vjust=-1)+theme_classic()+
  theme(axis.title = element_text(size = 14, colour = "black"))
ggplot(df_sum, aes(x=celltype, y=eGene, fill = celltype)) +
  geom_bar(stat="identity") + xlab("celltype")+ylab("Number of eGene in short-read") +
  scale_fill_manual(values = alpha(scale_fill_hue()$palette(8),.5))+
  geom_text(aes(label = eGene), vjust=-1)+theme_classic()
ggplot(df_sum, aes(x=celltype, y=eIsoform_mashr, fill = celltype)) +
  geom_bar(stat="identity") + xlab("celltype")+ylab("eIsoform by mashr") +
  scale_fill_manual(values = alpha(scale_fill_hue()$palette(12),.7))+
  geom_text(aes(label = eIsoform_mashr), vjust=-1)+theme_classic()

ggplot(df_sum, aes(x=celltype, y=overlapped, fill = celltype)) +
  geom_bar(stat="identity") + xlab("celltype")+ylab("shared eIsoform") +
  scale_fill_manual(values = alpha(scale_fill_hue()$palette(12),.5))+
  geom_text(aes(label = overlapped), vjust=-1)+theme_classic()
library(reshape2)
df_sum_long <- melt(df_sum)
df_sum_long$variable <- factor(df_sum_long$variable, c("eIsoform_tensorqtl","eIsoform_mashr_100K","eIsoform_mashr_300K",
                                                       "eIsoform_overlapped_100K","eIsoform_overlapped_300K","eIsoform_overlapped_300K_lfsr.1",
                                                       "eIsoform_100K_ratio","eIsoform_300K_ratio","eIsoform_300K_lfsr.1_ratio"))

df_sum_long$variable <- factor(df_sum_long$variable, c("eIsoform_tensorqtl","eIsoform_mashr_10K","eIsoform_mashr_50K","eIsoform_mashr_100K",
                                                       "eIsoform_mashr_200K","eIsoform_mashr_300K","eIsoform_mashr_1M","eIsoform_mashr_1M_all"))
library(ggbreak)

ggplot(subset(df_sum_long, variable!= "eIsoform_mashr_50K" ), 
       aes(x=celltype, y=value, fill = variable)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("covs")+ylab("Number of significant isoforms") +
  geom_text(aes(label = value), , vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_y_cut(breaks=c(600))+theme(axis.text = element_text(size = 14, colour = "black"))+
  theme_classic()
ggplot(subset(df_sum_long, variable %in% c("eIsoform_mashr_300K","eIsoform_mashr_100K","eIsoform_tensorqtl")), 
       aes(x=celltype, y=value, fill = variable)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("covs")+ylab("Number of significant isoforms") +
  geom_text(aes(label = value), , vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_y_cut(breaks=c(600))+theme(axis.text = element_text(size = 14, colour = "black"))+
  theme_classic()
ggplot(subset(df_sum_long, variable %in% c("eIsoform_tensorqtl", "eIsoform_mashr")), 
       aes(x=celltype, y=value, fill = variable)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("Cell type")+ylab("Number of significant isoforms") +
  geom_text(aes(label = value), , vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+ guides(fill=guide_legend(title=""))+
  scale_y_cut(breaks=c(600))+theme(axis.text = element_text(size = 14, colour = "black"))+
  theme_classic()
ggplot(df_sum_long, 
       aes(x=celltype, y=value, fill = variable)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("Cell type")+ylab("Number of significant isoforms") +
  geom_text(aes(label = value), , vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_y_cut(breaks=c(600))+theme(axis.text = element_text(size = 14, colour = "black"))+
  theme_classic()
ggplot(subset(df_sum_long, variable %in% c("eIsoform_overlapped_100K","eIsoform_overlapped_300K","eIsoform_overlapped_300K_lfsr.1")), 
       aes(x=celltype, y=value, fill = variable)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("covs")+ylab("Number of overlapped significant isoforms") +
  geom_text(aes(label = value), , vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  theme(axis.text = element_text(size = 14, colour = "black"))+
  theme_classic()+RotatedAxis()
ggplot(subset(df_sum_long, variable %in% c("eIsoform_100K_ratio","eIsoform_300K_ratio","eIsoform_300K_lfsr.1_ratio")), 
       aes(x=celltype, y=value, fill = variable)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("Cell type")+ylab("Percentage of overlapped significant isoforms (%)") +
  geom_text(aes(label = value), , vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  theme(axis.text = element_text(size = 14, colour = "black"))+
  theme_classic()+RotatedAxis()
dev.off()
saveRDS(eIso.sig.list,"/data/Choi_lung/scLongreads/tensorqtl/isoform_level/epi_sig_list.rds")
epi_p02_300K_mashr_sig_list <- epi_p02_300K_mashr_sig_list[gsub(" ", "_", names(epi_sig_list))]
epi_p02_1M_mashr_sig_list <- epi_p02_1M_mashr_sig_list[gsub(" ", "_", names(epi_sig_list))]
epi_p02_mashr_sig_list <- epi_p02_mashr_sig_list[gsub(" ", "_", names(epi_sig_list))]
mashr_rst.list <- mashr_rst.list[gsub(" ", "_", names(epi_sig_list))]
overlapped <- sapply(1:8, function(x){
  mashriso <- unique(mashr_rst.list[[x]]$phenotype_id)
  tensoriso <- epi_sig_list[[x]]$phenotype_id
  length(intersect(mashriso, tensoriso))
})
overlapped1 <- sapply(1:8, function(x){
  mashriso <- unique(epi_p02_1M_mashr_sig_list[[x]]$phenotype_id)
  tensoriso <- epi_sig_list[[x]]$phenotype_id
  length(intersect(mashriso, tensoriso))
})
overlapped2 <- sapply(1:8, function(x){
  mashriso <- unique(epi_p02_mashr_sig_list[[x]]$phenotype_id)
  tensoriso <- epi_sig_list[[x]]$phenotype_id
  length(intersect(mashriso, tensoriso))
})

#################################
## Nominal isoQTL
library(arrow)

chr1 <- read_parquet('/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2/output_nominal/Standardized_quantiled_age_batch_iso20_expr_30PF.cis_qtl_pairs.chr1.parquet')
rst.list <- c(rst.iso.ct.list, rst.iso.list)
shared_isos <- Reduce(intersect, lapply(rst.list, function(x)x$phenotype_id))

epi_p02_mashr_sig_list <- lapply(epi_p02_mashr_sig_list, function(x){
  x <- x %>% group_by(phenotype_id) %>% mutate(transcript_name = Search_transcript_name2(phenotype_id),
                                               gene_name = Search_transcript_gene_name(phenotype_id))
})
AT2_mashr <- epi_p02_mashr_sig_list$AT2
sig_isoGene_but_AT2 <- Reduce(union, lapply(epi_p02_mashr_sig_list[-3], function(x)unique(x$gene_name)))
sig_isoform_but_AT2 <- Reduce(union, lapply(epi_p02_mashr_sig_list[-3], function(x)unique(x$phenotype_id)))
sig_epi_isoform <- Reduce(union, lapply(epi_p02_mashr_sig_list, function(x)unique(x$gene_name)))
sig_epi_isoform <- Reduce(union, lapply(epi_p02_mashr_sig_list, function(x)unique(x$phenotype_id)))
AT2_tensor <- epi_sig_list$AT2

GTEx_sGene <- read.table(gzfile("/data/Choi_lung/scLongreads/GTEx/GTEx_Analysis_v8_sQTL/Lung.v8.sgenes.txt.gz"),sep="\t", header = TRUE)
AT2_tensor$transcript_name <- Search_transcript_name2(AT2_tensor$phenotype_id)
AT2_tensor$gene_name <- Search_transcript_gene_name(AT2_tensor$phenotype_id)

AT2_mashr <- AT2_mashr %>% group_by(phenotype_id) %>% mutate(transcript_name = Search_transcript_name2(phenotype_id),
                                                             gene_name = Search_transcript_gene_name(phenotype_id))
AT2_specific_isoGene <- setdiff(unique(AT2_mashr$gene_name), sig_isoGene_but_AT2)
AT2_specific_isoform <- unique(Search_transcript_gene_name(setdiff(unique(AT2_mashr$phenotype_id), sig_isoform_but_AT2)))
length(setdiff(AT2_mashr$gene_name, GTEx_sGene$gene_name))
length(setdiff(sig_epi_isoform, GTEx_sGene$gene_name))
intersect(AT2_specific_isoGene, setdiff(AT2_mashr$gene_name, GTEx_sGene$gene_name))

df_isoQTL <- data.frame(Aggregation = rep(c("sum w/ cn", "sum", "mean", "ratio"), each = 12),
                        covs_pf = rep(paste0(c(10,15,20,25,30,35,40,45,50,55,5,60),"PF"), 4),
                        Sig = sapply(rst.sig.iso.list, nrow))
df_isoQTL$Aggregation <- factor(df_isoQTL$Aggregation, levels = c("ratio", "mean", "sum", "sum w/ cn"))
ggplot(subset(df_isoQTL, covs_pf == "30PF"), aes(x=Aggregation, y=Sig, fill = Aggregation)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("covs")+ylab("Number of significant isoforms") +
  geom_text(aes(label = Sig), , vjust=-.3, color="black",
            position = position_dodge(0.9), size=3.5)+theme_bw()# +
theme(axis.text.x = element_text(size = 16, colour = "black"))



# Natri et al
df_sum <- data.frame(celltype = names(Top_sig_eQTL),
                     eGenes = sapply(Top_sig_eQTL, nrow))
library(ggplot2)
pdf("/data/lib14/project/scLongread/Natri_epi_eGenes.pdf", width = 10)
ggplot(df_sum, 
       aes(x=celltype, y=eGenes, fill = celltype)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("Cell type")+ylab("Number of eGenes") +
  geom_text(aes(label = eGenes), , vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  theme_classic(base_size = 14)+RotatedAxis()
dev.off()

natri_list <- lapply(Top_sig_eQTL, function(x)unique(x$feature_id))
Natri_eGenes <- Reduce(union, natri_list) # 321 eGenes in total


tensor_stat <- read.table("Tensor_stat.txt", sep = "\t")

ggplot(tensor_stat,aes(x=V1,y = V3))+geom_point() +
  stat_smooth(method=lm,formula=y~x) + 
  geom_text_repel(label=tensor_stat$V2)+
  stat_cor(label.y = 400,aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  ylab("Number of eIsoform")+
  xlab("Number of input")+theme_bw()

tensor_stat$V2 <- factor(tensor_stat$V2, levels = tensor_stat$V2)
tensor_stat$Categroy <- c(rep("EPI",8),rep("IMMUNE",15), rep("EC", 6), rep("STROMA",4))
tensor_stat$Categroy <- factor(tensor_stat$Categroy, levels = c("EPI", "IMMUNE", "EC", "STROMA"))
ggplot(tensor_stat, aes(x=V2, y=V1, fill = Categroy)) +
  geom_bar(stat="identity") + xlab("celltype")+ylab("Number of input isoforms") +
  geom_text(aes(label = V1), vjust=-1)+theme_classic()+RotatedAxis()
ggplot(tensor_stat, aes(x=V2, y=V3, fill = Categroy)) +
  geom_bar(stat="identity") + xlab("celltype")+ylab("Number of eIsoform") +
  geom_text(aes(label = V3), vjust=-1)+theme_classic()+RotatedAxis()

tmp.list <- lapply(eIso.sig.list, function(x) x$phenotype_id)
tmp <- Reduce(c, tmp.list)
length(unique(tmp))
write.table(tmp, "Tensor.sig.txt", sep = "\t")
write.table(gsub(" ", "_",celltypes), "Celltype_nPF", sep = "\t", quote = FALSE, row.names = FALSE)



rst.sig.iso.ct.list <- list()
rst.iso.ct.list <- list()
for (celltype in celltypes) {
  celltype_for_path <- gsub(" ", "_", celltype)
  output_path <- paste0("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/", celltype_for_path, "/output/")
  file_paths <- list.files(output_path, 
                           pattern = "cis_qtl.txt.gz")
  file_paths <- file_paths[grepl("Standardized_quantiled_expr_age_batch_3gpc", file_paths)]
  for (i in c(1:length(file_paths))) {
    file_path <- paste0(output_path, file_paths[i])
    covs_used <- str_split_fixed(str_split_fixed(file_paths[i], "_", n = 5)[,5], "\\.", n = 2)[,1]
    covs_used <- paste(celltype_for_path, covs_used, sep = "-")
    rst <- read.table(gzfile(file_path),sep="\t", header = TRUE)
    rst.iso.ct.list[[covs_used]] <- rst
    rst <- subset(rst, qval < 0.05)
    rst.sig.iso.ct.list[[covs_used]] <- rst
  }
}





rst.sig.iso.ct.list <- list()
rst.iso.ct.list <- list()
for (celltype in celltypes) {
  celltype_for_path <- gsub(" ", "_", celltype)
  output_path <- paste0("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/", celltype_for_path, "/output/")
  file_paths <- list.files(output_path, 
                           pattern = "cis_qtl.txt.gz")
  file_paths <- file_paths[grepl("Standardized_quantiled_expr_age_batch_gpc3", file_paths)]
  file_path <- paste0(output_path, file_paths)
  rst <- read.table(gzfile(file_path),sep="\t", header = TRUE)
  rst.iso.ct.list[[celltype_for_path]] <- rst
  rst <- subset(rst, qval < 0.05)
  rst.sig.iso.ct.list[[celltype_for_path]] <- rst
}
eIso.list <- rst.iso.ct.list[c(4,17,23,29,42,53,58,68,77,87,94,101,114,124,129,141,
                               153,154,168,175,181,193,200,213,218,233,240,251,
                               253,266,271,284,292)]
eIso.sig.list <- rst.sig.iso.ct.list[c(4,17,23,29,42,53,58,68,77,87,94,101,114,124,129,141,
                                       153,154,168,175,181,193,200,213,218,233,240,251,
                                       253,266,271,284,292)]
names(eIso.sig.list) <- celltypes
df_sum <- data.frame(celltype = names(eIso.sig.list),
                     eIsoform_tensorqtl= sapply(eIso.sig.list, nrow))
saveRDS(eIso.list, file = "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/isoQTL_list.rds")
saveRDS(eIso.sig.list, file = "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/isoQTL_sig_list.rds")

write.table(df_sum, "/data/lib14/R/Rscripts/Tensor.sig.txt", sep = "\t", quote = FALSE, row.names = FALSE)
eIso.sig.list <- readRDS("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/isoQTL_sig_list.rds")
mashr.eIso.sig.list <- readRDS("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/output/mashr_2M_sig_final_list.rds")
ct_specific_eIso <- NULL
ct_specific_eIso.list <- list()
for (i in c(1:33)) {
  ct <- eIso.sig.list[[i]]$phenotype_id
  rest.list <- lapply(eIso.sig.list[-i], function(x) x$phenotype_id)
  rest <- Reduce(union, rest.list)
  ct_iso <- setdiff(ct, rest)
  ct_specific_eIso.list[[i]] <- ct_iso
  nct <- length(setdiff(ct, rest))
  ct_specific_eIso <- c(ct_specific_eIso, nct)
}
ct_eIso <- Reduce(c,ct_specific_eIso.list)
length(unique(ct_eIso)) # 782 (71.22%)
all_eIso_tensor <- Reduce(union, lapply(eIso.sig.list, function(x) x$phenotype_id)) # 1098

df_sum <- data.frame(celltype = names(eIso.sig.list), 
                     ct_specific = ct_specific_eIso, 
                     eIsoform_tensorqtl= sapply(eIso.sig.list, nrow))
df_sum$pct_specific <- round((df_sum$ct_specific/df_sum$eIsoform_tensorqtl)*100, 2)
nshare.eIso <- c()
pctshare.eIso <- c()
for (i in c(1:33)) {
  NB <- NB.list[[i]]$phenotype_id[which(NB.list[[i]]$qval < 0.05)]
  tensor <- eIso.sig.list[[i]]$phenotype_id
  nshared <- length(intersect(NB, tensor))
  pctshared <- nshared/length(tensor)
  nshare.eIso <- c(nshare.eIso,nshared)
  pctshare.eIso <- c(pctshare.eIso,pctshared)
}

NB.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_list.rds")



df_sum$Categroy <- c(rep("EPI",8),rep("IMMUNE",15), rep("EC", 6), rep("STROMA",4))
df_sum$Categroy <- factor(df_sum$Categroy, levels = c("EPI", "IMMUNE", "EC", "STROMA"))
df_sum$celltype <- factor(df_sum$celltype, levels = names(eIso.sig.list))
ggplot(df_sum, aes(x=eIsoform_tensorqtl, y=celltype, fill = Categroy)) +
  geom_bar(stat="identity") + xlab("Number of eIsoforms")+ylab("Cell type") +
  scale_fill_manual(values = c("#ffd700", "#fa8775","#cd34b5", "#0000ff"))+
  geom_text(aes(label = eIsoform_tensorqtl), hjust=-.1)+theme_classic()+RotatedAxis()
ggsave("/data/lib14/project/scLongread/celltype_eIsoform.pdf", width = 8, height = 8)

mashr.eIso.top.list <- lapply(mashr.eIso.sig.list, function(x){
  result <- x %>% 
    group_by(phenotype_id) %>%
    filter(lfsr == max(lfsr)) %>%
    arrange(phenotype_id,variant_id,posterior_means,posterior_means,sds,lfsr,transcript_name,gene_name)
  result <- distinct(result, phenotype_id, .keep_all = TRUE)
  return(result)
})
ct_specific_eIso <- NULL
ct_specific_eIso.list <- list()
Tensor_sig_only <- list()
for (i in c(1:33)) {
  mashr.eIso.top.list[[i]]$Celltype <- names(mashr.eIso.top.list)[i]
  mashr.eIso.top.list[[i]]$CT_specific <- FALSE
  mashr.eIso.top.list[[i]]$Sig_in_tensor <- FALSE
  ct <- mashr.eIso.top.list[[i]]$phenotype_id
  rest.list <- lapply(mashr.eIso.top.list[-i], function(x) x$phenotype_id)
  rest <- Reduce(union, rest.list)
  ct_iso <- setdiff(ct, rest)
  ct_specific_eIso.list[[i]] <- ct_iso
  mashr.eIso.top.list[[i]]$CT_specific[which(mashr.eIso.top.list[[i]]$phenotype_id %in% ct_iso)] <- TRUE
  tensor_eIso <- eIso.sig.list[[i]]$phenotype_id
  mashr.eIso.top.list[[i]]$Sig_in_tensor[which(mashr.eIso.top.list[[i]]$phenotype_id %in% tensor_eIso)] <- TRUE
  Tensor_sig_only[[i]] <- setdiff(tensor_eIso, ct)
  nct <- length(setdiff(ct, rest))
  ct_specific_eIso <- c(ct_specific_eIso, nct)
}
mashr.eIso.top <- Reduce(rbind, mashr.eIso.top.list)
write.table(mashr.eIso.top, file = "/data/lib14/project/scLongread/Tables/Mashr_eIsoform_sum_top_ver.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


lung_GWAS <- read.table("/data/lib14/R/Rscripts/Lung_GWAS_51loci.txt", 
                        header = TRUE, sep = "\t")
rownames(lung_GWAS) <- lung_GWAS$RSNUM
Tensor_eIso_top <- read.table(file = "/data/lib14/project/scLongread/Tables/TensorQTL_eIsoform_sum.txt",
                              sep = "\t", header = TRUE)

Tensor_eIso_top_not_Mashr <- subset(Tensor_eIso_top, Sig_in_Mashr == FALSE)

table(Tensor_eIso_top_not_Mashr$Celltype.x)

file_paths <- list.dirs("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/output",)
file_paths <- file_paths[-1]
Tensor_all_qtl_list <- list()
for (ct_path in file_paths) {
  file_path <- paste0(ct_path, "/top_qtl_results_all_FDR0.05.txt")
  eQTL <- read.table(file_path, sep = "\t", header = TRUE)
  ct <- str_split_fixed(ct_path, "/", 8)[,8]
  lung_GWAS_sub <- lung_GWAS[,c(2,7)]
  colnames(lung_GWAS_sub) <- c("variant_id", "CCV")
  lung_GWAS_sub <- as.data.frame(lung_GWAS_sub)
  tmp1 <- left_join(eQTL, lung_GWAS_sub, by = "variant_id")
  tmp1$Celltype <- ct
  Tensor_all_qtl_list[[ct]] <- tmp1
}


Tensor_all_qtl <- Reduce(rbind, Tensor_all_qtl_list)
Tensor_sig_lungGWAS <- Tensor_all_qtl[!is.na(Tensor_all_qtl$CCV),]
length(unique(Tensor_sig_lungGWAS$variant_id))

intersect(lung_GWAS$RSNUM, Tensor_eIso_top$variant_id.x)
for (i in 1:33) {
  isoQTL_NB[[i]]$Celltype <- names(isoQTL_NB)[i]
  isoQTL_NB[[i]]$Celltype_specificity <- FALSE
  ct <- isoQTL_NB[[i]]$phenotype_id
  rest <- unique(Reduce(c, lapply(isoQTL_NB[-i], function(x) x$phenotype_id)))
  ct_specific <- setdiff(ct, rest)
  isoQTL_NB[[i]]$Celltype_specificity[which(isoQTL_NB[[i]]$phenotype_id %in% ct_specific)] <- TRUE
}
NB_eIso_top <- Reduce(rbind, isoQTL_NB)
intersect(lung_GWAS$RSNUM, NB_eIso_top$variant_id)

mashr_all <- readRDS("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/output/mashr_2M_seed15_qtl.rds")
Tensor_sig_only_all <- unique(Reduce(c, Tensor_sig_only))
View(Tensor_eIso_top[which(Tensor_eIso_top$phenotype_id.x %in% Tensor_sig_only_all),])
ids <- str_split_fixed(Tensor_eIso_top_not_Mashr$unique_id, "_", 2)[,1]
ids <- gsub(":", "|", ids)
posterior_means = mashr_all@assays@data$betas[ids,]
sds = mashr_all@assays@data$errors[ids,]
lfsr = mashr_all@assays@data$lfsrs[ids,]

