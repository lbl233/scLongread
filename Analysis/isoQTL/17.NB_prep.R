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

file_paths <- list.files("/data/Choi_lung/scLongreads/jaxqtl/AT2/permutation/TALONT002194518_seeds/",
                        pattern = "score.parquet")

file_paths <- file_paths[-11]
rst.list <- list()
for (file in file_paths) {
  file_path <- paste0("/data/Choi_lung/scLongreads/jaxqtl/AT2/permutation/TALONT002194518_seeds/",file)
  rst <- read_parquet(file_path)
  name <- str_split_fixed(file, "_", n = 4)[,3]
  rst.list[[name]] <- rst
}
# nominal mode
file_paths <- list.files("/data/Choi_lung/scLongreads/jaxqtl/AT2/",
                         pattern = "score.parquet")
file_paths <- file_paths[1:23]
chr.list <- list()
for (file in file_paths) {
  file_path <- paste0("/data/Choi_lung/scLongreads/jaxqtl/AT2/",file)
  rst <- read_parquet(file_path)
  name <- str_split_fixed(file, "_", n = 2)[,1]
  chr.list[[name]] <- rst
}

file_paths <- list.files("/data/Choi_lung/scLongreads/jaxqtl/AT2/",
                         pattern = "jaxqtl_nb.cis_score.tsv.gz")
chr.list <- list()
for (file in file_paths) {
  file_path <- paste0("/data/Choi_lung/scLongreads/jaxqtl/AT2/",file)
  rst <- read.table(gzfile(file_path), header = TRUE, sep = "\t")
  name <- str_split_fixed(file, "_", n = 2)[,1]
  chr.list[[name]] <- rst
}
chrs <- Reduce(rbind, chr.list)

rst_AT2 <- chrs
rst_perm <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/AT2/permutation/Nominal_permutation_AT2.rds")
rownames(rst_AT2) <- paste(rst_AT2$phenotype_id, rst_AT2$snp, sep = ":")
rownames(rst_AT2) <- paste(rst_AT2$phenotype_id, rst_AT2$variant_id, sep = ":")
rownames(rst_perm) <- paste(rst_perm$phenotype_id, rst_perm$snp, sep = ":")
rst_perm <- rst_perm[rownames(rst_AT2),]
pvals <- rst_AT2$pval_nominal
idx1 <- which(rst_AT2$converged == 1 & rst_AT2$alpha != 1e-8)
idx1 <- which(rst_AT2$beta_converged == 1 & rst_AT2$model_converged == 1 & rst_AT2$alpha_cov != 1e-8)
# pvals <- pvals[idx]

pvals_permutated <- rst_perm$pval_nominal
idx2 <- which(rst_perm$converged == 1 & rst_perm$alpha != 1e-8)
# pvals_permutated <- pvals_permutated[idx]
idx <- intersect(idx1, idx2)
pvals <- pvals[idx]
pvals_permutated <- pvals_permutated[idx]
library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")
set.seed(123)

idx <- sample(1:length(pvals), size = 100000, replace = FALSE)
pvals <- pvals[idx]
pvals_permutated <- pvals_permutated[idx]
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
  e = sort(logpermutatedp,decreasing=T)
  o = logpvector[order(logpermutatedp,decreasing=T)]
  # o = sort(logpvector,decreasing=T)
  # e = logpermutatedp[order(logpvector,decreasing=T)]
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

fbin <- character(0)
fN <- integer(0)
fx <- numeric(0)
fy <- numeric(0)
fcol <- character(0)
legendcol <- character(0)
conf <- list()
allcols <- brewer.pal(4,"Set1")
ycol <- "log10P"
fbin <- c(fbin,"AT2")
f = 1
plotdata <- qqplotdata(-log10(pvals))
plotdata <- qqplotdata(-log10(pvals_permutated))
plotdata <- qqplotdata_permutation(-log10(pvals), -log10(pvals_permutated))

fx <- c(fx,plotdata$e)
fy <- c(fy,plotdata$o)
fcol <- c(fcol,rep(allcols[f],length(plotdata$o)))
conf[[f]] <- data.frame('x'=c(plotdata$e,rev(plotdata$e)),
                        'y'=c(plotdata$c975,rev(plotdata$c025)))
legendcol <- c(legendcol,allcols[f])
opt =  list(break.top = 10,
            top.size = 0.125)


png(filename = "/data/lib14/R/Rscripts/AT2_obs_perm_qqplot.png", width = 8, height = 8, units = "in",res=300)
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
text(4,1,expression(paste(lambda[1000]," = ")),cex = 1.5)

x = pvals
z = qnorm(x / 2)
lambda = round(median(z^2) / qchisq(0.5,1), 3)
N.effect = 126
lambda_1000 = round(1+1000*(lambda-1)/N.effect  ,3)
text(4.7,1,paste(lambda_1000),cex = 1.5)

dev.off()

file_paths <- list.files("/data/Choi_lung/TTL/Files_for_Bolun_NB/Immune/Bcells/",
                         pattern = "jaxqtl_nb.cis_score.tsv.gz")
chr.list <- list()
for (file in file_paths) {
  file_path <- paste0("/data/Choi_lung/TTL/Files_for_Bolun_NB/Immune/Bcells/",file)
  rst <- read.table(gzfile(file_path), header = TRUE, sep = "\t")
  name <- str_split_fixed(file, "_", n = 2)[,1]
  chr.list[[name]] <- rst
}
sr_Bcells <- Reduce(rbind, chr.list)
sr_AT2 <- Reduce(rbind, chr.list)


bed <- read.table(gzfile("/data/Choi_lung/scLongreads/jaxqtl/AT2/phenotype.bed.gz"),
                  header = TRUE, comment.char = "", sep = "\t")
hist(as.numeric(bed[11686, c(5:130)]))
bed <- read.table(gzfile("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2//phenotype.bed.gz"),
                  header = TRUE, comment.char = "", sep = "\t")

expr <- bed[,c(5:130)]
sds <- apply(expr, 1, sd)
tmp <- apply(expr, 1, function(x){
  x <- as.numeric(x)
  #find Q3 and interquartile range for values in points column
  Q3 <- quantile(x, .75)
  IQR <- IQR(x)
  length(which(x > (Q3 + 3*IQR))) # test the threshold later
})

length(which(tmp < 8))


library(dplyr)
library(ggplot2)
library(stringr)
library(stringi)


file_paths <- list.files("/data/Choi_lung/scLongreads/jaxqtl/")
idx <- grep("swarm", file_paths)
file_paths <- file_paths[-idx]
file_paths <- file_paths[-c(9,12,25,26,27,35,37,39,40,43)]
NB.list <- list()
eIsoform <- c()
for (ct in file_paths) {
  files <- list.files(paste0("/data/Choi_lung/scLongreads/jaxqtl/",ct), 
                           pattern = "chunck_jaxqtl_nb.cis_score.tsv.gz")
  if(length(files) == 23){
    rst_1k <- NULL
    for (i in c(1:length(files))) {
      file_path <- paste0("/data/Choi_lung/scLongreads/jaxqtl/", ct, "/", files[i])
      rst <- read.table(gzfile(file_path),sep="\t", header = TRUE)
      rst_1k <- rbind(rst_1k, rst)
    }
    library(qvalue)
    rst_1k <- subset(rst_1k, (model_converged == 1 & alpha_cov != (1e-08)))
    qobj <- qvalue(rst_1k$pval_beta)
    rst_1k$qval <- qobj$qvalues
    message(ct, " model converged: ", nrow(rst_1k))
    message(ct, " eGenes: ",length(rst_1k$phenotype_id[which(rst_1k$qval < 0.05)]))
    eIsoform <- c(eIsoform, length(rst_1k$phenotype_id[which(rst_1k$qval < 0.05)]))
    NB.list[[ct]] <- rst_1k
  }else{
    message(ct, " not finished")
  }
}
df_sum <- data.frame(celltype = file_paths, # covar for sec_trans was not generated correctly
                     eIsoform = eIsoform)

isoQTL_tensor <- readRDS("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/isoQTL_sig_list.rds")
NB.list <- NB.list[gsub(" ", "_", names(isoQTL_tensor))]
df_sum <- data.frame(celltype = names(NB.list), # covar for sec_trans was not generated correctly
                     eIsoform = sapply(NB.list, function(x)length(x$phenotype_id[which(x$qval < 0.05)])))
df_sum$Categroy <- c(rep("EPI",8),rep("IMMUNE",15),rep("EC", 6), rep("STROMA",4))
df_sum$Categroy <- factor(df_sum$Categroy, levels = c("EPI", "IMMUNE", "EC", "STROMA"))
df_sum$celltype <- factor(df_sum$celltype, levels = df_sum$celltype)
ggplot(df_sum, aes(x=celltype, y=eIsoform, fill = Categroy)) +
  geom_bar(stat="identity") + xlab("Number of eIsoforms")+ylab("Cell type") +
  scale_fill_manual(values = c("#ffd700", "#fa8775","#cd34b5", "#0000ff"))+
  geom_text(aes(label = eIsoform), vjust=-.5)+theme_classic()+RotatedAxis()

NB.list <- NB.list[gsub(" ", "_", names(isoQTL_tensor))]
df_sum <- data.frame(celltype = names(NB.sig.list), # covar for sec_trans was not generated correctly
                     NB =  sapply(NB.sig.list, nrow),
                     tensor = sapply(isoQTL_tensor, nrow))
df_sum$Categroy <- c(rep("EPI",8),rep("IMMUNE",15),rep("EC", 6), rep("STROMA",4))
df_sum$Categroy <- factor(df_sum$Categroy, levels = c("EPI", "IMMUNE", "EC", "STROMA"))
df_sum$celltype <- factor(df_sum$celltype, levels = df_sum$celltype)
df_sum_long <- reshape2::melt(df_sum,value.name = "eIsoform")
df_sum_long$Categroy <- rep(c(rep("EPI",8),rep("IMMUNE",15),rep("EC", 6), rep("STROMA",4)), 2)
df_sum_long$Categroy <- factor(df_sum_long$Categroy, levels = c("EPI", "IMMUNE", "EC", "STROMA"))
df_sum_long$celltype <- factor(df_sum_long$celltype, levels = df_sum$celltype)

ggplot(df_sum_long, 
       aes(x=celltype, y=eIsoform, fill = variable)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("covs")+ylab("Number of significant isoforms") +
  geom_text(aes(label = eIsoform), , vjust=-.1, color="black",
            position = position_dodge(0.9), size=3)+
  theme(axis.text = element_text(size = 14, colour = "black"))+
  theme_classic()+RotatedAxis()


isoQTL_NB <- lapply(NB.list, function(x){
  x <- subset(x, qval < 0.05)
})


saveRDS(NB.list, file = "/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_list.rds")
saveRDS(isoQTL_NB, file = "/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")


nIso_shared <- c()
pctIso_shared <- c()
tensor_sig_nb_not <- list()
for (i in 1:length(isoQTL_NB)) {
  NB <- isoQTL_NB[[i]]$phenotype_id
  tensor <- isoQTL_tensor[[i]]$phenotype_id
  nIso_shared <- c(nIso_shared, length(intersect(NB,tensor)))
  tensor_sig_nb_not[[i]] <- (setdiff(tensor,NB))
  pctIso_shared <- c(pctIso_shared, length(intersect(NB,tensor))/length(tensor))
}

library(snpStats)
fam <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.fam"
bim <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bim"
bed <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bed"
sample <- read.plink(bed, bim, fam)

bed <- read.table(gzfile("/data/Choi_lung/scLongreads/jaxqtl/AT2/phenotype.bed.gz"),
                  header = TRUE, comment.char = "", sep = "\t")

bed <- read.table(gzfile("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/AT2//phenotype.bed.gz"),
                  header = TRUE, comment.char = "", sep = "\t")
isoQTLs <- c("rs2070686")
genotype_mtx <- sample$genotypes[,isoQTLs]
snps_mtx <- genotype_mtx@.Data
snps_mtx <- snps_mtx[colnames(bed)[c(5:130)],]

isoQTL_sum <- data.frame(rs2070686 = snps_mtx,
                         ENST00000318561 = as.numeric(bed[which(bed$trascript_id == "ENST00000318561"),c(5:130)]))
which(isoQTL_sum$rs2070686 == "01")
isoQTL_sum$rs2070686_TC <- "0|0"
isoQTL_sum$rs2070686_TC[which(isoQTL_sum$rs2070686 == "01")] <- "1|1"
isoQTL_sum$rs2070686_TC[which(isoQTL_sum$rs2070686 == "02")] <- "0|1"
table(isoQTL_sum$rs2070686_TC)
ggplot(isoQTL_sum, aes(x=rs2070686_TC, y=ENST00000318561, fill=rs2070686_TC)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="rs2070686", y = "Count expr of ENST00000318561")+
  theme_classic()
ggplot(isoQTL_sum, aes(x=rs2070686_TC, y=ENST00000318561, fill=rs2070686_TC)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="rs2070686", y = "Normalized expr of ENST00000318561")+
  theme_classic()

#########################################
### nominal p threshold calculation #####
#########################################
NB.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_list.rds")
for (i in 1:length(NB.list)) {
  NB.df <- NB.list[[i]]
  # determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
  ub <- sort(NB.df[NB.df$qval > 0.05, 'pval_beta'])[1]  # smallest p-value above FDR
  lb <- -sort(-NB.df[NB.df$qval <= 0.05, 'pval_beta'])[1]  # largest p-value below FDR
  pthreshold <- (lb+ub)/2
  cat("  * min p-value threshold @ FDR ", 0.05, ": ", pthreshold, "\n", sep="")
  NB.df[, 'pval_nominal_threshold'] <- signif(qbeta(pthreshold,
                                                    NB.df[, 'beta_shape1'], NB.df[, 'beta_shape2'], 
                                                    ncp=0, lower.tail=TRUE, log.p=FALSE), 6)
  NB.list[[i]] <- NB.df
}

saveRDS(NB.list, "/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_list.rds")
NB.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_list.rds")
# write.table(NB.df, gzfile(args$outfile), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
lung_GWAS <- read.table("/data/lib14/R/Rscripts/Lung_GWAS_51loci.txt", 
                        header = TRUE, sep = "\t")
rownames(lung_GWAS) <- lung_GWAS$RSNUM
library(arrow)
file_paths <- list.files("/data/Choi_lung/scLongreads/jaxqtl/")
idx <- grep("swarm", file_paths)
file_paths <- file_paths[-idx]
file_paths <- file_paths[-c(9,12,25:29,37,39,41,42,45)]

NB.nominal.list <- list()
NB.GWAS.list <- list()
for (ct in file_paths) {
  files <- list.files(paste0("/data/Choi_lung/scLongreads/jaxqtl/",ct), 
                      pattern = "wald.parquet")
  if(length(files) == 23){
    rst_1k <- NULL
    for (i in c(1:length(files))) {
      file_path <- paste0("/data/Choi_lung/scLongreads/jaxqtl/", ct, "/", files[i])
      rst <- read_parquet(file_path)
      rst_1k <- rbind(rst_1k, rst)
    }
    NB.nominal.list[[ct]] <- rst_1k
    rst_sub <- subset(rst_1k, snp %in% lung_GWAS$RSNUM)
    isoform_list <- unique(rst_sub$phenotype_id)
    rst_sub <- subset(rst_1k, phenotype_id %in% isoform_list)
    NB.GWAS.list[[ct]] <- rst_sub
  }else{
    message(ct, " not finished")
  }
}
saveRDS(NB.nominal.list, "/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_wald_list.rds")
saveRDS(NB.GWAS.list, "/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_wald_gwas_overlap_list.rds")


for (ct in file_paths) {
  files <- list.files(paste0("/data/Choi_lung/scLongreads/jaxqtl/",ct), 
                      pattern = "score.parquet")
  if(length(files) == 23){
    rst_1k <- NULL
    for (i in c(1:length(files))) {
      file_path <- paste0("/data/Choi_lung/scLongreads/jaxqtl/", ct, "/", files[i])
      rst <- read_parquet(file_path)
      rst_1k <- rbind(rst_1k, rst)
    }
    NB.nominal.list[[ct]] <- rst_1k
    cis <- NB.list[[ct]]
    cis <- subset(cis, qval < 0.05)
    rownames(cis) <- cis$phenotype_id
    tmp <- rst_1k %>% filter(phenotype_id %in% cis$phenotype_id) %>% 
      group_by(phenotype_id) %>% 
      mutate(significance = pval_nominal < cis[phenotype_id,"pval_nominal_threshold"])
    
    sig_qtl <- tmp[which(tmp$significance),]
    NB.nominal.sig.list[[i]] <- sig_qtl
    message(ct, " Number of SNPs overlapped with lung cancer GWAS signals: ", 
            length(intersect(unique(sig_qtl$snp), lung_GWAS$RSNUM)))
  }else{
    message(ct, " not finished")
  }
}
saveRDS(NB.nominal.list, "/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_list.rds")
saveRDS(NB.nominal.sig.list, "/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_sig_list.rds")
NB.nominal.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_list.rds")
NB.nominal.sig.list <- list()
for (ct in file_paths) {

    rst_1k <- NB.nominal.list[[ct]]
    cis <- NB.list[[ct]]
    cis <- subset(cis, qval < 0.05)
    rownames(cis) <- cis$phenotype_id
    tmp <- rst_1k %>% filter(phenotype_id %in% cis$phenotype_id) %>% 
      group_by(phenotype_id) %>% 
      mutate(significance = pval_nominal < cis[phenotype_id,"pval_nominal_threshold"])
    
    sig_qtl <- tmp[which(tmp$significance),]
    lung_GWAS_sub <- subset(lung_GWAS, RSNUM %in% unique(sig_qtl$snp))
    lung_GWAS_sub <- lung_GWAS_sub[,c(2,7)]
    colnames(lung_GWAS_sub) <- c("snp", "CCV")
    
    if(nrow(lung_GWAS_sub) == 0){
      sig_qtl$CCV <- NA
      message(ct, ": No lung GWAS signal")
      tmp1 <- sig_qtl
    }else{
      lung_GWAS_sub <- as.data.frame(lung_GWAS_sub)
      tmp1 <- left_join(sig_qtl, lung_GWAS_sub, by = "snp")
    }

    tmp1$Celltype <- ct
    NB.nominal.sig.list[[ct]] <- tmp1
    message(ct, ": ", 
            length(intersect(unique(sig_qtl$snp), lung_GWAS$RSNUM)))
}
for (i in 1:33) {
  rst <- NB.nominal.sig.list[[i]]
  file_path <- paste0("/data/lib14/project/scLongread/Tables/NB_all_sig_qtl_", names(NB.nominal.sig.list)[i], ".txt")
  write.table(rst, file_path, quote = FALSE, row.names = FALSE, sep = "\t")
}

NB_all_sig <- Reduce(rbind, NB.nominal.sig.list)
NB_sig_lungGWAS <- NB_all_sig[!is.na(NB_all_sig$CCV),]
length(unique(NB_sig_lungGWAS$snp))
length(unique(NB_sig_lungGWAS$phenotype_id))
NB_sig_lungGWAS$snp_specific <- FALSE

for (ct in unique(NB_sig_lungGWAS$Celltype)) {
  ct_tb <- subset(NB_sig_lungGWAS, Celltype == ct)
  rest <- subset(NB_sig_lungGWAS, Celltype != ct)
  ct_snp <- ct_tb$snp
  rest_snp <- rest$snp
  ct_snp_sp <- setdiff(ct_snp, rest_snp)
  NB_sig_lungGWAS$snp_specific[which(NB_sig_lungGWAS$snp %in% ct_snp_sp)] <- TRUE
}
final_table <- NB_sig_lungGWAS %>% group_by(phenotype_id) %>% mutate(transcript_name = Search_transcript_name2(phenotype_id),
                                             gene_name = Search_transcript_gene_name(phenotype_id), .after = phenotype_id)


NB.nominal.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_sig_list.rds")
isoQTLs_NB_List <- lapply(NB.nominal.sig.list, function(x)x$phenotype_id)
isoQTLs_NB_List <- unique(Reduce(c, isoQTLs_NB_List))

NB.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")
NB.sig.list.new <- list()
for (i in 1:33) {
  ct <- names(NB.sig.list)[i]
  ct_tb <- NB.sig.list[[ct]]
  ct_tb$Celltype <- ct
  ct_tb$Isoform_specific <- FALSE
  rest <- Reduce(rbind, NB.sig.list[-i])
  ct_iso <- ct_tb$phenotype_id
  rest_iso <- rest$phenotype_id
  ct_iso_sp <- setdiff(ct_iso, rest_iso)
  ct_tb$Isoform_specific[which(ct_tb$phenotype_id %in% ct_iso_sp)] <- TRUE
  NB.sig.list.new[[ct]] <- ct_tb
}
NB_top_isoQTL <- Reduce(rbind, NB.sig.list.new)

tmp <- distinct(NB_top_isoQTL, phenotype_id, .keep_all = TRUE)
rownames(tmp) <- tmp$phenotype_id
final_table <- final_table %>% group_by(phenotype_id) %>%
  mutate(Isoform_specific = tmp[phenotype_id,]$Isoform_specific)

write.table(final_table, file = "/data/lib14/project/scLongread/Tables/NB_sig_lungGWAS_sum.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
final_table <- read.table("/data/lib14/project/scLongread/Tables/NB_sig_lungGWAS_sum.txt",
                          sep = "\t", header = TRUE)
write.table(NB_top_isoQTL, file = "/data/lib14/project/scLongread/Tables/NB_eIsoform_sum.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
NB_top_isoQTL <- read.table("/data/lib14/project/scLongread/Tables/NB_eIsoform_sum.txt",
                            sep = "\t", header = TRUE)
Tensor_eIso_top <- read.table(file = "/data/lib14/project/scLongread/Tables/TensorQTL_eIsoform_sum.txt",
                              sep = "\t", header = TRUE)
length(unique(Tensor_eIso_top$phenotype_id.x))

Tensor_sig_nb_not <- data.frame(Transcript_id = setdiff(unique(Tensor_eIso_top$phenotype_id.x), unique(NB_top_isoQTL$phenotype_id)))
Tensor_sig_nb_not$Transcript_name <- Search_transcript_name2(Tensor_sig_nb_not$Transcript_id)
Tensor_sig_nb_not$Gene <- Search_transcript_gene_name(Tensor_sig_nb_not$Transcript_id)
write.table(Tensor_sig_nb_not, file = "/data/lib14/project/scLongread/Tables/Tensor_sig_nb_not_list.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
sr_Tensor_sig <- readRDS("/data/Choi_lung/TTL/tensor/Output_Sum_Final/full_sig_list.rds")
sr_Tensor_sig_long <- Reduce(rbind, sr_Tensor_sig)
length(unique(sr_Tensor_sig_long$phenotype_id))
NB_top_isoQTL <- NB_top_isoQTL %>% group_by(phenotype_id) %>% mutate(transcript_name = Search_transcript_name2(phenotype_id),
                                                                     gene_name = Search_transcript_gene_name(phenotype_id), 
                                                                     gene_id = Search_transcript_gene_id(phenotype_id), 
                                                                     .after = phenotype_id)
final_table <- final_table %>% group_by(phenotype_id) %>% mutate(gene_id = Search_transcript_gene_id(phenotype_id), 
                                                                .after = phenotype_id)

eGenes <- intersect(unique(sr_Tensor_sig_long$phenotype_id), unique(NB_top_isoQTL$gene_id))
setdiff(unique(final_table$gene_id), unique(sr_Tensor_sig_long$phenotype_id))
If_eGenes <- rep(TRUE, length(eGenes))
names(If_eGenes) <- eGenes
NB_top_isoQTL <- NB_top_isoQTL %>% group_by(gene_id) %>%
  mutate(If_eGenes = If_eGenes[gene_id])

NB_top_isoQTL <- NB_top_isoQTL %>% group_by(gene_name) %>%
  mutate(wrong_match = length(unique(gene_id)))

# compare with GTEx sGenes
GTEx_sGene <- read.table(gzfile("/data/Choi_lung/scLongreads/GTEx/GTEx_Analysis_v8_sQTL/Lung.v8.sgenes.txt.gz"),sep="\t", header = TRUE)
intersect(unique(GTEx_sGene$gene_name),NB_top_isoQTL$gene_name )
setdiff(NB_top_isoQTL$gene_name, unique(GTEx_sGene$gene_name) )
If_sGenes <- rep(TRUE, length(unique(GTEx_sGene$gene_name)))
names(If_sGenes) <- unique(GTEx_sGene$gene_name)
NB_top_isoQTL <- NB_top_isoQTL %>% group_by(gene_name) %>%
  mutate(If_sGenes = If_sGenes[gene_name])
write.table(NB_top_isoQTL, file = "/data/lib14/project/scLongread/Tables/NB_eIsoform_sum.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)



# check the distribution of ENST00000278282 (rs1434192708 )
bed <- read.table(gzfile("/data/Choi_lung/scLongreads/jaxqtl/Alveolar_macrophages_CCL3/phenotype.bed.gz"),
                  header = TRUE, comment.char = "", sep = "\t")
mtx <- apply((bed[,c(5:105)]), 2, as.numeric)
dim(mtx)
offset <- colSums(mtx)
hist(log(mtx[5942, ]/offset))
mu <- mean(as.numeric(bed[5942, c(5:105)]))
sigma <- var(as.numeric(bed[5942, c(5:105)]))
alpha = sigma/(mu - sigma)
p = mu/sigma
r = (mu^2)/(sigma - mu)
hist(rnbinom(101, r, p))
hist(rpois(101, mu))
pnorm(-0.4521367/0.02025496)*2
df <- data.frame(pois = rpois(5000, 1))

library(ggplot2)
library(ggh4x)

set.seed(42)
df <- data.frame(
  x = mtx[5942, ]
)

ggplot(df, aes(x)) +
  geom_bar(aes(y = after_stat(prop)),
           alpha = 0.5, width = 1, position = "identity") +
  stat_theodensity(distri = "nbinom")


library(snpStats)
fam <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.fam"
bim <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bim"
bed <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bed"
sample <- read.plink(bed, bim, fam)
mtx <- NormalizeData(mtx)
isoQTLs <- c("rs1434192708")
bed <- read.table(gzfile("/data/Choi_lung/scLongreads/jaxqtl/Alveolar_macrophages_CCL3/phenotype.bed.gz"),
                  header = TRUE, comment.char = "", sep = "\t")
genotype_mtx <- sample$genotypes[,isoQTLs]
snps_mtx <- genotype_mtx@.Data
snps_mtx <- snps_mtx[colnames(bed)[c(5:105)],]

isoQTL_sum <- data.frame(rs1434192708 = snps_mtx,
                         ENST00000278282 = mtx[5942,])
isoQTL_sum <- data.frame(rs1434192708 = snps_mtx,
                         ENST00000278282 = as.numeric(bed[which(bed$trascript_id == "ENST00000278282"),c(5:105)]))
which(isoQTL_sum$rs1434192708 == "01")
isoQTL_sum$rs1434192708_AAT <- "0|0"
isoQTL_sum$rs1434192708_AAT[which(isoQTL_sum$rs1434192708 == "01")] <- "1|1"
isoQTL_sum$rs1434192708_AAT[which(isoQTL_sum$rs1434192708 == "02")] <- "0|1"
table(isoQTL_sum$rs1434192708_AAT)
ggplot(isoQTL_sum, aes(x=rs1434192708_AAT, y=ENST00000278282, fill=rs1434192708_AAT)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="rs1434192708", y = "Count expr of ENST00000278282")+
  theme_classic()
tmp <- read.table("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/output/Alveolar_macrophages_CCL3/qtl_results_all.txt",
                  sep = "\t", header = TRUE)
