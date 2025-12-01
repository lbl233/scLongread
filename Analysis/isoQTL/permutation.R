# Permutation of DEI

# 


library(stringr)
library(stringi)
library(dplyr)
library(ggplot2)
library(Seurat)

files <- list.files("/data/Choi_lung/scLongreads/DEI/WholegroupPhaseI/output_permutation/UnfilteredResults/", pattern = ".csv")


celltypes <- str_split_fixed(files, "_", 5)[,5]
celltypes <- str_split_fixed(celltypes, "\\.", 2)[,1]
celltypes <- gsub(" ", "_", celltypes)

DEI_per.list <- list()
Sig_DEI <- NULL
for (i in 1:37) {
  DEI_rst <- read.csv(paste0("/data/Choi_lung/scLongreads/DEI/WholegroupPhaseI/output_permutation/UnfilteredResults/",files[i]))
  DEI_per.list[[celltypes[i]]] <- DEI_rst
  Sig_DEI <- c(length(which(DEI_rst$padj < 0.05 & abs(DEI_rst$log2FoldChange) > 1)), Sig_DEI)
}
names(DEI_per.list) <- celltypes

library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")

qqplotdata_permutation <- function(logpvector, logpermutatedp){
  o = sort(logpvector,decreasing=T)
  e = sort(logpermutatedp,decreasing=T)
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


output_path <- "/data/Choi_lung/scLongreads/DEI/WholegroupPhaseI/output_permutation/plot/"

for (celltype in celltypes) {
  
  message("Celltype: ", celltype)
  
  all_DEI <- DEI.list[[celltype]]
  perm_DEI <- DEI_per.list[[celltype]]
  
  pvals <- all_DEI$padj
  pvals_permutated <- perm_DEI$pvalue
  idx <- union(which(is.na(pvals)), which(is.na(pvals_permutated)))
  print(length(idx))
  pvals <- pvals[!is.na(pvals)]
  legendcol <- character(0)
  ycol <- "log10P"
  plotdata <- qqplotdata_permutation(-log10(pvals), -log10(pvals_permutated))
  plotdata <- plotdata[!(plotdata$o==Inf),]
  fx <- plotdata$e
  fy <- plotdata$o
  fcol <- rep("#E41A1C",length(plotdata$o))
  opt =  list(break.top = ceiling(max(fx)),
              top.size = 0.125)
  png(filename = paste0(output_path,"QQplot_permutation_vs_oberved_", celltype, ".png"), width = 8, height = 8, units = "in",res=300)
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
       xlab=expression(plain(Permutation)~~group("(",-log[10]*italic(P),")")),
       ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
       cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
       main=opt$maintitle,pch=19)
  
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
  title(paste0("The expected p-values under the permutation (x-axis) \n against the observed p-values (y-axis) in ",celltype))
  dev.off()
  # null vs observed p
  plotdata <- qqplotdata(-log10(pvals))
  plotdata <- plotdata[!(plotdata$o==Inf),]
  fx <- plotdata$e
  fy <- plotdata$o
  fcol <- rep("#E41A1C",length(plotdata$o))
  opt =  list(break.top =  ceiling(max(fx)),
              top.size = 0.125)
  png(filename = paste0(output_path,"QQplot_null_vs_oberved_", celltype, ".png"), width = 8, height = 8, units = "in",res=300)
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
       xlab=expression(plain(Null)~~group("(",-log[10]*italic(P),")")),
       ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
       cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
       main=opt$maintitle,pch=19)
  
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
  title(paste0("The expected p-values under the null hypothesis (x-axis) \n against the observed p-values (y-axis) in ",celltype))
  dev.off()
  
  # null vs permutation
  plotdata <- qqplotdata(-log10(pvals_permutated))
  plotdata <- plotdata[!(plotdata$o==Inf),]
  fx <- plotdata$e
  fy <- plotdata$o
  fcol <- rep("#E41A1C",length(plotdata$o))
  opt =  list(break.top = ceiling(max(fx,na.rm=T)),
              top.size = 0.125)
  png(filename = paste0(output_path,"QQplot_null_vs_permutation_", celltype, ".png"), width = 8, height = 8, units = "in",res=300)
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
       xlab=expression(plain(Null)~~group("(",-log[10]*italic(P),")")),
       ylab=expression(plain(Permutated)~~group("(",-log[10]*italic(P),")")),
       cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
       main=opt$maintitle,pch=19)
  
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
  title(paste0("The expected p-values under the null hypothesis (x-axis) \n against the permutated p-values (y-axis) in ",celltype))
  dev.off()
}

# FDR
df.list <- list()
for (celltype in celltypes) {
  rst <- DEI_per.list[[celltype]]
  df <- data.frame(x = as.numeric(c(0, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2)),
                   sig_pct = c(0, sum((rst$padj < 0.00001), na.rm = TRUE),
                               sum((rst$padj < 0.0001), na.rm = TRUE),
                               sum((rst$padj < 0.001), na.rm = TRUE),
                               sum((rst$padj < 0.01), na.rm = TRUE),
                               sum((rst$padj < 0.05), na.rm = TRUE),
                               sum((rst$padj < 0.1), na.rm = TRUE),
                               sum((rst$padj < 0.2), na.rm = TRUE)))
  #df$x <- factor(df$x, levels = as.numeric(c(0, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2)))
  df$sig_pct <- (df$sig_pct/nrow(rst))
  df$Celltype <- celltype
  df.list[[celltype]] <- df
}
df <- do.call(rbind, df.list)

df %>% group_by(x) %>% summarise(num_off = sum(sig_pct >= 0.05))

df$Celltype <- factor(df$Celltype, levels = gsub(" ", "_", levels(lr$Celltype)))
ggplot(df, aes(x = x, y = sig_pct, fill = Celltype)) + 
  geom_bar(stat="identity",position = "identity") +
  scale_fill_manual(values = cellcolors) + 
  geom_hline(yintercept = 5) + 
  scale_x_discrete(labels = paste0(round(c(0, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2)*100, digits=3), "%")) +
  xlab("FDR cutoff") + ylab("Percentage of significant isoform (%)") +
  facet_wrap(~Celltype, ncol = 8)+
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 0, face = "plain"),
        strip.text = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")
ggsave("/data/lib14/project/scLongread/FigSMethod.pdf", width = 20,height = 15)

