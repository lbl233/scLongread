library(vroom)
library(dplyr)
library(data.table)
library(collapse)
library(rhdf5)
library(mashr)
library(ashr)
#library(tidyverse)
source("/data/lib14/R/Rscripts/utilities.R")

mashData <- fun_h5_2_mashr("/data/Choi_lung/scLongreads/mashr/lead_isoQTLs_ez.h5", max.missing=0)

data <- mash_set_data(mashData$Bhat,
                      mashData$Shat,
                      V=mashData$V)
mashFit <- readRDS("mashr_all_1M_fit_ez.rds")
setwd("/data/Choi_lung/scLongreads/mashr/output/")
mashData_chunk <- mash(data, g=get_fitted_g(mashFit), fixg=TRUE, algorithm.version = 'R',
                       posterior_samples = 100)

# magnitude
msqe <- multistateQTL::callSignificance(mashData_chunk, assay="lfsrs", thresh=0.05)
beta_mtx <- msqe@assays@data$betas
dim(beta_mtx)
error <- h5read("lead_isoQTL.h5", "error")
row.names(error) <- rownames
colnames(error) <- colnames
betas <- h5read("lead_isoQTLs_ez.h5", "betas")
dim(betas)
error[betas==0] <- 0

# calculate beta using initial se
beta_final <- beta_mtx*error
lpairs_pm <- beta_final
rst.list <- list()
for (j in 1:33) {
  tmp1 <- rst_magnitude[which(rowSums(rst_magnitude) == j),]
  tmp2 <- lpairs_pm[which(rowSums(rst_magnitude) == j),]
  isoQTLin2.list <- list()
  for (i in 1:nrow(tmp1)) {
    idx <- which(tmp1[i,])
    isoQTLin2.list[[i]] <- tmp2[i,idx]
  }
  isoQTLin2 <- Reduce(rbind, isoQTLin2.list)
  rst <- rowSums(isoQTLin2 > 0)
  rst.list[[j]] <- rst
  if(length(which(!(rst %in% c(0,j)))) != 0){
    message(j)
  }
}

top_beta <- c()
beta_range <- c()
test <- apply(lpairs_pm, 1, function(x){
  min <- min(x)
  max <- max(x)
  beta_range <- abs(abs(max) - abs(min))
  # if all betas with the same sign
  if(min*max > 0){
    if(max > 0){
      x = x/max
      top_beta <- c(top_beta,max)
    }else{
      x = x/min
      top_beta <- c(top_beta,min)
    }
  }else{
    # if betas have opposite direction
    # choose the maximum of absolute value as the top beta for the pair
    if((max+min) > 0){
      x = x/max
      top_beta <- c(top_beta,max)
    }else{
      x = x/min
      top_beta <- c(top_beta,min)
    }
  }
  return(list((x > 0.5 & x <= 1),top_beta, beta_range))
})
rst_magnitude <- NULL
for (i in 1:length(test)) {
  rst_magnitude <- rbind(rst_magnitude, test[[i]][[1]])
}

# summarize isoQTL sharing across cell types
rownames(rst_magnitude) <- names(test)
magnitude_sum <- rowSums(rst_magnitude)
rst_sum_pairs <- as.data.frame(t(rst_magnitude))
rst_sum_pairs$Celltype <- rownames(rst_sum_pairs)
rst_sum_pairs$Celltype <- factor(rst_sum_pairs$Celltype, levels = names(isoQTL_sig_list))
rst_sum_pairs <- rst_sum_pairs[order(rst_sum_pairs$Celltype),]
rst_sum_pairs$Categroy <- c(rep("EPI",8),rep("IMMUNE",15), rep("EC", 6), rep("STROMA",4))
df <- data.frame(pairs = names(magnitude_sum),
                 shared_num = magnitude_sum)
df_sum <- as.data.frame(table(df$shared_num))
df_sum$aspect <- "magnitude"
# figure 4c
ggplot(df_sum, aes(x=Var1, y=Freq, fill = "#fa8775")) +
  geom_bar(stat="identity") + xlab("Number of shared cell type")+ylab("number of isoQTLs") +
  geom_text(aes(label = Freq),vjust = -.5)+theme_classic()

# Sign
top_beta <- NULL
for (i in 1:length(test)) {
  top_beta <- c(top_beta, test[[i]][[2]])
}
lpairs_sign <- apply(lpairs_pm, 2, function(x) (x*top_beta)>0)

sign_sum <- rowSums(lpairs_sign)
df1 <- data.frame(pairs = names(sign_sum),
                  shared_num = sign_sum)
df_sum1 <- as.data.frame(table(df1$shared_num))
df_sum1$aspect <- "sign"

ggplot(df_sum1, aes(x=Var1, y=Freq, fill = "#0000ff")) +
  geom_bar(stat="identity") + xlab("Number of shared cell type")+ylab("number of isoQTLs") +
  geom_text(aes(label = Freq),vjust = -.5)+theme_classic()+RotatedAxis()+NoLegend()

df_final <- rbind(df_sum, df_sum1)
df_final$Freq <- round((df_final$Freq/2842)*100, 1)
ggplot(df_final, 
       aes(x=Var1, y=Freq, fill = aspect)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("Number of shared cell type")+ylab("Number of lead isoQTL-eIsoform pairs") +
  geom_text(aes(label = Freq), , vjust=-.1, color="black",
            position = position_dodge(0.9), size=3)+
  theme(axis.text = element_text(size = 14, colour = "black"))+
  theme_classic()
