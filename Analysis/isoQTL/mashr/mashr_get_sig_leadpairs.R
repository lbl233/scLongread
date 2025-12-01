# Combine the results of multiple chunks
# Bolun Li 
# Jan 27th 2025


library(rhdf5)
library(data.table)
library(tidyverse)
library(argparse)
library(collapse)
# BiocManager::install("multistateQTL", force = TRUE,)
setwd("/data/Choi_lung/scLongreads/mashr/output")
res_list <- list.files("/data/Choi_lung/scLongreads/mashr/output/", recursive = TRUE, full.names = TRUE)
# res_list <- res_list[grepl("_200K.rds", res_list)]
res_list <- res_list[grepl("leadpairs", res_list)]
# res_list <- res_list[!grepl(".txt", res_list)]

for(r in 1:length(res_list)){
  tmp <- readRDS(res_list[[r]])
  # tmp <- mash(tmp, g=get_fitted_g(mashFit), fixg=TRUE) # run mash
  # saveRDS(tmp, res_list[[r]])
  tmp <- QTLExperiment::mash2qtle(tmp, sep="\\|")
  if(r == 1){
    msqe <- tmp
    message("Finished ", r, "/1000")
  } else{
    msqe <- QTLExperiment::rbind(msqe, tmp)
    message("Finished ", r, "/1000")
  }
}
saveRDS(msqe, "mashr_1M_leadpairs.rds")
msqe <- multistateQTL::callSignificance(msqe, assay="lfsrs", thresh=0.05)
msqe <- multistateQTL::getSignificant(msqe, assay="lfsrs")
message("# sig in : ", nrow(msqe))
saveRDS(msqe, "mashr_1M_leadpairs_sig.rds")
