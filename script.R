source("functions.R") 

#### Use the first function of functions file to obtain the dataframe with most significant SNPs ####
dataf1 <- ped_assocs(bed ="test_binary.bed", bim ="test_binary.bim", fam="test_binary.fam", MAF_filter = 0.01, mind_filter = 0.1, geno_filter = 0.1, hwe_filter = 0.01)

#### Use the second function of functions file to obtain the most significant SNPs image in tiff format, but you can use other format and resolution ####
tiff(paste0("manhattan",".tiff"), compression="lzw", width=1500, height=850, res=250)
  manhatan_plot(dataf1)
dev.off()
