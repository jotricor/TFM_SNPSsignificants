
library(hexbin)
library(snpStats)
library(tidyverse)
library(ggrepel)


ped_assocs <- function(bed,bim,fam, MAF_filter, mind_filter, geno_filter, hwe_filter){
  
  ##### Read data in plink format #####
  
  data_patients <- read.plink(bed, bim, fam)
  
  ##### Summaries and graphics respect to patients and genotypes #####
  
  snpsum <- col.summary(data_patients$genotypes)
  cat("total number of patients input: ", nrow(data_patients$fam), "\n")
  cat("total number of SNPs input: ", nrow(snpsum), "\n")
  png("plot_prefiltered.png")
  par(mfrow= c(1,2))
  hist(snpsum$MAF)
  hist(snpsum$z.HWE)
  dev.off()
  sample.qc <- row.summary(data_patients$genotypes)
  png("plot_qc.png")
  par(mfrow=c(1,1))
  plot(sample.qc)
  dev.off()
  
  ##### Add --mind filter that delete the patients with down call rate #####
  
  min_call_sample <- 1-mind_filter
  use_min <- sample.qc$Call.rate > min_call_sample 
  data_patients$genotypes <- data_patients$genotypes[use_min, ]
  data_patients$fam <- data_patients$fam[use_min, ]
  cat("total number of patients after mind filter: ", nrow(data_patients$fam), "\n")
  
  
  ##### Delete sex chromosomes for subsequent genotype filtering ######
  
  data_patients$genotypes <- data_patients$genotypes[, data_patients$map$chromosome %in% 1:22]
  data_patients$map <- data_patients$map[data_patients$map$chromosome %in% 1:22, ]
  snpsum <- col.summary(data_patients$genotypes)
  cat("total number of snps after delete of chromosomes: ", nrow(snpsum), "\n")
  
  ##### Add --geno, --maf and --hwe filters #####
  
  gen_call_snp <- 1-geno_filter
  use_geno <- snpsum$Call.rate > gen_call_snp
  cat("total number of snps deleted with geno filter: ", sum(!snpsum$Call.rate>gen_call_snp), "\n")
  
  snpsum$MAF[is.na(snpsum$MAF)] <- 0
  use_maf <- snpsum$MAF >= MAF_filter
  cat("total number of snps deleted with maf filter: ", sum(snpsum$Call.rate>gen_call_snp & !snpsum$MAF>=MAF_filter), "\n")
  
  zvalueHWEmin <- qnorm(hwe_filter)
  zvalueHWEmax <- qnorm(1-hwe_filter)
  cat("total number of snps deleted with hwe filter: ", sum(snpsum$Call.rate>gen_call_snp & snpsum$MAF>=MAF_filter & !(snpsum$z.HWE>zvalueHWEmin & snpsum$z.HWE<zvalueHWEmax)), "\n")
  
  use_genotypes <- snpsum$MAF >= MAF_filter & 
    snpsum$Call.rate> gen_call_snp &
    (use_hwe <- snpsum$z.HWE>zvalueHWEmin & snpsum$z.HWE<zvalueHWEmax)
  use_genotypes[is.na(use_genotypes)] <- FALSE
  data_patients$genotypes <- data_patients$genotypes[,use_genotypes]
  snpsum <- col.summary(data_patients$genotypes)
  cat("total number of SNPs output: ", nrow(snpsum), "\n")
  
  ##### Statistical tests #####
  
  tests <- single.snp.tests(affected, data = data_patients$fam, snp.data = data_patients$genotypes)
  
  
  chi2 <- chi.squared(tests, df=1)
  png("plot_chi2.png")
  par(mfrow=c(1,1))
  qq.chisq(chi2, df = 1)
  dev.off()
  
  position <- data_patients$map[use_genotypes, "position"] 
  chr<- data_patients$map[use_genotypes, "chromosome"]
  name_snp <-  data_patients$map[use_genotypes, "snp.name"]
  p1 <- p.value(tests, df = 1)
  dataf1 <- data.frame(SNP=name_snp,CHR=chr,BP=position,P=p1) 
  dataf1[is.na(dataf1)] <- 0 #reemplazar na por false.
  dataf1 <- dataf1[dataf1$CHR %in% 1:22,]
  
  
  dataf1$BP <- as.double(dataf1$BP)
  return(dataf1)
}


dataf1 <- ped_assocs(bed ="test_binary.bed", bim ="test_binary.bim", fam="test_binary.fam", MAF_filter = 0.01, mind_filter = 0.1, geno_filter = 0.1, hwe_filter = 0.01)


manhatan_plot <- function(dataf1){
  snps_int <- dataf1$SNP[dataf1$P < 5E-05]
  don <- dataf1 %>% 
    
    #### Calculate chromosome len #####
  group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    ##### Calculate cumulative position of each chromosome #####
  mutate(tot=cumsum(chr_len)-chr_len) %>% #Para que no se solapen los SNPs en una única columna
    select(-chr_len) %>%
    
    ##### Add this information to original data #####
  left_join(dataf1, ., by=c("CHR"="CHR")) %>%
    
    ##### Add a cumulative position of each SNP #####
  arrange(CHR, BP) %>% 
    mutate( BPcum=BP+tot) %>% 
    
    ##### Add highlight and annotation information #####
  mutate( is_highlight=ifelse(SNP %in% snps_int, "yes", "no"))  %>% 
    mutate( is_annotate=ifelse(-log10(P)>-log10(5E-05), "yes", "no"))
  ##### Prepare X axis #####
  dataf1$BPcum <- as.double(don$BPcum)
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  
  if(max(-log10(dataf1$P)) > 8) {
    max_y_lim <- max(-log10(gwasResults$P)) + 0.5 
  } else {
    max_y_lim <- 8.5
  }
  
  mp <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
    
    ##### Show all points #####
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    
    ##### Custom X axis #####:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    
    ##### NEW- to set the y-axis limit #####
  ylim(0,max_y_lim) +
    # Add highlighted points
    geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    
    ##### NEW - ADD LINES #####
  geom_hline(yintercept=-log10(5E-08), linetype="dashed", color = "red") +
    geom_hline(yintercept=-log10(5E-05), linetype="dotted", color = "blue") +
    
    ##### Add label using ggrepel to avoid overlapping #####
  geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
    xlab("Chromosome") + 
    ##### Custom the theme #####
  theme_bw() +
    theme( 
      legend.position="none",
      axis.text.x = element_text(angle = 90, size = 6),
      axis.text.y = element_text(size=6),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  tiff(paste0("manhattan_smaller_", ".tiff"), compression="lzw", width=1500, height=850, res=250); plot(mp); dev.off()
}
