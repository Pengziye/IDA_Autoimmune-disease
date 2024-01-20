###Install and load the required packages
library(TwoSampleMR)
library(ggplot2)
library(dplyr)
library(MRPRESSO)
library(TwoSampleMR)
library(LDlinkR)
###Read file and Select significant p values
setwd("your working path")
exposure<-read.csv("CRC.csv",header = T)
exposure_dat<-subset(exposure,exposure$pvalue <5e-08)  
write.csv(exposure_dat,file="exposure_sign_p.csv",quote=FALSE)

#Read exposure data and perform clumping with standard parameters
exposure_after <- system.file("extdata/exposure_sign_p.csv",package = "TwoSampleMR")
#Note that the exposure_after.csv file needs to be placed in the extdata folder under the R installation path.
exposure_dat <- read_exposure_data(
  filename = exposure_after,
  clump = TRUE,
  sep= ",",           
  snp_col = "rsID",   
  beta_col = "effect",   
  se_col = "se",    
  effect_allele_col ="EA",
  other_allele_col = "NEA",
  eaf_col = "EAF",
  pval_col = "pvalue",
  samplesize_col= "SampleSize",
  ncase_col= "N_cases", 
  ncontrol_col= "N_controls"
) 


###Extract outcome SNPs

outcome = fread("outcome.txt", sep = "\t")
outcome_dat <- read_exposure_data(
  snps = exposure_dat$SNP,
  filename = outcome,
  sep= ",",           
  snp_col = "rsID",   
  beta_col = "effect",   
  se_col = "se",    
  effect_allele_col ="EA",
  other_allele_col = "NEA",
  eaf_col = "EAF",
  pval_col = "pvalue",
  samplesize_col= "SampleSize",
  ncase_col= "N_cases", 
  ncontrol_col= "N_controls"
) 

###Harmonization of data
dat<-harmonise_data(exposure_dat =exposure_dat,outcome_dat = outcome_dat )

#Calculate R2 per SNP
dat$rq=2*(dat$beta.exposure)^2*(dat$eaf.exposure)*(1-dat$eaf.exposure)

#Calculate R2 total
weight_mean_R2_PSC<-sum(dat$rq * (1/dat$se.exposure)/sum(1/dat$se.exposure))

#Calculate a per-SNP F statistic 
dat$EAF2 <- (1 - dat$eaf.exposure)
dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2) 

PVEfx <- function(BETA, MAF, SE, N){      
  pve <- (2*(BETA^2)*MAF*(1 - MAF))/
    ((2*(BETA^2)*MAF*(1 - MAF)) + ((SE^2)*2*N*MAF*(1 - MAF))) 
  return(pve) 
} 
N <- sample_size   
dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N ) 
dat$FSTAT <- ((N - 1 - 1)/1)*(dat$PVE/(1 - dat$PVE)) 

# Calculate a total instrument F statistic 
k <- nrow(subset(dat, dat$ambiguous == FALSE))
F <- ((N - 1 - k)/k)*((sum(dat$PVE))/(1 - sum(dat$PVE))) 

# Calculate I2
res_dat<-mr_singlesnp(dat,all_method = c("mr_ivw")) 
res_dat2<-res_dat[grep("^rs",res_dat$SNP),]  
rea_meta_dat<-metafor::rma(yi=res_dat2$b,
                           sei = res_dat2$se,
                           weights = 1/ dat$se.outcome^2,
                           data = res_dat2,
                           method = "FE")


#mr_heterogeneity to detect heterogeneity test
heterogeneity_Anaemia_PSC<-mr_heterogeneity(dat)

#MR PRESSO to detect horizontal pleiotropy
mr_presso_PSC<-run_mr_presso(dat,NbDistribution = 3000) 

mr_intercept_PSC<-mr_pleiotropy_test(dat)


#Perform MR
out_dat<-generate_odds_ratios(mr_res = mr(dat)) 
#Leave-one-out plot
dat_leave<-mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))


#Multivariable MR
library(ggplot2)
library(dplyr)
library(MRPRESSO)
library(TwoSampleMR)
library(readr) 
library(MRInstruments)

#Read exposures and outcome
id_exposure <- c("exposure1.csv", "exposure2.csv")
id_outcome <- "outcome.csv"

#Extract and clump significant IVs for 2 different exposures
 
exposure_dat <- mv_extract_exposures_local(
  id_exposure,
  sep= ",",           
  snp_col = "rsID",   
  beta_col = "effect",   
  se_col = "se",    
  effect_allele_col ="EA",
  other_allele_col = "NEA",
  eaf_col = "EAF",
  pval_col = "pvalue",
  samplesize_col= "SampleSize",
  ncase_col= "N_cases", 
  ncontrol_col= "N_controls", 
  log_pval =FALSE, min_pval =1e-200, 
  pval_threshold =5e-08, 
  clump_r2 =0.001, 
  clump_kb =10000, 
  harmonise_strictness =2)

outcome_dat<-read_outcome_data(
  snps =exposure_dat$SNP, 
  filename ="outcome.txt", 
  sep =",", 
  snp_col = "rsID",   
  beta_col = "effect",   
  se_col = "se",    
  effect_allele_col ="EA",
  other_allele_col = "NEA",
  eaf_col = "EAF",
  pval_col = "pvalue",
  samplesize_col= "SampleSize",)

mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

res <- mv_multiple(mvdat)



