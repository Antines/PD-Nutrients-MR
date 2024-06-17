## Mendelian randomization between nutrients and PD

remotes::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("MRCIEU/MRInstruments")
devtools::install_github("rondolab/MR-PRESSO")
library(MRInstruments)
library(TwoSampleMR)
require(ggplot2)	
library(dplyr)
library(MRPRESSO)
library(data.table)
library(httr)

setwd("<work/dictionary>")
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')  #clumping's html options required for function to work

######################################################################################################################
#                               Search GWAS datasets by phenotype title
######################################################################################################################

not_found <- c() #vector for undiscovered phenotypes 
prepare_GWAS_datasets <- function(phenotypes, catalog_data) {
  for (phenotype in phenotypes) {
    pattern = paste0("\\b(", paste(phenotype, collapse = "|"), ")\\b") #search by regex
    
    # Subset data by phenotype and ancestry
    subset_data <- catalog_data %>%
      filter(grepl(pattern, Phenotype_simple, ignore.case = T) | 
               grepl(pattern, MAPPED_TRAIT_EFO, ignore.case = T) |
               grepl(pattern, Phenotype, ignore.case = T))%>%
      filter(pval <= 5e-08 & grepl("(?i)European|Irish|Caucasian", Initial_sample_description))
    
    # Check data is available
    if (nrow(subset_data) == 0) {
      warning(paste("No data available for phenotype:", phenotype))
      not_found <<- c(not_found, phenotype)
      next  
    }
    
    # New column with counted sample size (init + replication) -> filter the biggest
    subset_data <- subset_data %>%
      mutate(Sample_Size = sapply(stringr::str_extract_all(Initial_sample_description, "\\b\\d{1,3},?\\d{1,3}\\b"), 
                                  function(x) sum(as.numeric(gsub(",", "", x))))) %>%
      mutate(Repl_sample_Size = if (!all(is.na(Replication_sample_description))) {
        sapply(stringr::str_extract_all(Replication_sample_description, "\\b\\d{1,3},?\\d{1,3}\\b"), 
               function(x) sum(as.numeric(gsub(",", "", x))))
      } else {
        NA
      }) %>%
      mutate(Total_Sample_Size = Sample_Size + ifelse(is.na(Repl_sample_Size), 0, Repl_sample_Size)) %>%
      filter(Total_Sample_Size == max(Total_Sample_Size))
    
    # Format data for MR
    formatted_data <- format_data(subset_data,  type = "exposure", snps = NULL, snp_col = "SNP", beta_col = "beta", se_col = "se",  eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "pval", log_pval = FALSE)
    formatted_data$samplesize.exposure <- subset_data$Sample_Size[1]
    
    # Handle SE calculation if not present
    if (!"se" %in% colnames(formatted_data)) {
      formatted_data$se.exposure <- TwoSampleMR:::get_se(formatted_data$beta.exposure, formatted_data$pval.exposure)
    }
    
    #Clump with defaults parameters
    formatted_data <- clump_data(formatted_data, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, pop = "EUR") #set options() if error occured
    
    # Write to CSV
    file_name <- paste0(gsub(" ", "_", phenotype), "_GWAS.csv")
    write.csv(formatted_data, file = file_name, quote = FALSE)
    
    # Create variables in the global environment
    assign(paste0(gsub(" ", "_", phenotype), "_init"), subset_data, envir = .GlobalEnv)
    assign(paste0(gsub(" ", "_", phenotype), "_dat"), formatted_data, envir = .GlobalEnv)
  }
}


phenotypes <- c("dairy", "phosphate", "calcium", "milk", "vitamin a", "thiamine", "vitamin b2 riboflavin", 
                "vitamin b3 niacin", "vitamin b5 pantothenic acid", "vitamin b6 pyridoxine", "vitamin b7 biotin", 
                "vitamin b9 folate", "vitamin b12 cobalamin", "vitamin c ascorbic acid", "vitamin k", 
                "folic acid", "calciferol", "phylloquinone", "alpha-tocopherol", "vitamin e", 
                "vitamin d", "thiamin", "fluoride", "mineral", "copper", "iodine", "iron", "magnesium", "manganese", 
                "molybdenum", "chromium", "phosphorus", "potassium", "selenium", "sodium", "sulfur", "zinc",
                "fe level", "Mg", "protein intake", "fat intake", "carbohydrates intake")
prepare_GWAS_datasets(phenotypes, gwas_catalog)

not_found


######################################################################################################################
#                               MR; PD as outcome, nutrients as exposure
######################################################################################################################

#Create variable with file and dictionary names
file <- list.files()[grepl("\\.csv$", list.files())]
folder <- tools::file_path_sans_ext(file)
for (chr in folder) {
  dir.create(chr)
}

#Run MR
for(i in 1:length(file)){
  
  exp_data <- read_exposure_data(file[i], sep=",", snp_col = "SNP", beta_col = "beta.exposure",  eaf_col = "eaf.exposure", se_col="se.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col= "other_allele.exposure", pval_col = "pval.exposure", samplesize_col = "samplesize.exposure", clump=FALSE)
  if(all(is.na(exp_data$effect_allele.exposure)) && all(is.na(exp_data$other_allele.exposure)) || all(is.na(exp_data$eaf.exposure))) {
    next
  }

  out_data <- read_outcome_data(snps = exp_data$SNP,filename = "PD_GWAS_noUKBB.txt", sep="\t",  snp_col = "SNP", beta_col = "Effect", se_col = "StdErr", eaf_col = "Freq1", effect_allele_col = "Allele1", other_allele_col = "Allele2", pval_col = "P-value", ncase_col = "ncase", ncontrol_col = "ncontrol",samplesize_col="samplesize.outcome")
  out_data$r.outcome <- get_r_from_lor(out_data$beta.outcome, out_data$eaf.outcome, out_data$ncase.outcome, out_data$ncontrol.outcome, 0.01,  model = "logit")
  dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
  dat$units.outcome<-"log odds" #since PD dataset contain binary trait
  #dat$units.exposure<-"log odds" #all exposure are continious traits
  
  dat1<-subset(dat, dat$eaf.exposure!="NA")
  dat1$r.exposure <- get_r_from_pn(dat1$pval.exposure, dat1$samplesize.exposure)
  #dat1$r.exposure<- get_r_from_lor(dat1$beta.exposure, dat1$eaf.exposure, dat1$ncase.exposure, dat1$ncontrol.exposure, 0.01,  model = "logit")
  
  steiger <- steiger_filtering(dat1)
  sig<-subset(steiger, steiger$steiger_dir==TRUE)
  
  if (nrow(sig) == 0) {next}
  
  tryCatch({
    presso <- mr_presso(BetaOutcome = "beta.outcome", 
                        BetaExposure = "beta.exposure", 
                        SdOutcome = "se.outcome", 
                        SdExposure = "se.exposure", 
                        OUTLIERtest = TRUE, 
                        DISTORTIONtest = TRUE, 
                        data = sig, 
                        NbDistribution = 1000, 
                        SignifThreshold = 0.05)
    
    capture.output(print(presso), file = paste(folder[i], "presso.txt", sep = "/"))
  }, error = function(e) {
    write("Not enough IV occurance", file = paste(folder[i], "Fail_to_presso.txt", sep = "/"))
  })
  mr(sig)
  
  #F-statitistics and R2 for each exposure
  R2<-mean(sig$r.exposure)
  capture.output(print(R2), file = paste(folder[i], "r2.txt", sep = "/"))
  n <- mean(sig$samplesize.exposure)
  k <- nrow(subset(sig, sig$ambiguous == FALSE))
  F<-(R2*(n-1-k))/((1-R2)*k)
  capture.output(print(F), file = paste(folder[i], "f.txt", sep = "/"))
  
  #Output tables and graphs summarizing the results
  mr_report(sig, study = folder[i],output_path = folder[i]) 
  res_single <- mr_singlesnp(sig)
  p5 <- mr_forest_plot(res_single)
  p5[[1]]
  ggsave(p5[[1]], file= "plot.jpg", path=folder[i] , width=7, height=12)
}


######################################################################################################################
                                   #Reverse MR; PD as exposure, nutrients as outcome
######################################################################################################################

setwd("<work/dictionary>")
file <- list.files()[grepl("\\.out.tsv$", list.files())]
folder <- tools::file_path_sans_ext(file)
for (chr in folder) {
  dir.create(chr)
}

# Reverse MR; PD as exposure, nutrients as outcome
data <- data.frame(fread("../PD_GWAS_noUKBB.txt", sep = "\t"))
PD_GWAS <- filter(data, P.value < 5e-8 & P.value != 0)
PD_GWAS <- format_data(
  PD_GWAS,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  eaf_col = "Freq1",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P-value",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize"
)

write.csv(PD_GWAS, file = "PD_as_exposure", quote = FALSE)
pd_exp_data <- read_exposure_data("PD_as_exposure", sep=",", snp_col = "SNP", beta_col = "beta.exposure", 
                                  eaf_col = "eaf.exposure", se_col="se.exposure", 
                                  effect_allele_col = "effect_allele.exposure", 
                                  other_allele_col= "other_allele.exposure", pval_col = "pval.exposure",
                                  samplesize_col = "samplesize.exposure", ncase_col = "ncase.exposure", 
                                  ncontrol_col = "ncontrol.exposure", log_pval = FALSE)

pd_clump <- clump_data(pd_exp_data, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, pop = "EUR")

for(i in 1:length(file)){
  exp_data <- pd_clump 
  out_data <- read_outcome_data(snps = exp_data$SNP, filename = file[i], sep="\t", 
                                snp_col = "SNP", beta_col = "beta",  eaf_col = "EAF",
                                se_col="se", effect_allele_col = "A1", 
                                other_allele_col= "A2", pval_col = "Pvalue",
                                samplesize_col = "N")

  out_data$r.outcome<- get_r_from_pn(out_data$pval.outcome, out_data$samplesize.outcome)
  dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
  #dat$units.outcome<-"log odds"
  dat$units.exposure<-"log odds"
  dat1 <- subset(dat, dat$eaf.exposure!="NA")
  dat1$r.exposure<- get_r_from_lor(dat1$beta.exposure, dat1$eaf.exposure, dat1$ncase.exposure, dat1$ncontrol.exposure, 0.01,  model = "logit")
  
  steiger <- steiger_filtering(dat1)
  sig<-subset(steiger, steiger$steiger_dir==TRUE)
  
  if (nrow(sig) == 0) {next}
  
  tryCatch({
    presso <- mr_presso(BetaOutcome = "beta.outcome", 
                        BetaExposure = "beta.exposure", 
                        SdOutcome = "se.outcome", 
                        SdExposure = "se.exposure", 
                        OUTLIERtest = TRUE, 
                        DISTORTIONtest = TRUE, 
                        data = sig, 
                        NbDistribution = 1000, 
                        SignifThreshold = 0.05)
    
    capture.output(print(presso), file = paste(folder[i], "presso_rev.txt", sep = "/"))
  }, error = function(e) {
    write("Not enough IV occurance", file = paste(folder[i], "Fail_to_presso_rev.txt", sep = "/"))
  })
  
  mr(sig)
  mr_report(sig, study = folder[i], output_path = folder[i], output_type = "md") 
  res_single <- mr_singlesnp(dat)
  p5 <- mr_forest_plot(res_single)
  p5[[1]]
  ggsave(p5[[1]], file= "plot_rev.jpg", path=folder[i] , width=7, height=12)
}


