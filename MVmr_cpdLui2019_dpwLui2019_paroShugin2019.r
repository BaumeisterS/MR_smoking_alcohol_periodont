#+ echo = TRUE, warning = FALSE, message = FALSE
#' --------- 
#' title: "MR Cigarettes per day (Lui2019) and Peridontitis  (SHRINE2019) - adjusted for Drinks per week (Liu 2019)
#' Author: "SE Baumeister"
#' Version: "16.02.2021"
#' output: html_document
#' date: "`r Sys.Date()`"
-------------------------------
Sys.setenv(LANG = "en")
.libPaths("c:/R/R-4.0.3/library/")
.libPaths()
#update.packages(ask=FALSE) 
#remove.packages("tibble")
x<-c("xfun", "digest","rlang","tidyverse","devtools","scales","data.table","survey","knitr",
     "flextable","officer","ggpubr","EValue","meta","survival", "acepack","R.utils","LDlinkR",
     "MendelianRandomization","ashr")
#install.packages(x)
lapply(x, require, character.only = TRUE)
(.packages())
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("biomaRt")
# devtools::install_github("jdstorey/qvalue")
# devtools::install_github("slowkow/proxysnps")
# BiocManager::install("myvariant")
# devtools::install_github("MRCIEU/TwoSampleMR", upgrade ="always", force = TRUE)
# install.packages("http://cnsgenomics.com/software/gsmr/static/gsmr_1.0.9.tar.gz",repos=NULL,type="source")
# #devtools::install_github("jean997/cause")
#  install.packages("remotes")
# remotes::install_github("gqi/MRMix", force=TRUE)
#devtools::install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
y<-c("biomaRt","qvalue","proxysnps","myvariant","TwoSampleMR","MRInstruments","MRPRESSO","gsmr","MRMix", "MVMR","cause")
lapply(y, require, character.only = TRUE)
(.packages())
packageVersion("TwoSampleMR")
sessionInfo()
# Load MR functions
source("c:/R_fct/mr_fcts.R")
ls()
setwd("c:/arbeit/Publikationen/mr_smoking_alc_paro/data/MV_cpdLui2019noUKBB_dpwLui2019noUKBB_paroShugin2019/")
getwd()
mydir <- getwd()
#--------------------
# Exposure data: Cigarettes per day, Lui 2019, PMID  30643251, 
#--------------------
cpd<-fread("c:/gwas_summary/smoking/Liu2019/CigarettesPerDay.WithoutUKB.txt.gz") %>%
  format_data(. , type="outcome", chr_col="CHROM", phenotype_col="exposure", snp_col="RSID",
              beta_col="BETA", se_col="SE", eaf_col="AF",
              effect_allele_col="ALT",
              other_allele_col="REF", pval="PVALUE") %>% 
  mutate(outcome="cpd", id.outcome="cpd",Phenotype="cpd") %>% 
  dplyr::rename(beta=beta.outcome,se=se.outcome,
                effect_allele=effect_allele.outcome,
                other_allele=other_allele.outcome,
                eaf=eaf.outcome,
                pval=pval.outcome) %>% 
  dplyr::select(SNP,outcome,Phenotype,id.outcome,beta,se,pval,effect_allele,other_allele,eaf)
#--------------------
# 1.0) Exposure data: Drinks per week, Lui 2019, PMID  30643251, 
#--------------------
dpw<-fread("c:/gwas_summary/alc/Liu2019/DrinksPerWeek.WithoutUKB.txt.gz") %>%
  format_data(. , type="outcome", chr_col="CHROM", phenotype_col="exposure", snp_col="RSID",
              beta_col="BETA", se_col="SE", eaf_col="AF",
              effect_allele_col="ALT",
              other_allele_col="REF", pval="PVALUE") %>% 
  mutate(outcome="dpw", id.outcome="dpw",Phenotype="dpw") %>% 
  dplyr::rename(beta=beta.outcome,se=se.outcome,
                effect_allele=effect_allele.outcome,
                other_allele=other_allele.outcome,
                eaf=eaf.outcome,
                pval=pval.outcome) %>% 
  dplyr::select(SNP,outcome,Phenotype,id.outcome,beta,se,pval,effect_allele,other_allele,eaf)
nrow(dpw)
head(dpw)

exposure_dat <- mvmr_extract_exposures_local(exposure_data_list=list(cpd,dpw),
                                             pval_threshold = 5e-06,
                                             clump_r2 = 0.1)
#--------------------
# 1.0) Outcome data: Paro Shungin2019
#--------------------
pos.rsid<-fread("c:/gwas_summary/1000G-MAR2012-B37-ALL_HG19.annotation-UCSC_Function-uniqueIDs.txt.gz",
                sep="auto",header="auto") %>%
  dplyr::mutate(SNP=RSID,CHR=as.numeric(Chr),POS=as.numeric(Pos)) %>% 
  dplyr::select(SNP,CHR,POS) %>% 
  dplyr::arrange(.,POS) %>% 
  as_tibble() 
head(pos.rsid)
mean(pos.rsid$POS)
nrow(pos.rsid)
paro<-fread("c:/gwas_summary/periodontitis/Shungin2019NatureComms/EUR_perio_incl_HCHSSOL.txt",
            fill=TRUE) %>% 
  separate(.,MarkerName,sep=":",c("CHR","POS")) %>% 
  dplyr::mutate(CHR=as.numeric(CHR),POS=as.numeric(POS))
head(paro)
paro2<- right_join(paro,pos.rsid, by=c("POS","CHR")) %>% 
  dplyr::distinct(SNP,.keep_all=TRUE) %>% 
  dplyr::mutate(Allele1=toupper(Allele1), Allele2=toupper(Allele2),outcome="paro",rsid=SNP,pval=`P-value`) %>% 
  dplyr::mutate(chr.outcome=CHR,
                effect_allele.outcome=Allele1, 
                other_allele.outcome=Allele2,
                beta.outcome=Effect,
                se.outcome=StdErr,
                pval.outcome=`P-value`, 
                samplesize.outcome=N,
                ncase.outcome=n_cases,
                ncontrol.outcome=n_controls,
                outcome=outcome,
                id.outcome=outcome)  %>% 
  as_tibble()
rm(pos.rsid)
rm(paro)
table(exposure_dat$SNP %in% paro2$SNP)
fwrite(paro2, "paro2.csv", row.names=TRUE, quote=TRUE)
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename="paro2.csv",
  sep=",",
  phenotype_col="outcome",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col="se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col= "other_allele.outcome",
  eaf_col = "eaf.outcome",
  pval_col = "pval.outcome"
) 
#----------------------------------------------
# 1.2) Harmonise the exposure and outcome data
#----------------------------------------------
table(exposure_dat$SNP %in% outcome_dat$SNP)
dat <- mv_harmonise_data(exposure_dat, outcome_dat) 
#-----------------
# 1.3) Perform MR
#-----------------
res <- mv_multiple(dat)  
res<-res %>% 
  as.data.frame() %>% 
  dplyr::rename(b=result.b,se=result.se,pval=result.pval,
                exposure=result.exposure,outcome=result.outcome) %>% 
  mutate(OR=round(exp(b),2),pval=round(pval,4)) %>% 
  filter(exposure=="cpd" | exposure=="dpw") 
res$CI<-paste("(",round(exp(res$b-1.96*res$se),3),";",round(exp(res$b+1.96*res$se),3),")",sep="")
res<-res %>% 
  dplyr::select(exposure,outcome,OR,CI,pval)
res

mvmr_cpd_dpw_paro<-res %>%  as.data.frame()  
save(mvmr_cpd_dpw_paro,file="mvmr_cpd_dpw_paro.RData")

mvmr_cpd_dpw_paro<-mvmr_cpd_dpw_paro
mvmr_cpd_dpw_paro<-flextable(mvmr_cpd_dpw_paro) %>%
  theme_booktabs() %>%
  autofit()
print(mvmr_cpd_dpw_paro)
mvmr_cpd_dpw_paro <- read_docx() %>%
  body_add_flextable(value = mvmr_cpd_dpw_paro)
fileout3 <- tempfile(fileext = ".docx")
fileout3 <- "mvmr_cpd_dpw_paro.docx"
print(mvmr_cpd_dpw_paro, target = fileout3)


# # # Remove created outcome datasets
delfiles <- dir(path=mydir ,pattern="*.csv")
file.remove(file.path(mydir, delfiles))
