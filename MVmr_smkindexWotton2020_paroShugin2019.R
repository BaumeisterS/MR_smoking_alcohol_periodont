#+ echo = TRUE, warning = FALSE, message = FALSE
#' --------- 
#' title: "MVMR Smoking index (Wotton 2020) and Peridontitis  (SHRINE2019) - adjusted for bmi and education
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
sessionInfo()
# Load MR functions
source("c:/R_fct/mr_fcts.R")
ls()
setwd("c:/arbeit/Publikationen/mr_smoking_alc_paro/data/smkindexWotton2020_paroShugin2019_MV/")
mydir <- getwd()
#--------------------
# 1.0) Exposure data: Smoking index, Wotton 2020, PMID  31689377,  https://data.bris.ac.uk/data/dataset/10i96zb8gm0j81yz0q6ztei23d
#--------------------
smidx<-fread("c:/gwas_summary/smoking/Wootton2020/2019.10.02 Lifetime Smoking GWAS Data Sheet 1.txt.gz") %>% 
  format_data(. , type="outcome", chr_col="CHR", phenotype_col="exposure", snp_col="SNP",
                           beta_col="BETA", se_col="SE", eaf_col="EAF",
                           effect_allele_col="EFFECT_ALLELE",
                           other_allele_col="OTHER_ALLELE", pval="P") %>% 
  mutate(outcome="smidx", id.outcome="smidx",Phenotype="smidx") %>% 
  dplyr::rename(beta=beta.outcome,se=se.outcome,
                effect_allele=effect_allele.outcome,
                other_allele=other_allele.outcome,
                eaf=eaf.outcome,
                pval=pval.outcome) %>% 
  dplyr::select(SNP,outcome,Phenotype,id.outcome,beta,se,pval,effect_allele,other_allele,eaf)
smidx
#----------------------------------------------
# Education Lee 2019
#----------------------------------------------
edu<-fread("c:/gwas_summary/educ/GWAS_EA_excl23andMe.txt")  %>%
  format_data(. , type="outcome", chr_col="CHR", phenotype_col="exposure", snp_col="MarkerName",
              beta_col="Beta", se_col="SE", eaf_col="EAF",
              effect_allele_col="A1",
              other_allele_col="A2", pval="Pval") %>% 
  mutate(outcome="edu", id.outcome="edu",Phenotype="edu") %>% 
  dplyr::rename(beta=beta.outcome,se=se.outcome,
                effect_allele=effect_allele.outcome,
                other_allele=other_allele.outcome,
                eaf=eaf.outcome,
                pval=pval.outcome) %>% 
  dplyr::select(SNP,outcome,Phenotype,id.outcome,beta,se,pval,effect_allele,other_allele,eaf)
edu
colnames(edu)
class(edu$beta)


#----------------------------------------------
# BMI  Pulit 2019
#----------------------------------------------
# bmi<-fread("c:/gwas_summary/bmi/pulit2019/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz") %>%
#   separate(.,SNP,sep=":",c("SNP","A1","A2")) %>% 
 load("c:/gwas_summary/bmi/pulit2019/bmi.Rdata")  
 bmi<-bmi %>% 
  format_data(. , type="outcome", chr_col="CHR", phenotype_col="exposure", snp_col="SNP",
              beta_col="BETA", se_col="SE", eaf_col="Freq_Tested_Allele",
              effect_allele_col="Tested_Allele",
              other_allele_col="Other_Allele", pval="P") %>% 
  mutate(outcome="bmi", id.outcome="bmi",Phenotype="bmi") %>% 
  dplyr::rename(beta=beta.outcome,se=se.outcome,
                effect_allele=effect_allele.outcome,
                other_allele=other_allele.outcome,
                eaf=eaf.outcome,
                pval=pval.outcome) %>% 
  dplyr::select(SNP,outcome,Phenotype,id.outcome,beta,se,pval,effect_allele,other_allele,eaf)
bmi
colnames(bmi)
class(bmi$se)


exposure_dat <- mvmr_extract_exposures_local(exposure_data_list=list(smidx,edu,bmi),clump_r2 = 0.01)
head(exposure_dat)
nrow(exposure_dat)
#rm(list=c("smidx","edu","bmi2"))
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
res
res<-res %>% 
  as.data.frame() %>% 
  dplyr::rename(b=result.b,se=result.se,pval=result.pval,
                exposure=result.exposure,outcome=result.outcome) %>% 
  mutate(OR=round(exp(b),2),pval=round(pval,4)) %>% 
  filter(exposure=="smidx") 
res$CI<-paste("(",round(exp(res$b-1.96*res$se),3),";",round(exp(res$b+1.96*res$se),3),")",sep="")
res<-res %>% 
  dplyr::select(exposure,outcome,OR,CI,pval)
res
mvmr_smidx_paro<-res %>%  as.data.frame()  
save(mvmr_smidx_paro,file="mvmr_smidx_paro.RData")
mvmr_smidx_paro<-mvmr_smidx_paro
mvmr_smidx_paro<-flextable(mvmr_smidx_paro) %>%
  theme_booktabs() %>%
  autofit()
print(mvmr_smidx_paro)
mvmr_smidx_paro <- read_docx() %>%
  body_add_flextable(value = mvmr_smidx_paro)
fileout3 <- tempfile(fileext = ".docx")
fileout3 <- "mvmr_smidx_paro.docx"
print(mvmr_smidx_paro, target = fileout3)



# # # Remove created outcome datasets
delfiles <- dir(path=mydir ,pattern="*.csv")
file.remove(file.path(mydir, delfiles))
