#+ echo = TRUE, warning = FALSE, message = FALSE
#' --------- 
#' title: "MR Cigarettes per day (Lui2019) and Peridontitis  (SHRINE2019)
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
setwd("c:/arbeit/Publikationen/mr_smoking_alc_paro/data/dpwLui2019noUKBB_paroShugin2019/")
mydir <- getwd()
#--------------------
# 1.0) Exposure data: Drinks per week, Lui 2019, PMID  30643251, 
#--------------------
dpw<-fread("c:/gwas_summary/alc/Liu2019/DrinksPerWeek.txt.gz") %>%
  mutate(exposure="dpw") %>%
  dplyr::rename(SNP=RSID)
dpw2<-dpw %>%  format_data(. , type="exposure", chr_col="CHROM", phenotype_col="exposure", snp_col="SNP",
                           beta_col="BETA", se_col="SE", eaf_col="AF",
                           effect_allele_col="ALT",
                           other_allele_col="REF", pval="PVALUE") %>% 
  mutate(rsid=SNP,pval=pval.exposure)
exposure_dat_sel<-subset(dpw2,pval.exposure<5e-6) %>% 
  ieugwasr::ld_clump(clump_r2 = 0.1)
nrow(exposure_dat_sel)
head(exposure_dat_sel)
exposure_dat_sel_snp<- exposure_dat_sel %>% 
  dplyr::select(SNP)
nrow(exposure_dat_sel_snp)
save(exposure_dat_sel_snp,file="exposure_dat_sel_snp.Rdata")
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
table(exposure_dat_sel$SNP %in% paro2$SNP)
n<-max(paro2$N, na.rm=TRUE)
cases<-max(paro2$n_cases, na.rm=TRUE)
controls<-max(paro2$n_controls, na.rm=TRUE)
k<-cases/n
ratio <-cases/controls
# # # # Get eaf
# mart<-useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
# snp_attributes<-c("refsnp_id", "chr_name", "minor_allele", "minor_allele_freq")
# eaf<-getBM(attributes=snp_attributes, filters="snp_filter", values=exposure_dat_sel_snp$SNP,
#            mart=mart)  %>%
#   dplyr::rename(SNP=refsnp_id) %>%
#   dplyr::select(SNP,minor_allele_freq)
# head(eaf)
# paro3 <- paro2 %>%
#   dplyr::select(!(eaf.outcome))
# paro4 <- left_join(paro3,eaf,by="SNP")  %>%
#   dplyr::rename(eaf.outcome=minor_allele_freq) %>%
#   distinct()
#' # calculate R2 and F-statistic per SNP (IV)
colnames(exposure_dat_sel)
R2_F<-calculate_r2_F(exposure_dat_sel$eaf.exposure,exposure_dat_sel$beta.exposure,exposure_dat_sel$se.exposure,n)
R2_F
exposure_dat_sel$F<-R2_F$Fstat
exposure_dat_sel$R2<-R2_F$R2
min(exposure_dat_sel$F)
#' #' #' Power calculation according to Brion MJ, Int J Epidemiol 2013 (http://cnsgenomics.com/shiny/mRnd/); PMID=24159078
r2.xz<-sum(R2_F$R2)
r2.xz*100
expected.or<-c(0.90,0.85,0.8,0.78,0.75,0.7,0.65)
pow<-sapply(expected.or,function(x){calculate_power_binary(x,k,r2.xz,n)})
data.frame(list(expected.OR=expected.or,Power=round(pow,2)))
fwrite(paro2, "paro2.csv", row.names=TRUE, quote=TRUE)
outcome_dat <- read_outcome_data(
  snps = exposure_dat_sel$SNP,
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
table(exposure_dat_sel$SNP %in% outcome_dat$SNP)
head(outcome_dat)
#outcome_dat<-outcome_dat %>%   filter(pval.outcome>0.05) 
dat <- harmonise_data(exposure_dat=exposure_dat_sel, outcome_dat=outcome_dat,action=2) %>% 
  mutate(id.outcome = "paro", id.exposure="dpw")
# # # # Export Supplementary Table 1	Association of instruments with drinks per week
s_table1_snp_dpw<-dat %>%
  dplyr::select(SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,
                beta.exposure,se.exposure,pval.exposure,F) %>%
  dplyr::rename(EA=effect_allele.exposure,
                OA=other_allele.exposure,
                EAF=eaf.exposure,
                BETA=beta.exposure,
                SE=se.exposure,
                P=pval.exposure) %>%
  dplyr::mutate(EAF=round(EAF,3),
                BETA=round(BETA,3),
                SE=round(SE,3),
                P=formatC(P,format="e",digits=3),
                F=round(F,1))
nrow(s_table1_snp_dpw)
head(s_table1_snp_dpw)
s_table1_snp_dpw
save(s_table1_snp_dpw,file="s_table1_snp_dpw.RData")
s_table1_snp_dpw_X<-s_table1_snp_dpw
s_table1_snp_dpw<-flextable(s_table1_snp_dpw) %>%
  theme_booktabs() %>%
  autofit()
print(s_table1_snp_dpw)
s_table1_snp_dpw <- read_docx() %>%
  body_add_flextable(value = s_table1_snp_dpw)
fileout3 <- tempfile(fileext = ".docx")
fileout3 <- "s_table1_snp_dpw.docx"
print(s_table1_snp_dpw, target = fileout3)
## Export Supplementary Table 2	Association of instruments with paro
s_table2_dpw_snp_paro<-s_table1_snp_dpw_X %>%
  dplyr::select(SNP) %>%
  inner_join(.,dat,by="SNP") %>%
  dplyr::rename(BETA=beta.outcome,
                SE=se.outcome,
                Pvalue=pval.outcome) %>%
  dplyr::select(SNP,BETA,SE,Pvalue) %>%
  dplyr::mutate(BETA=round(BETA,3),
                SE=round(SE,3),
                Pvalue=round(Pvalue,3))
head(s_table2_dpw_snp_paro)
nrow(s_table2_dpw_snp_paro)
s_table2_dpw_snp_paro<-flextable(s_table2_dpw_snp_paro) %>%
  theme_booktabs() %>%
  autofit()
print(s_table2_dpw_snp_paro)
s_table2_dpw_snp_paro <- read_docx() %>%
  body_add_flextable(value = s_table2_dpw_snp_paro)
fileout3 <- tempfile(fileext = ".docx")
fileout3 <- "s_table2_dpw_snp_paro.docx"
print(s_table2_dpw_snp_paro, target = fileout3)
#-----------------
# 1.3) Perform MR
#-----------------
table(dat$palindromic); table(dat$ambiguous); table(dat$remove); table(dat$mr_keep)
table(sign(dat$beta.exposure)==sign(dat$beta.outcome))
res <- mr(dat, method_list=c("mr_ivw_mre","mr_penalised_weighted_median","mr_raps", "mr_ivw_radial"))
res$method<-as.character(res$method)
signthresh<-0.05
resPRESSO<-run_mr_presso(dat,NbDistribution = 1000, SignifThreshold	= signthresh)
if(is.na(resPRESSO[[1]]$"Main MR results"[2,3])){
  buf<-resPRESSO[[1]]$"Main MR results"[1,c("MR Analysis","Causal Estimate","Sd","P-value")]
  buf$nsnp<-nrow(dat)
}else{
  buf<-resPRESSO[[1]]$"Main MR results"[2,c("MR Analysis","Causal Estimate","Sd","P-value")]
  buf$nsnp<-nrow(resPRESSO[[1]]$"MR-PRESSO results"$"Outlier Test")-length(resPRESSO[[1]]$"MR-PRESSO results"$"Distortion Test"$"Outliers Indices")
}
names(buf)<-c("method","b","se","pval","nsnp")
buf$method<-paste("MR PRESSO:",buf$method,sep=" ")
buf<-cbind(buf,res[1,c("id.exposure","id.outcome")])
res <- res %>%  dplyr::select(method, b, se, pval, nsnp, id.exposure, id.outcome)
res2<-rbind(res,buf)
res<-res2
res
res$OR<-round(exp(res$b),2)
res$CI<-paste("(",round(exp(res$b-1.96*res$se),3),";",round(exp(res$b+1.96*res$se),3),")",sep="")
res
mr_dpw_paro<-res %>%  as.data.frame()  
save(mr_dpw_paro,file="mr_dpw_paro.RData")
##################################
# heterogeneity statistics
##################################
res_single <- mr_singlesnp(dat, all_method = c("mr_egger_regression")) 
I2<-round(Isq(res_single$b,res_single$se),digits=2) %>% as.data.frame()
het<-mr_heterogeneity(dat, method_list=c("mr_ivw_mre")) %>% 
  dplyr::rename(Exposure=id.exposure,Outcome=id.outcome) %>%
  dplyr::mutate(Q=round(Q,digits=1), Q_pval=round(Q_pval,3))  %>%  
  bind_cols(.,I2) %>% 
  dplyr::rename(I2=".") %>% 
  dplyr::select(Exposure,Outcome,method,Q,Q_df,Q_pval,I2) %>%  
  as.data.frame()
het
het_dpw_paro<-het %>%  as.data.frame()  
save(het_dpw_paro,file="het_dpw_paro.RData")
# directional pleiotropy
plei<-mr_pleiotropy_test(dat) %>% 
  dplyr::rename(Exposure=id.exposure,Outcome=id.outcome) %>% 
  mutate(Intercept=formatC(egger_intercept,format="e",digits=3),
         SE=round(se,digits=3),Pval=formatC(pval,format="e",digits=3)) %>%  
  dplyr::select(Outcome,Exposure,Intercept,SE,Pval) %>% 
  as.data.frame()
plei
plei_dpw_paro<-plei %>%  as.data.frame()  
save(plei_dpw_paro,file="plei_dpw_paro.RData")
# leave one out
res_loo <- mr_leaveoneout(dat,   method=mr_ivw_mre)
res_loo$OR<-round(exp(res_loo$b),2)
res_loo$CI<-paste("(",round(exp(res_loo$b-1.96*res_loo$se),2),";",round(exp(res_loo$b+1.96*res_loo$se),2),")",sep="")
res_loo<-res_loo %>% dplyr::rename(Pvalue=p)
res_loo$Qvalue<-p.adjust(res_loo$Pvalue,method="BH")
res_loo
res_loo<-res_loo  %>%  
  dplyr::select(exposure,outcome,SNP,OR,CI,Pvalue,Qvalue) %>% 
  mutate(OR=round(OR,3), 
         Pvalue=formatC(Pvalue,format="e",digits=3),
         Qvalue=formatC(Qvalue,format="e",digits=3)) %>% 
  setDT()
res_loo
loo_dpw_paro<-res_loo %>%  as.data.frame()  
save(loo_dpw_paro,file="loo_dpw_paro.RData")
# # # Remove created outcome datasets
delfiles <- dir(path=mydir ,pattern="*.csv")
file.remove(file.path(mydir, delfiles))
