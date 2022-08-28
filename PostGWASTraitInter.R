library(tidyverse)
library(dplyr)
library(TwoSampleMR)
library(readr)
library(readxl)

#-----FUNCTIONS------
unify_sumstat_colnames<-function(sumstat, outcomename=gwasname,rsid_coln=rsid_col, effect_coln=effect_col, se_coln=se_col, pval_coln=pval_col, effect_allele_coln=effect_allele_col, other_allele_coln=other_allele_col, n=npersons){
  print('Unifying summary statistic column names...')
  sumstat$rsid=unlist(sumstat[rsid_coln])
  sumstat$effect=unlist(sumstat[effect_coln])
  sumstat$se=unlist(sumstat[se_coln])
  sumstat$pval=unlist(sumstat[pval_coln])
  sumstat$A1=unlist(sumstat[effect_allele_coln])
  sumstat$A2=unlist(sumstat[other_allele_coln])
  sumstat %>% 
    mutate(npersons=n, phenotype=gwasname) %>% 
    dplyr::select(rsid, effect, se, pval, A1, A2, npersons, phenotype)
}

two_sample_mr.backward.singular<-function(outcome_sumstat=processed_sumstat_name, exposure_id, exposure_name, exposure_threshold_p=MR.BACKWARD.THRESHOLD_PVAL, eval_singlesnp=TRUE){
  print('Performing two sample MR for:')
  print(paste('Exposure:',exposure_name))
  
  instruments<-extract_instruments(exposure_id, p1=exposure_threshold_p)
  outcomes<- read_outcome_data(
    snps = instruments$SNP,
    filename = processed_sumstat_name,
    sep = "\t",
    snp_col = "rsid",
    beta_col = "effect",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "pval",
    samplesize_col = "npersons",
    phenotype_col='phenotype'
  )
  dat <- harmonise_data(instruments, outcomes) 
  res <- generate_odds_ratios(mr(dat))
  scatter<-mr_scatter_plot(res, dat)
  pleiotropy_test<-mr_pleiotropy_test(dat)
  pleiotropy_test_p=pleiotropy_test$pval
  heterogeneity_test<-mr_heterogeneity(dat)
  
  res<-res %>% 
    mutate(exposure.pipename=exposure_name, exposure_threshold_pval=exposure_threshold_p, pleiotropy_test_pval=pleiotropy_test_p)
  res<-left_join(res, heterogeneity_test)
  result<-list()
  result$dat<-dat
  result$res<-res
  result$scatter<-list(scatter)
  result$pleiotropy_test<-pleiotropy_test
  if (eval_singlesnp==TRUE){
    res_single <- mr_singlesnp(dat)
    res_single_forestplot<- mr_forest_plot(res_single)
    result$single_res<-res_single
    result$single_forest<-list(res_single_forestplot)
  }
  result
}



two_sample_mr.backward.serial<-function(outcome_sumstat=processed_sumstat_name, exposure_ids=MR.COMPARISON_STUDIES, exposure_names=MR.COMPARISON_STUDIES.NAMES, exposure_threshold_p=MR.BACKWARD.THRESHOLD_PVAL, eval_singlesnp=TRUE){
  print('Starting MR evaluation...')
  allres<-list()
  result_fulldf<-'Undefined'
  print('Starting two sample MR experiment series...')
  for (i in c(1:length(exposure_ids))){
    exposure_id=exposure_ids[i]
    exposure_name<-exposure_names[i]
    cur_result<-two_sample_mr.backward.singular(exposure_id=exposure_id, exposure_name=exposure_name, exposure_threshold_p = exposure_threshold_p, eval_singlesnp = eval_singlesnp)
    cur_result_df=cur_result$res
    if (i==1){
      result_fulldf = cur_result_df
    } else{
      result_fulldf = rbind(c(result_fulldf, cur_result_df))
    }
    allres[exposure_id]=list(cur_result)
  }
  allres$full_dataframe=result_fulldf
}

#-----SETTINGS-----
gwas_sumstatfile_path<-'/home/biorp/Gitrepos/Psychiatry/SUMSTATS/anged/anged/dop_pheno_all.anged.glm.logistic.hybrid'
annotationfile_path<-'/home/biorp/Gitrepos/Psychiatry/SUMSTATS/psychiatric_genomics_chip_mapping.tsv' #file with mapping columns: CHR, POS, rsid
gwasname<-'anhedonia'
chr_col='#CHROM'
pos_col='POS'
effect_col='LOG_OR'
effect_allele_col='ALT'
other_allele_col='REF'
se_col='LOG(OR)_SE'
pval_col='P'
or_col='OR'
rsid_col='rsid'
npersons=4520
processed_sumstat_name=paste(gwasname, '_preprocessed.tsv', sep='')
mr_ts_backward_resname=paste(gwasname, '_twosample_backward.rds')

#Mendelian Randomization
#Two Sample
MR.FORWARD=FALSE
MR.BACKWARD=TRUE
ao<-available_outcomes()
MR.COMPARISON_STUDIES=c('ieu-a-1188','ieu-a-1187','ieu-b-102','ieu-b-41','ieu-b-42','ieu-b-5070','finn-b-F5_GAD','finn-b-KRA_PSY_ANXIETY','ukb-d-KRA_PSY_ANXIETY','finn-b-F5_ALLANXIOUS','ieu-a-118','ebi-a-GCST002920','ebi-a-GCST003770','ebi-a-GCST005327') #indexes of studies from ao
MR.COMPARISON_STUDIES.NAMES=c('Depression_Wray_2018_ieu-a-1188','Depression_Wray_2018_ieu-a-1187', 'Depression_Howard_2019_ieu-b-102','Bipolar_Stahl_2019_ieu-b-41','Schizophrenia_Ripke_2014_ieu-b-42', 'Schizophrenia_Lam_2019_ieu-b-5070','GAD_NA_2021_finn-b-F5_GAD','Anxiety_NA_2021_finn-b-KRA_PSY_ANXIETY', 'Anxiety_Neale_2018_ukb-d-KRA_PSY_ANXIETY', 'Anxiety_NA_2021_finn-b-F5_ALLANXIOUS','Neuroticism_deMoor_2014_ieu-a-118','Neuroticism_deMoor_2015_ebi-a-GCST002920','Neuroticism_Okbay_2016_ebi-a-GCST003770','Neuroticism_Turley_2018_ebi-a-GCST005327','Neuroticism_Luciano_2017_ebi-a-GCST005232')
MR.FORWARD.THRESHOLD_PVALS=c(5e-8) #threshold p-value for instrument extraction from the studied sumstat.
MR.BACKWARD.THRESHOLD_PVALS=c(5e-7,5e-8)
#Multivariate
MR.MULTIVARIATE=TRUE
MR.MULTIVARIATE.COMPARISON_STUDIES.SUMSTAT_PATHS=c()
MR.MULTIVARIATE.COMPARISON_STUDIES.NAMES=c()

#Local genetic covariance
LOCGCOVAR=TRUE
LOCGCOVAR.COMPARISON_STUDIES.SUMSTAT_PATHS=c()
LOCGCOVAR.COMPARISON_STUDIES.NAMES=c()

#Pathway enrichment
ENRICHMENT.PATHWAY=TRUE
ENRICHMENT.PATHWAY.THRESHOLD_PVAL=5e-8
ENRICHMENT.PATHWAY.LOCGCOVAR.GENERAL=TRUE
ENRICHMENT.PATHWAY.LOCGCOVAR.DIRECTIONAL=TRUE

#Brain region enrichment
ENRICHMENT.ABA=TRUE
ENRICHMENT.ABA.THRESHOLD_PVAL=5e-8
ENRICHMENT.ABA.LOCGCOVAR.GENERAL=TRUE
ENRICHMENT.ABA.LOCGCOVAR.DIRECTIONAL=TRUE

#Drug search
FIND.DRUGS=TRUE
FIND.DRUGS.LOCGCOVAR.GENERAL=TRUE
FIND.DRUGS.LOCGCOVAR.DIRECTIONAL=TRUE

#Signaling molecule search
FIND.SIGNALINGMOLS=TRUE
FIND.SIGNALINGMOLS.LOCGCOVAR.GENERAL=TRUE
FIND.SIGNALINGMOLS.LOCGCOVAR.DIRECTIONAL=TRUE


#-----RUN-----
#----1. Preprocessing----
if (!file.exists(processed_sumstat_name)){
  gwas<-read_tsv(gwas_sumstatfile_path)
  if (!'rsid' %in% colnames(gwas)){
    print('Annotating the summary statistics...')
    annotation<-read_tsv(annotationfile_path)
    gwas$CHR=unlist(gwas[chr_col])
    gwas$POS=unlist(gwas[pos_col])
    gwas<-left_join(gwas, annotation)
    print('Count of unmapped rsids:')
    print(sum(is.na(gwas$rsid)))
  } else {
    print('Using rs IDs from the summary statistics...')
  }
  
  if ((!'LOG_OR' %in% colnames(gwas)) & (effect_col=='LOG_OR')){
    print('Calculating log(OR)...')
    gwas$LOG_OR=log(unlist(gwas[or_col]))
  }
  
  if (!other_allele_col %in% colnames(gwas)){
    print('Getting the other allele...')
    gwas$A1=gwas[effect_allele_col]
    gwas<-gwas %>% 
      mutate(A2=ifelse(A1==REF, ALT, REF))
  }
  
  gwas<-unify_sumstat_colnames(gwas)
  gwas %>% 
    write_tsv(processed_sumstat_name)
}

#----2. Mendelian randomization (two sample)----
if (MR.BACKWARD==TRUE){
  if (!file.exists(mr_ts_backward_resname)){
    allres<-list()
    allresdf<-'Undefined'
    for (ind in c(1:length(MR.BACKWARD.THRESHOLD_PVAL))){
      MR.BACKWARD.THRESHOLD_PVAL=MR.BACKWARD.THRESHOLD_PVALS[ind]
      print(paste('Exposure p-value:', MR.BACKWARD.THRESHOLD_PVAL))
      sectionname<-paste('exposure_p_threshold', MR.BACKWARD.THRESHOLD_PVAL, sep='_')
      mr_results<-two_sample_mr.backward.serial(outcome_sumstat=processed_sumstat_name, exposure_ids=MR.COMPARISON_STUDIES, exposure_names=MR.COMPARISON_STUDIES.NAMES, exposure_threshold_p=MR.BACKWARD.THRESHOLD_PVAL, eval_singlesnp=TRUE)
      allres[sectionname]<-list(mr_results)
      if (ind==1){
        allresdf=allresdf
      } else {
        allresdf<-rbind(c(allresdf, mr_results$full_dataframe))
      }
    }
    allres$resdf<-allresdf
    write_rds(allres, file=mr_ts_backward_resname)
  }
}

