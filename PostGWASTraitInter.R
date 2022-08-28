library(tidyverse)
library(dplyr)
library(TwoSampleMR)
library(readr)
#----FUNCTIONALITY----
#Adds rsid to summary statistics if a mapping file with chromosomal positions is supplied
#Performs bulk MR analysis in both directions with specific summary statistics from TwoSampleMR package or with summary statistics for phenotypes containing specific substrings in their names
#Does LDSC preprocessing
#Performs local genetic covariation analysis with SUPERGNOVA

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

two_sample_mr.backward.singular<-function(outcome_sumstat=processed_sumstat_name, exposure_id, exposure_name, exposure_threshold_p=MR.BACKWARD.THRESHOLD_PVAL, eval_singlesnp=TRUE, onlysum=MR.SAVE_ONLYSUMMARY){
  print(paste('Exposure:',exposure_name))
  
  proceed=FALSE
  tryCatch({
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
  proceed=TRUE
  }, error = function(e){
    print('Network error...')
  }
  )
  
  if (proceed==TRUE){
    result<-list()
    tryCatch({
      dat <- harmonise_data(instruments, outcomes) 
      res <- generate_odds_ratios(mr(dat))
      
      
      if (onlysum==FALSE){
        scatter<-mr_scatter_plot(res, dat)
        pleiotropy_test<-mr_pleiotropy_test(dat)
        pleiotropy_test_p=pleiotropy_test$pval
        heterogeneity_test<-mr_heterogeneity(dat)
        directionality_testres<-directionality_test(dat)
        
        res<-res %>% 
          mutate(exposure.pipename=exposure_name, exposure_threshold_pval=exposure_threshold_p, pleiotropy_test_pval=pleiotropy_test_p, directionality_test_pval=directionality_testres$pval)
        res<-left_join(res, heterogeneity_test)
        result$scatter<-list(scatter)
        result$pleiotropy_test<-pleiotropy_test
        res$directionality_test<-directionality_testres
        eval_singlesnp=FALSE
      } else{
        res<-res %>% 
          mutate(exposure.pipename=exposure_name, exposure_threshold_pval=exposure_threshold_p)
        
      }
      
      result$dat<-dat
      result$res<-res
      
      if (eval_singlesnp==TRUE){
        res_single <- mr_singlesnp(dat)
        res_single_forestplot<- mr_forest_plot(res_single)
        result$single_res<-res_single
        result$single_forest<-list(res_single_forestplot)
      }
   
    }, error = function(e){
      print('Not enough SNPs...')
      }
      )
    
    if (length(result)==0){
      'Undefined'
    } else {
      result
    }
  } else {
    'Undefined'
  }
}

two_sample_mr.forward.singular<-function(exposure_sumstat=processed_sumstat_name, exposure_sumstat_mode='obj', outcome_id, outcome_name, exposure_threshold_p=MR.FORWARD.THRESHOLD_PVAL, eval_singlesnp=TRUE, onlysum=MR.SAVE_ONLYSUMMARY){
  if (!exposure_sumstat_mode=='obj'){
    instruments<-read_exposure_data(filename=exposure_sumstat,
                                        sep='\t',
                                        snp_col = 'rsid',
                                        beta_col = "effect",
                                        se_col='se',
                                        effect_allele_col='A1',
                                        pval_col = 'pval',
                                        samplesize_col = 'npersons',
                                        other_allele_col = 'A2')
    instruments<-instruments %>% 
      filter(pval.exposure<exposure_threshold_p)
  } else {
    instruments<-exposure_sumstat
  }
  outcomes<-extract_outcome_data(snps=instruments$SNP, outcomes = outcome_id)
  #print(outcomes)
  result<-list()
  tryCatch({
    dat <- harmonise_data(instruments, outcomes) 
    res <- generate_odds_ratios(mr(dat))
    #print(res)
    if (onlysum==FALSE){
      scatter<-mr_scatter_plot(res, dat)
      pleiotropy_test<-mr_pleiotropy_test(dat)
      pleiotropy_test_p=pleiotropy_test$pval
      heterogeneity_test<-mr_heterogeneity(dat)
      directionality_testres<-directionality_test(dat)
      res<-res %>% 
        mutate(outcome.pipename=outcome_name, exposure_threshold_pval=exposure_threshold_p, pleiotropy_test_pval=pleiotropy_test_p, directionality_test_pval=directionality_testres$pval)
      res<-left_join(res, heterogeneity_test)
      result$scatter<-list(scatter)
      result$pleiotropy_test<-pleiotropy_test
      result$directionality_test<-directionality_testres
      eval_singlesnp=FALSE
    } else {
    res<-res %>% 
      mutate(outcome.pipename=outcome_name, exposure_threshold_pval=exposure_threshold_p)
    }
    
    result$dat<-dat
    result$res<-res
    
    
    if (eval_singlesnp==TRUE){
      res_single <- mr_singlesnp(dat)
      res_single_forestplot<- mr_forest_plot(res_single)
      result$single_res<-res_single
      result$single_forest<-list(res_single_forestplot)
    }
    
  }, error = function(e){
    print('Not enouht SNPs...')
  }
  )
  
  if (length(result)==0){
    'Undefined'
  } else {
    result
  }
}



two_sample_mr.backward.serial<-function(outcome_sumstat=processed_sumstat_name, exposure_ids=MR.COMPARISON_STUDIES, exposure_names=MR.COMPARISON_STUDIES.NAMES, exposure_threshold_p=MR.BACKWARD.THRESHOLD_PVAL, eval_singlesnp=TRUE, onlysum=MR.SAVE_ONLYSUMMARY){
  print('Starting MR evaluation...')
  allres<-list()
  result_fulldf<-tibble()
  print('Starting two sample MR experiment series...')
  
  for (i in c(1:length(exposure_ids))){
    if (i%%10==0){
      print(paste((i/total_nsstats)*100, '% complete...', sep=''))
    }
    exposure_id=exposure_ids[i]
    exposure_name<-exposure_names[i]
    cur_result<-two_sample_mr.backward.singular(exposure_id=exposure_id, exposure_name=exposure_name, exposure_threshold_p = exposure_threshold_p, eval_singlesnp = eval_singlesnp, onlysum=onlysum)
    if (!cur_result=='Undefined'){
      cur_result_df=cur_result$res
      if (length(allres)==0){
        result_fulldf = cur_result_df
      } else {
        result_fulldf = rbind(result_fulldf, cur_result_df)
      }
      allres[exposure_id]=list(cur_result)
    } 
  }
  allres$full_dataframe=result_fulldf
  allres
}




two_sample_mr.forward.serial<-function(exposure_sumstat=processed_sumstat_name, outcome_ids=MR.COMPARISON_STUDIES, outcome_names=MR.COMPARISON_STUDIES.NAMES, exposure_threshold_p=MR.FORWARD.THRESHOLD_PVAL, eval_singlesnp=TRUE, onlysum=MR.SAVE_ONLYSUMMARY){
  print('Starting MR evaluation...')
  allres<-list()
  result_fulldf<-tibble()
  print('Starting two sample MR experiment series...')
  
  instruments<-read_exposure_data(filename=exposure_sumstat,
                                  sep='\t',
                                  snp_col = 'rsid',
                                  beta_col = "effect",
                                  se_col='se',
                                  effect_allele_col='A1',
                                  pval_col = 'pval',
                                  samplesize_col = 'npersons',
                                  other_allele_col = 'A2')
  instruments<-instruments %>% 
    filter(pval.exposure<exposure_threshold_p)
  
  print(paste('Instrument count:', length(instruments$SNP)))
  
  
  for (i in c(1:length(outcome_ids))){
    if (i%%10==0){
      print(paste((i/total_nsstats)*100, '% complete...', sep=''))
    }
    outcome_id=outcome_ids[i]
    outcome_name<-outcome_names[i]
    cur_result<-two_sample_mr.forward.singular(exposure_sumstat = instruments, exposure_sumstat_mode = 'obj',outcome_id=outcome_id, outcome_name=outcome_name, exposure_threshold_p = exposure_threshold_p, eval_singlesnp = eval_singlesnp, onlysum=onlysum)
    if (!cur_result=='Undefined'){
      cur_result_df=cur_result$res
      if (length(allres)==0){
        result_fulldf = cur_result_df
      } else {
        result_fulldf = rbind(result_fulldf, cur_result_df)
      }
      allres[outcome_id]=list(cur_result)
    } 
  }
  allres$full_dataframe=result_fulldf
  allres
}







#-----SETTINGS-----
#----SETTINGS: General----
gwas_sumstatfile_path<-'/home/biorp/Gitrepos/Psychiatry/SUMSTATS/anged/anged/dop_pheno_all.anged.glm.logistic.hybrid'
annotationfile_path<-'/home/biorp/Gitrepos/Psychiatry/SUMSTATS/psychiatric_genomics_chip_mapping.tsv' #file with mapping columns: CHR, POS, rsid
ldsc_munge_sumstats_path<-'/home/biorp/TOOLS/ldsc/munge_sumstats.py'
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





#----SETTINGS: Mendelian Randomization----
#Two Sample
mr_ts_backward_resname=paste(gwasname, '_twosample_backward.rds', sep='')
mr_ts_forward_resname=paste(gwasname, '_twosample_forward.rds', sep='')
MR.FORWARD=TRUE
MR.BACKWARD=TRUE
ao<-available_outcomes()
MR.SAVE_ONLYSUMMARY=TRUE
MR.COMPARISON_STUDIES.MODE='PHENONAMES' #'PHENONAMES' or 'SPECIFIC'


MR.COMPARISON_STUDIES.RUN_PHENONAMES<-c('schiz','depress','bipol','cholester','lipo','diab','fatty','covid','respir','serotonin','cholin','dopamin','anxiety','neuro','mental','brain')
MR.COMPARISON_STUDIES.RUN_SUBCATEGORIES<-c('Fatty acid','Cofactors and vitamins','Lipid','Amino acid','Carbohydrate','Metal','Personality','Unknown metabolite','Peptide','Protein','Behavioural','Sleeping','Diabetes','Psychiatric / neurological','Autoimmune / inflammatory','Immune system','Education','Metabolites ratio','Biomarker')

MR.COMPARISON_STUDIES.RUN_SPECIFIC=c() #indexes of studies from ao
MR.PHENONAMES.NCASE_THRESHOLD=4000
MR.FORWARD.THRESHOLD_PVALS=c(5e-5) #threshold p-value for instrument extraction from the studied sumstat.
MR.BACKWARD.THRESHOLD_PVALS=c(5e-8)

unique(ao$subcategory)

#-----SETTINGS: Local Covariance (SUPERGNOVA)-----
#Local genetic covariance
LOCGCOVAR=TRUE
LOCGCOVAR.COMPARISON_STUDIES.SUMSTAT_PATHS=c()
LOCGCOVAR.COMPARISON_STUDIES.NAMES=c()
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
if (MR.COMPARISON_STUDIES.MODE=='PHENONAMES'){
  MR.COMPARISON_STUDIES<-c()
  for (substring in MR.COMPARISON_STUDIES.RUN_PHENONAMES){
    MR.COMPARISON_STUDIES_PT=ao %>%
      filter(grepl(substring, trait, ignore.case = TRUE)) %>% 
      filter(ncase>MR.PHENONAMES.NCASE_THRESHOLD) %>% 
      select(id) %>% 
      pull()
    MR.COMPARISON_STUDIES=c(MR.COMPARISON_STUDIES, MR.COMPARISON_STUDIES_PT)
  }
  for (subcat in MR.COMPARISON_STUDIES.RUN_SUBCATEGORIES){
    MR.COMPARISON_STUDIES_PT=ao %>%
      filter(subcategory==subcat) %>% 
      filter(ncase>MR.PHENONAMES.NCASE_THRESHOLD) %>% 
      select(id) %>% 
      pull()
    MR.COMPARISON_STUDIES=c(MR.COMPARISON_STUDIES, MR.COMPARISON_STUDIES_PT)
  }
  MR.COMPARISON_STUDIES=unique(MR.COMPARISON_STUDIES)
  } else {
  MR.COMPARISON_STUDIES=MR.COMPARISON_STUDIES.RUN_SPECIFIC
}
MR.COMPARISON_STUDIES.NAMES=rep('default', length(MR.COMPARISON_STUDIES))


if (MR.BACKWARD==TRUE){
  if (!file.exists(mr_ts_backward_resname)){
    print('Performing two sample MR analysis in backward direction (comparison phenotypes->target summary statistic)')
    allres<-list()
    allresdf<-'Undefined'
    total_nsstats=length(MR.BACKWARD.THRESHOLD_PVALS)
    for (ind in c(1:length(MR.BACKWARD.THRESHOLD_PVALS))){
      
      MR.BACKWARD.THRESHOLD_PVAL=MR.BACKWARD.THRESHOLD_PVALS[ind]
      print(paste('Exposure p-value:', MR.BACKWARD.THRESHOLD_PVAL))
      sectionname<-paste('exposure_p_threshold', MR.BACKWARD.THRESHOLD_PVAL, sep='_')
      mr_results<-two_sample_mr.backward.serial(outcome_sumstat=processed_sumstat_name, exposure_ids=MR.COMPARISON_STUDIES, exposure_names=MR.COMPARISON_STUDIES.NAMES, exposure_threshold_p=MR.BACKWARD.THRESHOLD_PVAL, eval_singlesnp=FALSE, onlysum=MR.SAVE_ONLYSUMMARY)
      if (length(allres)==0){
        allresdf=mr_results$full_dataframe
      } else {
        allresdf<-rbind(allresdf, mr_results$full_dataframe)
      }
      if (MR.SAVE_ONLYSUMMARY==FALSE){
        allres[sectionname]<-list(mr_results)
      }
    }
    allres$resdf<-list(allresdf)
    print('Significant interactions:')
    print(allres$resdf[[1]]%>% 
      filter(pval<0.05) %>% 
      filter(method=='Inverse variance weighted'))
    saveRDS(allres, file=mr_ts_backward_resname)
  }
}



if (MR.FORWARD==TRUE){
  if (!file.exists(mr_ts_forward_resname)){
    print('Performing two sample MR analysis in backward direction (comparison phenotypes->target summary statistic)')
    allres<-list()
    allresdf<-'Undefined'
    total_nsstats=length(MR.FORWARD.THRESHOLD_PVALS)
    for (ind in c(1:length(MR.FORWARD.THRESHOLD_PVALS))){
      
      MR.FORWARD.THRESHOLD_PVAL=MR.FORWARD.THRESHOLD_PVALS[ind]
      print(paste('Exposure p-value:', MR.FORWARD.THRESHOLD_PVAL))
      sectionname<-paste('exposure_p_threshold', MR.FORWARD.THRESHOLD_PVAL, sep='_')
      mr_results<-two_sample_mr.forward.serial(exposure_sumstat=processed_sumstat_name, outcome_ids=MR.COMPARISON_STUDIES, outcome_names=MR.COMPARISON_STUDIES.NAMES, exposure_threshold_p=MR.FORWARD.THRESHOLD_PVAL, eval_singlesnp=FALSE)
      if (length(allres)==0){
        allresdf=mr_results$full_dataframe
      } else {
        allresdf<-rbind(allresdf, mr_results$full_dataframe)
      }
      if (MR.SAVE_ONLYSUMMARY==FALSE){
        allres[sectionname]<-list(mr_results)
      }
    }
    allres$resdf<-list(allresdf)
    print('Significant interactions:')
    print(allres$resdf[[1]]%>% 
      filter(pval<0.05) %>% 
      filter(method=='Inverse variance weighted'))
    saveRDS(allres, file=mr_ts_forward_resname)
  }
}




#-----3. Preprocessing summary statistics with LDSC-----
system('python ')

python munge_sumstats.py --sumstats /home/biorp/Gitrepos/DepressionRevision/SUMSTATS/ieu_a_31.summary.tsv --out ieu_a_31_IBD --merge-alleles w_hm3.snplist --snp rsid --N-col n --a1 A1 --a2 A2 --p P
ldsc_munge_sumstats_path

