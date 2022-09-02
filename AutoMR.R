suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(TwoSampleMR))
suppressMessages(library(readr))
suppressMessages(library(igraph))
suppressMessages(library(ggsci))
suppressMessages(library(NetPathMiner))
suppressMessages(library(optparse))

#----FUNCTIONALITY----
#Adds rsid to summary statistics if a mapping file with chromosomal positions is supplied
#Performs bulk MR analysis in both directions with specific summary statistics from TwoSampleMR package or with summary statistics for phenotypes containing specific substrings in their names

#-----FUNCTIONS------
calculate_se_frombeta<-function(dfpath, betacol, tstatcol, savepath, sep='\t'){ #also works for logOR
  df<-read_delim(dfpath, delim=sep)
  df$se_calc<-unlist(df[betacol]/df[tstatcol])
  df %>% 
    write_delim(savepath, delim=sep)
}

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
          mutate(exposure.pipename=exposure_name, exposure_threshold_pval=exposure_threshold_p, pleiotropy_test_pval=pleiotropy_test_p, directionality_test_pval=directionality_testres$steiger_pval)
        res<-left_join(res, heterogeneity_test)
        result$scatter<-list(scatter)
        result$pleiotropy_test<-pleiotropy_test
        result$directionality_test<-directionality_testres
        result$res<-res
        
      } else{
        res<-res %>% 
          mutate(exposure.pipename=exposure_name, exposure_threshold_pval=exposure_threshold_p)
          eval_singlesnp=FALSE
        
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
    instruments<-clump_data(instruments)
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
        mutate(outcome.pipename=outcome_name, exposure_threshold_pval=exposure_threshold_p, pleiotropy_test_pval=pleiotropy_test_p, directionality_test_pval=directionality_testres$steiger_pval)
      res<-left_join(res, heterogeneity_test)
      result$res<-res
      result$scatter<-list(scatter)
      result$pleiotropy_test<-pleiotropy_test
      result$directionality_test<-directionality_testres
      
    } else {
    res<-res %>% 
      mutate(outcome.pipename=outcome_name, exposure_threshold_pval=exposure_threshold_p)
      eval_singlesnp=FALSE
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
  total_nsstats=length(exposure_ids)
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
  total_nsstats=length(outcome_ids)
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
  instruments<-clump_data(instruments)
  
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

GroupByVertex01 = function(Groups, spacing = 5) {
  Position = (order(Groups) + spacing*Groups)
  Angle    = Position * 2 * pi / max(Position)
  matrix(c(cos(Angle), sin(Angle)), ncol=2)
}

GroupByVertex02 = function(Groups) {
  numGroups = length(unique(Groups))
  GAngle    = (1:numGroups) * 2 * pi / numGroups
  Centers   = matrix(c(cos(GAngle), sin(GAngle)), ncol=2)
  x = y = c()
  for(i in 1:numGroups) {
    curGroup = which(Groups == unique(Groups)[i])
    VAngle = (1:length(curGroup)) * 2 * pi / length(curGroup)
    x = c(x, Centers[i,1] + cos(VAngle) / numGroups )
    y = c(y, Centers[i,2] + sin(VAngle) / numGroups)
  }
  matrix(c(x, y), ncol=2)
}


#-----ARGUMENT PARSING-----
option_list = list(
  make_option(c("-i", "--sumstat_file"), type="character", 
              help="Path to the analysed summary statistic in .tsv format.", metavar="character"),
  make_option(c("-r", "--rsid_column"), type="character", default="rsid", 
               help="Column with rs ids.", metavar="character"),
  make_option(c("-e", "--effect_column"), type="character", default="LOG_OR", 
               help="Column with effect measurements(beta or logOR).", metavar="character"),
  make_option(c("-o", "--odds_ratio_column"), type="character", default="OR", 
                help="Column with odds ratios if present and logOR are not available", metavar="character"),
  make_option(c("-s", "--standard_error_column"), type="character", default="SE", 
                help="Column with standard errors for effects.", metavar="character"),
  make_option(c("-p", "--pval_column"), type="character", default="P", 
                help="Column with p-values.", metavar="character"),
  make_option(c("", "--effect_allele_column"), type="character", default='A1', 
                help="Column with effect alleles.", metavar="character"),
  make_option(c("", "--other_allele_column"), type="character", default='A2', 
                help="Column with other (non-effect) alleles.", metavar="character"),
  make_option(c("-n", "--number_of_individuals"), type="integer", default=5000, 
                 help="Count of individuals.", metavar="integer"),
  make_option(c("-t", "--stat_column"), type="character", default="STAT", 
                 help="Column with test statistic values if standard error calculation is required", metavar="character"),
  make_option(c("-f", "--forward_direction"), type="logical", default=TRUE, 
                 help="Perform MR in forward direction (summary statistic-based exposure -> comparison studies outcomes)", metavar="logical"),
  make_option(c("-b", "--backward_direction"), type="logical", default=TRUE, 
                 help="Perform MR in backward direction (comparison studies exposures->summary statistic-based outcome). Note: this option performes two-sample MR, not multivariate one.", metavar="logical"),
  make_option(c("", "--grep_phenotypes"), type="character", default='schiz,depress,bipol,cholester,lipo,diab,fatty,covid,respir,serotonin,cholin,dopamin,anxiety,neuro,mental,brain', 
                 help="Select comparison studies with traits including the provided substrings. (For example, the value can be 'schiz,depress,bipol,cholester').", metavar="character"),
  make_option(c("", "--grep_subcategories"), type="character", default='Fatty acid,Cofactors and vitamins,Lipid,Amino acid,Carbohydrate,Metal,Personality,Unknown metabolite,Peptide,Protein,Behavioural,Sleeping,Diabetes,Psychiatric / neurological,Autoimmune / inflammatory,Immune system,Education,Metabolites ratio,Biomarker,Diabetes,Keto acid,Cancer,Glycemic,Immune system,Immune cell subset frequency,Education,Hormone,Blood pressure,Anthropometric', 
                 help="Select comparison studies in the specified subcategories. (For example, 'Fatty acid,Cofactors and vitamins,Lipid')", metavar="character"),
  make_option(c("", "--specific_studies"), type="character", default='', 
                help="Select studies with specific IEU ids. (For example, 'ieu-b-5075,ieu-b-5064,eqtl-a-ENSG00000184100')", metavar="character"),
  make_option(c("-d", "--comparison_study_n_cases_threshold"), type="integer", default=4000, 
                 help="Minimal count of cases for a comparison study.", metavar="integer"),
  make_option(c("","--forward_pvals"), type="character", default='5e-5', 
                  help="Threshold p-value(s) for instrument extraction from the input summary statistic. (For example, '5e-5,5e-8')", metavar="character"),
  make_option(c("","--backward_pvals"), type="character", default='5e-8', 
                  help="Threshold p-value(s) for instrument extraction from the comparison studies. (For example, '5e-7,5e-8')", metavar="character"),
  make_option(c("-v", "--visualize"), type="logical", default=TRUE, 
                 help="Visualize the results of the MR study in a graph", metavar="logical"),
  make_option(c("-x", "--gwas_name"), type="character", default='Default,current study', 
                 help="Name of the analysed summary statistic to use in file naming. Supply GWAS source after the name itseld. (For example, 'Anhedonia,current study')", metavar="character"),
  make_option(c("", "--log_odds_ratio_column"), type="character", default='LOG_OR', 
              help="Name of the log(OR) column (if present)", metavar="character")
  );




opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print(opt)

#-----SETTINGS-----
gwas_sumstatfile_path<-opt$sumstat_file
#annotationfile_path<-'/home/biorp/Gitrepos/Psychiatry/SUMSTATS/psychiatric_genomics_chip_mapping.tsv' #file with mapping columns: CHR, POS, rsid
gwasname<-word(opt$gwas_name, 1, sep=',')
gwassource<-word(opt$gwas_name, 2, sep=',')
effect_col=opt$effect_column
effect_allele_col=opt$effect_allele_column
other_allele_col=opt$other_allele_column
se_col=opt$standard_error_column
stat_col=opt$stat_column
pval_col=opt$pval_column
or_col=opt$odds_ratio_column
log_or_col<-opt$log_odds_ratio_column
rsid_col=opt$rsid_column
npersons=opt$number_of_individuals
processed_sumstat_name=paste(gwasname, '_preprocessed.tsv', sep='')






#Two Sample MR settings
mr_ts_backward_resname=paste(gwasname, '_twosample_backward.rds', sep='')
mr_ts_forward_resname=paste(gwasname, '_twosample_forward.rds', sep='')
mr_ts_backward_signif_resname=paste(gwasname, '_twosample_backward_signif.rds', sep='')
mr_ts_forward_signif_resname=paste(gwasname, '_twosample_forward_signif.rds', sep='')
mr_ts_resdfname=paste(gwasname, '_all_mr_res.tsv', sep='')
mr_ts_signif_resdfname=paste(gwasname, '_all_mr_res_signif.tsv', sep='')
mr_ts_graphdata<-paste(gwasname, '_all_mr_res_signif_graph.rds', sep='')

MR.FORWARD=opt$forward_direction
MR.BACKWARD=opt$backward_direction
ao<-available_outcomes()
MR.SAVE_ONLYSUMMARY=TRUE
#MR.COMPARISON_STUDIES.MODE='PHENONAMES' #'PHENONAMES' or 'SPECIFIC'
MR.RESULTS.PVAL_THRESHOLD=0.05


MR.COMPARISON_STUDIES.RUN_PHENONAMES<-strsplit(opt$grep_phenotypes, split=',', fixed=TRUE)[[1]]
MR.COMPARISON_STUDIES.RUN_SUBCATEGORIES<-strsplit(opt$grep_subcategories, split=',', fixed=TRUE)[[1]]
#print(MR.COMPARISON_STUDIES.RUN_SUBCATEGORIES)

tibble(GWAS_subcategory=MR.COMPARISON_STUDIES.RUN_SUBCATEGORIES) %>% 
  write_csv(paste(gwasname, 'target_GWAS_subcategories.csv', sep=''))


MR.COMPARISON_STUDIES.RUN_SPECIFIC=strsplit(opt$specific_studies, split=',', fixed=TRUE)[[1]] #indexes of studies from ao
MR.PHENONAMES.NCASE_THRESHOLD=opt$comparison_study_n_cases_threshold
MR.FORWARD.THRESHOLD_PVALS=as.numeric(strsplit(opt$forward_pvals, split=',', fixed=TRUE)[[1]])#threshold p-value for instrument extraction from the studied sumstat.
MR.BACKWARD.THRESHOLD_PVALS=as.numeric(strsplit(opt$backward_pvals, split=',', fixed=TRUE)[[1]])
MR.VISUALIZE=opt$visualize

#---***----
#-----RUN-----
#----1. Preprocessing----
if (!file.exists(processed_sumstat_name)){
  gwas<-read_tsv(gwas_sumstatfile_path)
  # if (!'rsid' %in% colnames(gwas)){
  #   print('Annotating the summary statistics...')
  #   annotation<-read_tsv(annotationfile_path)
  #   gwas$CHR=unlist(gwas[chr_col])
  #   gwas$POS=unlist(gwas[pos_col])
  #   gwas<-left_join(gwas, annotation)
  #   print('Count of unmapped rsids:')
  #   print(sum(is.na(gwas$rsid)))
  # } else {
  #   print('Using rs IDs from the summary statistics...')
  # }
  
  if ((!log_or_col %in% colnames(gwas)) & (effect_col=='LOG_OR')){
    print('Calculating log(OR)...')
    gwas$LOG_OR=log(unlist(gwas[or_col]))
  }
  
  if ((log_or_col %in% colnames(gwas)) & (effect_col=='LOG_OR')){
    gwas$LOG_OR=log(unlist(gwas[log_or_col]))
  }
  
  if (!other_allele_col %in% colnames(gwas)){
    print('Getting the other allele...')
    gwas$A1=gwas[effect_allele_col]
    gwas<-gwas %>% 
      mutate(A2=ifelse(A1==REF, ALT, REF))
  }
  
  if (!se_col %in% colnames(gwas)){
    print('Calculating SE...')
    gwas[se_col]<-unlist(gwas[effect_col]/gwas[stat_column])
  }
  
  gwas<-unify_sumstat_colnames(gwas)
  gwas %>% 
    write_tsv(processed_sumstat_name)
}

#----***----
#----2. Mendelian randomization (two sample) run----
#----2.1 Selection of comparison studies----
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
MR.COMPARISON_STUDIES=c(MR.COMPARISON_STUDIES, MR.COMPARISON_STUDIES.RUN_SPECIFIC)
MR.COMPARISON_STUDIES=unique(MR.COMPARISON_STUDIES)
MR.COMPARISON_STUDIES.NAMES=rep('default', length(MR.COMPARISON_STUDIES))
print('Count of comparison GWAS studies:')
print(length(MR.COMPARISON_STUDIES.NAMES))



#----2.2 MR run in backward direction-----
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


#-----2.3 MR run in foreward direction-----
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



#----2.4 Saving nominally significant results----
if (MR.BACKWARD==TRUE){
  if (!file.exists(mr_ts_backward_signif_resname)){
    print('Gathering data for significant MR hits in backward direction...')
    bckwd<-readRDS(mr_ts_backward_resname)
    significant_expids<-bckwd$resdf[[1]] %>% 
      filter(pval<MR.RESULTS.PVAL_THRESHOLD) %>% 
      filter(method=='Inverse variance weighted') %>% 
      select(id.exposure, exposure_threshold_pval)
    pvals<-significant_expids$exposure_threshold_pval
    significant_expids<-significant_expids$id.exposure
    bckwd_signif_res<-list()
    for (i in c(1:length(significant_expids))){
      expid=significant_expids[i]
      threshp=pvals[i]
      cur_result<-two_sample_mr.backward.singular(outcome_sumstat = processed_sumstat_name, exposure_id = expid, exposure_name = 'default', exposure_threshold_p = threshp, eval_singlesnp = TRUE,onlysum = FALSE)
      bckwd_signif_res[expid]<-list(cur_result)
    }
    bckwd_signif_res$resdf<-bckwd$resdf[[1]]
    saveRDS(bckwd_signif_res, file=mr_ts_backward_signif_resname)
  }
}



if (MR.FORWARD==TRUE){
  if (!file.exists(mr_ts_forward_signif_resname)){
    print('Gathering data for significant MR hits in forward direction...')
    fwd<-readRDS(mr_ts_forward_resname)
    significant_outids<-fwd$resdf[[1]] %>% 
      filter(pval<MR.RESULTS.PVAL_THRESHOLD) %>% 
      filter(method=='Inverse variance weighted') %>% 
      select(id.outcome, exposure_threshold_pval)
    pvals<-significant_outids$exposure_threshold_pval
    significant_outids<-significant_outids$id.outcome
    fwd_signif_res<-list()
    for (i in c(1:length(significant_outids))){
      outid=significant_outids[i]
      threshp=pvals[i]
      cur_result<-two_sample_mr.forward.singular(exposure_sumstat = processed_sumstat_name, exposure_sumstat_mode = 'path',outcome_id = outid,outcome_name = 'default', exposure_threshold_p = threshp, eval_singlesnp = TRUE, onlysum = FALSE)
      fwd_signif_res[outid]<-list(cur_result)
    }
    fwd_signif_res$resdf<-fwd$resdf[[1]]
    saveRDS(bckwd_signif_res, file=mr_ts_forward_signif_resname)
  }
}


#----2.5 Saving resulting dataframes-----
if (!file.exists(mr_ts_resdfname)){
  print('Saving MR results...')
  if (MR.FORWARD==TRUE && MR.BACKWARD==TRUE){
    bckwd<-readRDS(mr_ts_backward_resname)$resdf[[1]]
    fwd<-readRDS(mr_ts_forward_resname)$resdf[[1]]
    all_results<-full_join(bckwd, fwd)
    all_results$exposure<-gsub('exposure', gwasname, all_results$exposure)
    
    all_results %>% 
      write_tsv(mr_ts_resdfname)
  }
  if (MR.FORWARD==TRUE && MR.BACKWARD==FALSE){
    fwd<-readRDS(mr_ts_forward_resname)$resdf[[1]]
    fwd$exposure<-gsub('exposure', gwasname, fwd$exposure)
    fwd %>% 
      write_tsv(mr_ts_resdfname)
  }
  if (MR.FORWARD==FALSE && MR.BACKWARD==TRUE){
    bckwd<-readRDS(mr_ts_backward_resname)$resdf[[1]]
    bckwd %>% 
      write_tsv(mr_ts_resdfname)
  }
}

if (!file.exists(mr_ts_signif_resdfname)){
  print('Saving significant MR results...')
  mrres<-read_tsv(mr_ts_resdfname)
  mrres_signif<-mrres %>% 
    filter(pval<MR.RESULTS.PVAL_THRESHOLD) %>% 
    filter(method=='Inverse variance weighted')
  mrres_signif %>% 
    write_tsv(mr_ts_signif_resdfname)
}


#-----***-----
#-----3. Visualization of MR results-----
if (MR.VISUALIZE==TRUE){
  if (!file.exists(mr_ts_graphdata)){
  print('Plotting the results...')
  mrres<-read_tsv(mr_ts_resdfname)
  mrres_signif<-read_tsv(mr_ts_signif_resdfname)
  
  nodes<-c()
  nodes<-c(nodes, mrres_signif$exposure)
  nodes<-c(nodes, mrres_signif$outcome)
  nodes<-unique(nodes)
  other_nodes=nodes[nodes!=gwasname]
  
  mr_data_forvis<-mrres %>% 
    filter(method=='Inverse variance weighted') %>% 
    filter((exposure %in% other_nodes) | (outcome %in% other_nodes))
  mr_data_forvis<-mr_data_forvis %>% 
    mutate(is_significant_interaction=as.numeric(pval<MR.RESULTS.PVAL_THRESHOLD)) %>% 
    mutate(neglgp=-log10(pval)) %>% 
    mutate(increase=as.numeric(or>1)+1)
  mr_data_forvis_graph<-mr_data_forvis %>% 
    select(-c(id.exposure, id.outcome, method, exposure.pipename, outcome.pipename)) %>% 
    mutate(color_code=increase*is_significant_interaction)
  
  ao_ann<-available_outcomes()
  ao_ann<-ao_ann%>% 
    mutate(vertex=paste(trait, '||', paste('id:',id, sep='')))
    
  nodedata<-tibble(vertex=nodes)
  nodedata<-left_join(nodedata, ao_ann)
  nodedata$subcategory[is.na(nodedata$subcategory)]<-'Target'
  nodedata$trait[is.na(nodedata$trait)]<-gwasname
  nodedata$subcategory[nodedata$subcategory=='NA']<-'Other'
  nodedata$id[nodedata$vertex==gwasname]<-gwassource
  
  #----->Manual subcategory setting----
  #Here subcategories can be written manually
  #nodedata$subcategory<-c('Psychiatric / Neurological', 'Psychiatric / Neurological', 'Diabetes','Diabetes','Personality','Target','Psychiatric / Neurological','Psychiatric / Neurological', 'Psychiatric / Neurological','Psychiatric / Neurological','Psychiatric / Neurological','Psychiatric / Neurological','Lipids','Lipids','Lipids','Respiratory system', 'Respiratory system','Psychiatric / Neurological','Psychiatric / Neurological','Inflammatory GI diseases','Inflammatory GI diseases','Inflammatory GI diseases', 'Cancer','Cancer','Cancer')
  mr_data_forvis_graph <- mr_data_forvis_graph[, c("exposure", "outcome", colnames(mr_data_forvis_graph)[c(3:length(colnames(mr_data_forvis_graph)))])]
  
  
  g = graph_from_data_frame(mr_data_forvis_graph, directed = TRUE, vertices = nodedata)
  vertex_groups<-V(g)[c(1:length(V(g)))]$subcategory
  vertex_groups<-factor(vertex_groups)
  vertex_groups_numeric<-unclass(vertex_groups)
  
  
  node_colors<-pal_aaas(palette = c("default"))(max(vertex_groups_numeric))
  edge_colors<-c('slategray','midnightblue','darkorchid')[E(g)$color_code+1]
  V(g)$labels<-paste(paste(V(g)$trait, sep=''),'\n', V(g)$id)
  V(g)$size<-as.numeric(V(g)$name==gwasname)*15+5
  g<-setAttribute(g, 'subcategory', V(g)$subcategory)
  l = layoutVertexByAttr(g, "subcategory", cluster.strength=7)
  
  
  
  pdf(file=paste(gwasname,"_mrplot.pdf"))
  plot(g,vertex.size=V(g)$size,
       vertex.label=V(g)$labels, 
       vertex.frame.color = "white",
       vertex.label.dist=1.2,
       vertex.label.color='gray10',
       vertex.color=adjustcolor(node_colors[vertex_groups_numeric], alpha.f = .9),
       vertex.label.family="Helvetica",
       vertex.label.cex=0.5, 
       edge.curved=TRUE,
       edge.curved=0.05,
       edge.arrow.size=1,
       edge.arrow.width=0.6,
       edge.color=adjustcolor(edge_colors, alpha.f=.7), 
       edge.width=as.integer(cut(abs(E(g)$neglgp), breaks = 5))*1.3, 
       ylim=c(-1,1),xlim=c(-1,1.2), asp = .8, 
       layout=l)
  dev.off()

  }
}

