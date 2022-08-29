library(tidyverse)
library(dplyr)
library(TwoSampleMR)
library(readr)
library(igraph)
library(ggsci)
library(NetPathMiner)

#----FUNCTIONALITY----
#Adds rsid to summary statistics if a mapping file with chromosomal positions is supplied
#Performs bulk MR analysis in both directions with specific summary statistics from TwoSampleMR package or with summary statistics for phenotypes containing specific substrings in their names

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




#-----SETTINGS-----
gwas_sumstatfile_path<-'/home/biorp/Gitrepos/Psychiatry/SUMSTATS/anged/anged/dop_pheno_all.anged.glm.logistic.hybrid'
annotationfile_path<-'/home/biorp/Gitrepos/Psychiatry/SUMSTATS/psychiatric_genomics_chip_mapping.tsv' #file with mapping columns: CHR, POS, rsid
ldsc_munge_sumstats_path<-'/home/biorp/TOOLS/ldsc/munge_sumstats.py'
gwasname<-'Anhedonia'
gwassource<-'current study'
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






#Two Sample MR settings
mr_ts_backward_resname=paste(gwasname, '_twosample_backward.rds', sep='')
mr_ts_forward_resname=paste(gwasname, '_twosample_forward.rds', sep='')
mr_ts_backward_signif_resname=paste(gwasname, '_twosample_backward_signif.rds', sep='')
mr_ts_forward_signif_resname=paste(gwasname, '_twosample_forward_signif.rds', sep='')
mr_ts_resdfname=paste(gwasname, '_all_mr_res.tsv', sep='')
mr_ts_signif_resdfname=paste(gwasname, '_all_mr_res_signif.tsv', sep='')
mr_ts_graphdata<-paste(gwasname, '_all_mr_res_signif_graph.rds', sep='')

MR.FORWARD=TRUE
MR.BACKWARD=TRUE
ao<-available_outcomes()
MR.SAVE_ONLYSUMMARY=TRUE
MR.COMPARISON_STUDIES.MODE='PHENONAMES' #'PHENONAMES' or 'SPECIFIC'
MR.RESULTS.PVAL_THRESHOLD=0.05

MR.COMPARISON_STUDIES.RUN_PHENONAMES<-c('schiz','depress','bipol','cholester','lipo','diab','fatty','covid','respir','serotonin','cholin','dopamin','anxiety','neuro','mental','brain')
MR.COMPARISON_STUDIES.RUN_SUBCATEGORIES<-c('Fatty acid','Cofactors and vitamins','Lipid','Amino acid','Carbohydrate','Metal','Personality','Unknown metabolite','Peptide','Protein','Behavioural','Sleeping','Diabetes','Psychiatric / neurological','Autoimmune / inflammatory','Immune system','Education','Metabolites ratio','Biomarker', "Diabetes", "Keto acid", "Cancer", "Glycemic", "Immune system", "Immune cell subset frequency","Education","Hormone","Blood pressure","Anthropometric")
MR.COMPARISON_STUDIES.RUN_SPECIFIC=c() #indexes of studies from ao
MR.PHENONAMES.NCASE_THRESHOLD=4000
MR.FORWARD.THRESHOLD_PVALS=c(5e-5) #threshold p-value for instrument extraction from the studied sumstat.
MR.BACKWARD.THRESHOLD_PVALS=c(5e-8)
MR.VISUALIZE=TRUE

#---***----
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

#----***----
#----2. Mendelian randomization (two sample) run----
#----2.1 Selection of comparison studies----
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
    
    
  }
}
mrres<-read_tsv(mr_ts_resdfname)
mrres_signif<-read_tsv(mr_ts_signif_resdfname)

nodes<-c()
nodes<-c(nodes, mrres_signif$exposure)
nodes<-c(nodes, mrres_signif$outcome)
nodes<-unique(nodes)
nodes
other_nodes=nodes[nodes!='anhedonia']

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
mr_data_forvis_graph

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
nodedata$subcategory<-c('Psychiatric / Neurological', 'Psychiatric / Neurological', 'Diabetes','Diabetes','Personality','Target','Psychiatric / Neurological','Psychiatric / Neurological', 'Psychiatric / Neurological','Psychiatric / Neurological','Psychiatric / Neurological','Psychiatric / Neurological','Lipids','Lipids','Lipids','Respiratory system', 'Respiratory system','Psychiatric / Neurological','Psychiatric / Neurological','Inflammatory GI diseases','Inflammatory GI diseases','Inflammatory GI diseases', 'Cancer','Cancer','Cancer')



g = graph_from_data_frame(mr_data_forvis_graph, directed = TRUE, vertices = nodedata)
vertex_groups<-V(g)[c(1:length(V(g)))]$subcategory
vertex_groups<-factor(vertex_groups)
vertex_groups_numeric<-unclass(vertex_groups)


node_colors<-pal_aaas(palette = c("default"))(max(vertex_groups_numeric))

#GBV1 = GroupByVertex01(vertex_groups_numeric)
#GBV2 = GroupByVertex02(vertex_groups)
edge_colors<-c('slategray','blueviolet','darkred')[E(g)$color_code+1]
edge_colors<-c('slategray','midnightblue','darkorchid')[E(g)$color_code+1]

V(g)$labels<-paste(paste(V(g)$trait, sep=''),'\n', V(g)$id)

V(g)$size<-as.numeric(V(g)$name=='Anhedonia')*15+5

g<-setAttribute(g, 'subcategory', V(g)$subcategory)
l = layoutVertexByAttr(g, "subcategory", cluster.strength=7)

plot(g,vertex.size=V(g)$size,
     vertex.label=V(g)$labels, 
     vertex.frame.color = "white",
     vertex.label.dist=1.2,
     vertex.label.color='gray10',
     vertex.color=adjustcolor(node_colors[vertex_groups_numeric], alpha.f = .9),
     vertex.label.family="Helvetica",
     vertex.label.cex=0.7, 
     edge.curved=TRUE,
     edge.curved=0.05,
     edge.arrow.size=1,
     edge.arrow.width=0.6,
     edge.color=adjustcolor(edge_colors, alpha.f=.7), 
     edge.width=as.integer(cut(abs(E(g)$neglgp), breaks = 5))*1.3, 
     ylim=c(-1,1),xlim=c(-0.7,1.2), asp = 0.8, 
     layout=l)

