getBinMatMIC <- function(pheno_table, geno_table, antibiotic, drug_class,
                         geno_sample_col=NULL, pheno_sample_col=NULL,
                         plot_cols=c("resistant"="IndianRed", 
                                     "intermediate"="orange", 
                                     "nonsusceptible"="orange", 
                                     "susceptible"="lightgrey")) {
  
  # subset pheno & geno dataframes to those samples with overlap
  overlap <- compare_geno_pheno_id(pheno_table, geno_table)
  pheno_matched <- overlap$pheno_matched
  geno_matched <- overlap$geno_matched
    
  # extract list of relevant drug markers
  
  # for now, user must specify the 'drug_class' as a parameter and we assume genotype_table has a 'drug_class'
  # to do: determine the appropriate drug_class value to match to, from the 'antibiotic' parameter; 
  #        and same if geno_table has a drug_agent entry with no drug_class
  #if(is.null(drug_class)) <- getDrugClass(antibiotic)
  
  # extract geno rows for relevant markers
  if (!("drug_class" %in% colnames(geno_table))) {
    stop("input geno_table must have a column labelled `drug_class`")
  }
  else if (!(drug_class %in% geno_matched$drug_class)) {
    stop(paste("No markers in matching drug class", drug_class, "were identified in input geno_table"))
  }
  else {
    markers <- geno_matched %>% filter(drug_class %in% drug_class) %>% pull(marker)
  }
  
  geno_binary_markers <- geno_matched %>% 
    filter(marker %in% markers) %>% 
    group_by(geno_id_col, marker) %>% count() %>%
    ungroup() %>%
    right_join(pheno_matched %>% filter(Antibiotic==antibiotic) %>% select(BioSample), 
               join_by("Name"=="BioSample")) %>%
    distinct() %>%
    pivot_wider(id_cols=Name, names_from=`Gene symbol`, values_from=n, values_fill=0) 
  
  if ("NA" %in% colnames(afp_binary_markers)) {
    afp_binary_markers <- afp_binary_markers %>%select(-c("NA"))
  }
  
  afp_binary_markers <- afp_binary_markers %>% 
    select(Name, sort(names(.))) %>%
    inner_join(pheno_matched %>% filter(Antibiotic==antibiotic) %>% select(BioSample, `Resistance phenotype`, `MIC (mg/L)`, `Disk diffusion (mm)`), 
               join_by("Name"=="BioSample")) %>%
    mutate(resistant=case_when(`Resistance phenotype`=="resistant" ~ 1,
                               `Resistance phenotype` %in% c("intermediate", "susceptible") ~ 0,
                               TRUE ~ NA)) %>%
    select(-`Resistance phenotype`, -Name)
  
  return(afp_binary_markers)
  
}
