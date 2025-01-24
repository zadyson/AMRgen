#input = pheno data frame, name of the ID column (otherwise take first column, or try to find 'sample' or 'biosample'; same for geno data frame
#compare the sample IDs, make a list of unique entries that appear in both
#report the number of overlaps and uniques
#return copies of geno and pheno data frames that each contain the overlapping samples only

compare_geno_pheno_id <- function(geno_data, pheno_data, geno_sample_col=NULL, pheno_sample_col=NULL, rename_id_cols=F){
  if (is.null(geno_sample_col)){
    geno_sample_col <- colnames(geno_data)[1]
  }
  if(is.null(pheno_sample_col)){
    pheno_sample_col <- colnames(pheno_data)[1]
  }

  # get sample ids from both files, extract as vector
  geno_sample_ids <- unique(pull(geno_data, geno_sample_col))
  pheno_sample_ids <- unique(pull(pheno_data, pheno_sample_col))

  # compare the two 
  # Unique values in pheno and geno
  unique_in_pheno <- setdiff(pheno_sample_ids, geno_sample_ids)
  unique_in_geno <- setdiff(geno_sample_ids, pheno_sample_ids)

  # Overlapping values
  overlapping_values <- intersect(pheno_sample_ids, geno_sample_ids)
  
  geno_matched <- geno_data %>% filter(get(geno_sample_col) %in% overlapping_values)
  pheno_matched <- pheno_data %>% filter(get(pheno_sample_col) %in% overlapping_values)
  
  if (rename_id_cols) {
    geno_matched <- geno_matched %>% rename(id=any_of(geno_sample_col))
    pheno_matched <- pheno_matched %>% rename(id=any_of(pheno_sample_col))
  }

  return(list(pheno_unique = unique_in_pheno, 
    geno_unique = unique_in_geno, 
    overlap_ids = overlapping_values,
    geno_matched = geno_matched,
    pheno_matched = pheno_matched))

}