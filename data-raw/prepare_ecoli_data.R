# subset of NCBI AST data that has already been re-interpreted using import_ncbi_ast
# provided for testing solo_ppv_analysis and amr_upset functions
ecoli_ast <- read_tsv("ecoli_pheno.tsv.gz") %>%
  import_ncbi_ast(interpret = TRUE, ecoff = TRUE, default_guideline = "CLSI") %>%
  filter(`Scientific name` = "Escherichia coli") %>%
  filter(pheno != "NI")

# small dataframe that mimics a raw import of NCBI AST tab-delim file
# provided for testing import_ncbi_ast, including re-interpreting (so only 10 rows)
ecoli_ast_raw <- ecoli_ast %>%
  select(-c(spp_pheno:guideline)) %>%
  rename(`#BioSample` = id) %>%
  head(n = 10)

# remove unneeded column headers from AST object
ecoli_ast <- ecoli_ast %>%
  select(id, spp_pheno, pheno, drug_agent, mic, disk, guideline)

# set of AMRFinderPlus (v3.12.8, DB 2024-01-31.1) genotype results sourced from AllTheBacteria
# for the same biosamples as ecoli_ast
# provided for testing import_amrfp, solo_ppv_analysis and amr_upset functions
ecoli_geno_raw <- read_tsv("ecoli_geno.tsv.gz")

usethis::use_data(ecoli_ast, internal = FALSE, overwrite = TRUE)
usethis::use_data(ecoli_ast_raw, internal = FALSE, overwrite = TRUE)
usethis::use_data(ecoli_geno_raw, internal = FALSE, overwrite = TRUE)
