# subset of NCBI AST data that has already been re-interpreted using import_ncbi_ast
# provided for testing solo_ppv_analysis and amr_upset functions
ecoli_ast <- read_tsv("ecoli_pheno.tsv.gz")

# small dataframe that mimics a raw import of NCBI AST tab-delim file
# provided for testing import_ncbi_ast, including re-interpreting (so only 10 rows)
ecoli_ast_raw <- ecoli_ast %>% select(-c(spp_pheno:guideline)) %>% rename(`#BioSample`=id) %>% head(n=10)

# set of AMRfinderplus (v3.12.8, DB 2024-01-31.1) genotype results sourced from AllTheBacteria
# for the same biosamples as ecoli_ast
# provided for testing solo_ppv_analysis and amr_upset functions
ecoli_geno <- read_tsv("ecoli_geno.tsv.gz")

# dataframe that mimics a raw import of AMRfinderplus tab-delim file
# provided for testing import_amrfp
ecoli_geno_raw <- ecoli_geno %>% select(-c(marker:drug_class))

usethis::use_data(ecoli_ast, internal = TRUE)
usethis::use_data(ecoli_ast_raw, internal = TRUE)
usethis::use_data(ecoli_geno, internal = TRUE)
usethis::use_data(ecoli_geno_raw, internal = TRUE)
