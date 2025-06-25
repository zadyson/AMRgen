# subset of NCBI AST data that has already been re-interpreted using import_ncbi_ast
# provided for testing solo_ppv_analysis and amr_upset functions
ecoli_ast <- read_tsv("ecoli_pheno.tsv.gz") %>%
  import_ncbi_ast() %>%
  mutate(drug_agent = as.ab(drug_agent)) %>%
  mutate(mic = as.mic(mic)) %>%
  mutate(disk = as.disk(disk)) %>%
  mutate(spp_pheno = as.mo(spp_pheno)) %>%
  mutate(pheno = as.sir(mic, ab = drug_agent, mo = spp_pheno, guideline = "CLSI 2025")) %>%
  mutate(pheno = as.sir(pheno)) %>% select(id, spp_pheno, pheno, drug_agent, mic, disk, guideline)

# small dataframe that mimics a raw import of NCBI AST tab-delim file
# provided for testing import_ncbi_ast, including re-interpreting (so only 10 rows)
ecoli_ast_raw <- ecoli_ast %>%
  select(-c(spp_pheno:guideline)) %>%
  rename(`#BioSample` = id) %>%
  head(n = 10)

# set of AMRfinderplus (v3.12.8, DB 2024-01-31.1) genotype results sourced from AllTheBacteria
# for the same biosamples as ecoli_ast
# provided for testing import_amrfp, solo_ppv_analysis and amr_upset functions
ecoli_geno_raw <- read_tsv("ecoli_geno.tsv.gz")

usethis::use_data(ecoli_ast, internal = F, overwrite = TRUE)
usethis::use_data(ecoli_ast_raw, internal = F, overwrite = TRUE)
usethis::use_data(ecoli_geno_raw, internal = F, overwrite = TRUE)
