## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 10,
  out.width = "100%"
)

## ----setup--------------------------------------------------------------------
library(AMRgen)
library(tidyverse)

## ----import_amrfp-------------------------------------------------------------
# Example AMRfinderplus genotyping output
ecoli_geno_raw

# Load AMRfinderplus output 
#    (replace 'ecoli_geno_raw' with the filepath for any AMRfinderplus output)
ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")

# Check the format of the processed genotype table
head(ecoli_geno)

## ----import_ncbi_ast----------------------------------------------------------
# Example E. coli AST data from NCBI 
# This one has already been imported and phenotypes interpreted from assay data
# You can make your own from NCBI format data using: 
#    import_ncbi_ast("filepath/AST.tsv", interpret=T)
ecoli_ast

head(ecoli_ast)

## ----get_binary_matrix--------------------------------------------------------
# Get matrix combining phenotype data for ciprofloxacin, binary calls for R/NWT phenotype,
#    and genotype presence/absence data for all markers associated with the relevant drug 
#    class (which are labelled "Quinolones" in AMRfinderplus).
cip_bin <- get_binary_matrix(ecoli_geno, ecoli_ast, antibiotic="Ciprofloxacin", drug_class_list=c("Quinolones"), sir_col="pheno", keep_assay_values=T, keep_assay_values_from = "mic")

# check format
head(cip_bin)

# list colnames, to see full list of quinolone markers included
colnames(cip_bin)

## ----manual_logistic, fig.height=8--------------------------------------------
# Manually run Firth's logistic regression model using the binary matrix produced above
modelR <- logistf::logistf(R ~ ., data=cip_bin %>% select(-c(id,pheno,mic,NWT)))

summary(modelR)

# Extract model summary details using `logistf_details()`
modelR_summary <- logistf_details(modelR)

modelR_summary

# Plot the point estimates and 95% confidence intervals of the model
plot_estimates(modelR_summary)


## ----amr_logistic, fig.height=8-----------------------------------------------
# Alternatively, use the amr_logistic() function to model R and NWT and plot the results together
models <- amr_logistic(geno_table = import_amrfp(ecoli_geno_raw, "Name"), pheno_table = ecoli_ast, 
                       antibiotic = "Ciprofloxacin", drug_class_list = c("Quinolones"), maf=10)

# Output tables
models$modelR

models$modelNWT

# Note the matrix output is the same as cip_bin, but without the MIC data as this is not required
#    for logistic regression.
models$bin_mat

## ----solo_ppv_analysis, fig.height=8------------------------------------------
# Run a solo PPV analysis
soloPPV_cipro <- solo_ppv_analysis(ecoli_geno, ecoli_ast, antibiotic="Ciprofloxacin", drug_class_list=c("Quinolones"), sir_col="pheno")

# Output table
soloPPV_cipro$solo_stats

# Interim matrices with data used to compute stats and plots
soloPPV_cipro$solo_binary

soloPPV_cipro$amr_binary

## ----amr_upset, fig.height=8--------------------------------------------------
# Compare ciprofloxacin MIC data with quinolone marker combinations,
#    using the binary matrix we constructed earlier via get_binary_matrix()
cipro_mic_upset <- amr_upset(cip_bin, min_set_size=2, assay="mic", order="value")

# Output table
cipro_mic_upset$summary

## ----get_eucast_distribution--------------------------------------------------
# get MIC distribution for ciprofloxacin, for all organisms
get_eucast_mic_distribution("cipro")

# specify microorganism to only get results for that pathogen
ecoli_cip_mic_data <- get_eucast_mic_distribution("cipro", "E. coli")

# get disk diffusion data instead
ecoli_cip_disk_data <- get_eucast_disk_distribution("cipro", "E. coli")

# plot the MIC data 
mics <- rep(ecoli_cip_mic_data$mic, ecoli_cip_mic_data$count)
ggplot2::autoplot(mics, ab = "cipro", mo = "E. coli", title = "E. coli cipro reference distribution")


## ----compare_mic_with_eucast--------------------------------------------------
# Compare reference distribution to random test data
my_mic_values <- AMR::random_mic(500)
comparison <- compare_mic_with_eucast(my_mic_values, ab = "cipro", mo = "E. coli")
comparison
ggplot2::autoplot(comparison)


# Compare reference distribution to example E. coli data
comparison <- compare_mic_with_eucast(ecoli_ast %>% filter(drug_agent=="CIP") %>% pull(mic), ab = "cipro", mo = "E. coli")
comparison
ggplot2::autoplot(comparison) + ggtitle("E. coli - Ciprofloxacin")

