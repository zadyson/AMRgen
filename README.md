# AMRgen

R Package for Genetic and Phenotypic Resistance Interpretation 

## Overview

**AMRgen** is an open-source R package designed to bridge the gap between genotypic and phenotypic antimicrobial resistance (AMR) data. Developed as an extension to the [AMR R package](https://github.com/msberends/AMR), it provides tools to interpret AMR genes, integrate these findings with antimicrobial susceptibility test (AST) data, and calculate genotype-phenotype associations.

This package is developed in collaboration with the ESGEM-AMR Working Group and is tailored for researchers and healthcare professionals tackling AMR globally.

------------------------------------------------------------------------

## Key Features

-   **Import Genotype and Phenotype Data**: Import from common formats (NCBI AST for phenotypes; AMRfinderplus and hAMRonization for genotypes
-   **Genotype-Phenotype Integration**: Links AMR gene presence with phenotypic resistance profiles, enabling deeper insights into resistance mechanisms.
-   **Automated EUCAST MIC Distribution Integration**: Fetch MIC distribution data directly from [EUCAST](https://mic.eucast.org) for seamless comparison with local susceptibility data.
-   **Visualisation**: Generate powerful UpSet plots to identify intersections of AMR gene presence and phenotypic resistance, highlighting multidrug resistance patterns.
-   **Modular and Extensible**: Leverages the robust foundation of the AMR package, including antibiotic selectors and clinical breakpoint interpretations.

> Planned for development
-   **NCBI-Compliant Export**: Export phenotype data to NCBI-compliant antibiogram format.
-   **Expanded Data Import**: Import and parse phenotype data from other tools (e.g. CARD, ResFinder).

------------------------------------------------------------------------

## Getting Started

To install and explore the package, follow the instructions below:

### Installation
Note that this package requires the latest version of the `AMR` package (still in beta).

Install the latest version of the `AMR` package with:
```r
install.packages("remotes") # if you haven't already
remotes::install_github("msberends/AMR")
```

Then install this package
```r
# Install from GitHub
remotes::install_github("interpretAMR/AMRgen")
```

It is best to restart R before running the installation. If you didn't do this and/or you encounter issues with the examples below after install, it may help to also restart after the install and start fresh with the examples below.

## Usage Examples

```r
library(AMRgen)
library(tidyverse)
```

### Import pheno data (from NCBI AST) and geno data (AMRfinderplus output), and compare geno/pheno for drugs of interest

```r
# small example E. coli AST data from NCBI
ecoli_ast_raw

# import without re-interpreting resistance
pheno <- import_ncbi_ast(ecoli_ast_raw)

# small example E. coli AMRfinderplus data
ecoli_geno_raw

# import AMRfinderplus data
geno <- import_amrfp(ecoli_geno_raw %>% head(n=10), "Name")

```


### Import small example E. coli AST data from NCBI, and re-interpret resistance and ECOFF using AMR package

_**(WARNING: phenotype interpretation can take a few minutes)**_

```r
pheno <- import_ncbi_ast(ecoli_ast_raw, interpret = T, ecoff=T)
```

### Investigate ciprofloxacin resistance vs quinolone genotype markers, via solo PPV and upset plots

```r
# larger example E. coli AST data from NCBI, already imported (with import_ncbi_ast) and re-interpreted (with as.sir)
ecoli_ast

# import matching AMRfinderplus data
ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")

# find genomes with just one quinolone resistance marker, then estimate and plot positive predictive value (PPV) for ciprofloxacin resistance/NWT
# (using phenotypes interpreted with AMR package; alternatively set sir_col="Resistance phenotype" to use the classifications from the raw NCBI AST file)
soloPPV_cipro <- solo_ppv_analysis(ecoli_geno, ecoli_ast, antibiotic="Ciprofloxacin", drug_class_list=c("Quinolones"), sir_col="pheno")

# view PPV summary statistics:
soloPPV_cipro$solo_stats

# get matrix combining data on ciprofloxacin phenotype (MIC, plus binary R and NWT) and genotype (binary presence/absence for quinolone resistance markers)
cip_bin<- get_binary_matrix(ecoli_geno, ecoli_ast, antibiotic="Ciprofloxacin", drug_class_list=c("Quinolones"), sir_col="pheno", keep_assay_values=T, keep_assay_values_from = "mic")

# do upset plot of MIC vs genotype marker combinations (using complexUpset)
amr_complexUpset(cip_bin)

# do upset plot of MIC vs genotype marker combinations, for combinations observed at least 5 times (using AMRgen function, not requiring complexUpset)
cip_upset <- amr_upset(cip_bin, min_set_size=5, order="mic")

# check the summary information for each marker combination (number of genomes, MIC median and range, number and fraction resistant)
cip_upset$summary

```


### Explore logistic regression models of genotype vs phenotype

```r

# get binary matrix
ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")
cip_bin<- get_binary_matrix(ecoli_geno, ecoli_ast, antibiotic="Ciprofloxacin", drug_class_list=c("Quinolones"), sir_col="pheno")

# logistic regression model for ciprofloxacin resistance (R vs S/I), predicted from all detected markers associated with quinolone resistance
model <- logistf::logistf(R ~ ., data=cip_bin %>% select(-c(id,pheno,NWT)))
model_summary <- logistf_details(model)
plot_estimates(model_summary)

# include only markers observed in at least 10 samples
model <- logistf::logistf(R ~ ., data=cip_bin %>% select(-c(id,pheno,NWT)) %>% select_if(funs(sum(.)>10)))
model_summary <- logistf_details(model)
plot_estimates(model_summary, title="Logistic regression on Cipro R")

# predict NWT (defined by ECOFF) rather than R
model_NWT <- logistf::logistf(NWT ~ ., data=cip_bin %>% select(-c(id,pheno,R)) %>% select_if(funs(sum(.)>10)))
model_NWT_summary <- logistf_details(model_NWT)
plot_estimates(model_NWT_summary, title="Logistic regression on Cipro NWT")

# compare estimates for R and NWT (on a single plot)
compare_estimates(model_summary, model_NWT_summary, single_plot = T, title1="R", title2="NWT", title="R and NWT for Cipro")

# compare estimates for R and NWT (two plots, side-by-side)
compare_estimates(model_summary, model_NWT_summary, single_plot = F, title1="R", title2="NWT", title="R and NWT for Cipro")

# organise layout using patchwork
library(patchwork)
compare_estimates(model_summary, model_NWT_summary, single_plot = F, title1="R", title2="NWT", title="R and NWT for Cipro") + plot_layout(guides="collect", axes="collect")

# do logistic regression analysis for ciprofloxacin R and NWT vs quinolone associated markers present in ≥50 samples and do combined plot
logistic_cipro <- amr_logistic(ecoli_geno, ecoli_ast, "Ciprofloxacin", c("Quinolones"), maf=5)

logistic_cipro$plot
logistic_cipro$bin_mat
logistic_cipro$modelR
logistic_cipro$modelNWT

```
### Explore evidence for geno/pheno association from solo PPV, logistic regression, and upset plots

```r
ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")

soloPPV_cipro <- solo_ppv_analysis(ecoli_geno, ecoli_ast, antibiotic="Ciprofloxacin", drug_class_list=c("Quinolones"), sir_col="pheno")

logistic_cipro <- amr_logistic(ecoli_geno, ecoli_ast, "Ciprofloxacin", c("Quinolones"), maf=5)

# combine solo PPV and logistic regression coefficients
allstatsR <- merge_logreg_soloppv(logistic_cipro$modelR, soloPPV_cipro$solo_stats %>% filter(category=="R"), title="Quinolone markers vs Cip R")
allstatsNWT <- merge_logreg_soloppv(logistic_cipro$modelNWT, soloPPV_cipro$solo_stats %>% filter(category=="NWT"), title="Quinolone markers vs Cip NWT")

allstatsR$plot + allstatsNWT$plot + plot_layout(guides="collect")

# plots show no markers have solo PPV >50% for R, but many have solo PPV >70% for NWT
# all these are also significant in logistic regression models 
# which markers are predictive?
allstatsNWT$combined %>% filter(ppv>0.5) %>% arrange(-ppv)

# which markers are significant in regression but have low solo PPV (or no solo genomes to provide evidence)?
allstatsNWT$combined %>% filter(pval<0.05) %>% arrange(-ppv)

# do upset plot of MIC vs genotype marker combinations, ordering combinations by mean MIC
# shows MIC distribution for solo markers is lower than combinations of markers
amr_complexUpset(soloPPV_cipro$amr_binary, min_set_size=5, intersect_counts=F)

# using the other upset function
amr_upset(soloPPV_cipro$amr_binary, min_set_size=2, order="mic")


```

### Download and plot reference MIC distribution from eucast.org

```r
# get MIC distribution for ciprofloxacin, for all organisms
get_eucast_mic_distribution("cipro")

# specify microorganism to only get results for that pathogen
kleb_cip_mic_data <- get_eucast_mic_distribution("cipro", "K. pneumoniae")

# get disk diffusion data instead
kleb_cip_disk_data <- get_eucast_disk_distribution("cipro", "K. pneumoniae")

# plot the MIC data 
mics <- rep(kleb_cip_mic_data$mic, kleb_cip_mic_data$count)
ggplot2::autoplot(mics, ab = "cipro", mo = "K. pneumoniae", title = "Look at my MICs!")

```


### Compare own MIC data vs reference distribution from EUCAST

```r

# random test
my_mic_values <- AMR::random_mic(500)
comparison <- compare_mic_with_eucast(my_mic_values, ab = "cipro", mo = "K. pneumoniae")
comparison
ggplot2::autoplot(comparison)


# example E. coli data
comparison <- compare_mic_with_eucast(ecoli_ast %>% filter(drug_agent=="CIP") %>% pull(mic), ab = "cipro", mo = "E. coli")
comparison
ggplot2::autoplot(comparison) + ggtitle("E. coli - Ciprofloxacin")
```

## Contributions

Contributions are welcome! If you encounter issues or wish to suggest new features, please open an issue or submit a pull request.

## Licence

This package is distributed under the GNU GPL-3.0 Licence. See `LICENSE` for details.



## Dev notes:

Comnon formats generated by import functions, and expected by analysis/plotting functions

### Definition: 'genotype' dataframe

required fields:
- a column indicating the sample ID (any string; to be matched to phenotype data)
- a column named 'marker' indicating the label for a specific gene or mutation (S3 class 'gene')

at least one of:
- a column named 'drug_agent' (S3 class 'ab'; this will usually be generated by a genotype parser function applying as.ab)
- a column named 'drug_class' indicating a drug class associated with this marker (controlled vocab string; allowed values are those in antbiotics$group, this will usually be generated by a genotype parser function, if not provided it will be generated from 'drug_agent') [might develop S3 class? e.g. needs to include 'efflux']

optionally:
- a column indicating the species (S3 class mo; to facilitate interpretation)
  
genotype parsers should be generating all 4 fields

### Definition: 'phenotype' dataframe

required fields:
- a column indicating the sample ID (any string; to be matched to genotype data)
- a column named 'drug_agent' (S3 class 'ab'; this will usually be generated by a phenotype parser function applying as.ab)

at least one of:
- a column of class 'disk' or 'mic'
- a column of class 'sir'

optionally:
- a column indicating the species (S3 class mo; to facilitate interpretation)

### Expected workflow (target for dev)

* import genotype data -> genotype dataframe (e.g. `import_amrfp`)
* import phenotype data -> phenotype dataframe (e.g. `import_ncbi_ast`)
  - interpret SIR if required (as.sir; requires either a species column, or that all rows are a single species)
* optionally: filter both files to the desired sample sets (e.g. filter on species, check common sample identifiers exist)
* pass filtered genotype & phenotype objects (which have common sample identifiers) to functions for
  - generating binary matrix of SIR vs marker presence/absence suitable for regression modelling (`get_binary_matrix`)
  - cross-tabulating SIR vs marker presence/absence, calculating & plotting PPV (`solo_ppv_analysis`)
  - upset plots showing MIC/DD distribution stratified by genotype profile (`amr_complexUpset` or `amr_upset`)

### Code for testing harmonize_data - but note the required input files are not in this repo

```r
# test code amrfinder plus
# note both amrfinderplus test files appear malformed according to hamronize 
# and produce errors but the code works
user_software_name <- "amrfinderplus"
user_software_version <- "3.12.8"
user_input_filename <- "/Users/lshzd1/Desktop/ATB_Achromobacter_AFP.tsv"
user_database_version <- "2024-01-31.1"

test_data <- harmonize_data(user_software_name, user_software_version, 
                            user_database_version, user_input_filename)


# test code rgi - same arguments as for amrfinder plus
user_software_name <- "rgi"
user_software_version <- "version x"
user_input_filename <- "/Users/lshzd1/Desktop/2025-01-14_11:06:52.908_KPN2009.fasta.txt"
user_database_version <- "database y"

test_data <- harmonize_data(user_software_name, user_software_version, 
                            user_database_version, user_input_filename)


# test code - resfinder - must be json (hamronize can't do txt file)
user_software_name <- "resfinder"
user_software_version <- "4.6.0"
user_input_filename <- "/Users/lshzd1/Desktop/KPN2214.json"
user_database_version <- "2024-08-06"

test_data <- harmonize_data(user_software_name, user_software_version, 
                            user_database_version, user_input_filename)
```
