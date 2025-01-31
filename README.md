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
``` r
# Install from GitHub
remotes::install_github("interpretAMR/AMRgen")
```

## Usage Examples

```
library(AMRgen)
```

### Import pheno data (from NCBI AST) and geno data (AMRfinderplus output), and compare geno/pheno for drugs of interest

```
# import small example E. coli AST data from NCBI (without re-interpreting resistance)
pheno <- import_ncbi_ast("testdata/Ecoli_AST_NCBI_n50.tsv")

# import small example E. coli AMRfinderplus data
geno <- import_amrfp("testdata/Ecoli_AMRfinderplus_n50.tsv", "Name")

# get matrix combining data on phenotype (binary R and NWT for a drug) and genotype (binary presence/absence for markers for relevant drug class/es)

# 1: meropenem phenotype vs carbapenem/cephalosporin genetic markers
mero_vs_blaGenes <- get_binary_matrix(geno, pheno, antibiotic="Meropenem", drug_class_list=c("Carbapenems", "Cephalosporins"), sir_col="Resistance phenotype")

# 2: ciprofloxacin phenotype vs quinolone genetic markers
cipro_vs_quinoloneMarkers <- get_binary_matrix(geno, pheno, "Ciprofloxacin", c("Quinolones"), sir_col="Resistance phenotype")
```


### Import small example E. coli AST data from NCBI, and re-interpret resistance and ECOFF using AMR package

_**(WARNING: phenotype interpretation can take a few minutes)**_

```
pheno <- import_ncbi_ast("testdata/Ecoli_AST_NCBI_n50.tsv", interpret = T, ecoff=T)
```

### Investigate ciprofloxacin resistance vs quinolone genotype markers, via solo PPV and upset plots
_(Using the resistance interpretations as imported from NCBI AST file, not reinterpreted from assay values)_

```
# import larger example E. coli AST data from NCBI (without re-interpreting resistance)
pheno <- import_ncbi_ast("testdata/Ecoli_AST_NCBI.tsv.gz")

# import larger example E. coli AMRfinderplus data
geno <- import_amrfp("testdata/Ecoli_AMRfinderplus.tsv.gz", "Name")

# find genomes with just one quinolone resistance marker, then estimate and plot positive predictive value (PPV) for ciprofloxacin resistance/NWT
soloPPV_cipro <- solo_ppv_analysis(geno, pheno, antibiotic="Ciprofloxacin", drug_class_list=c("Quinolones"), sir_col="Resistance phenotype")

# view PPV summary statistics:
soloPPV_cipro$solo_stats

# plot PPV summary:
soloPPV_cipro$combined_plot

# get matrix combining data on ciprofloxacin phenotype (MIC, plus binary R and NWT) and genotype (binary presence/absence for quinolone resistance markers)
cip_bin<- get_binary_matrix(geno, pheno, antibiotic="Ciprofloxacin", drug_class_list=c("Quinolones"), sir_col="Resistance phenotype", keep_assay_values=T, keep_assay_values_from = "mic")

# do upset plot of MIC vs genotype marker combinations (using AMRgen function, not requiring complexUpset)
amr_upset(cip_bin, min_set_size=2, order="mic")

# do upset plot of MIC vs genotype marker combinations (using complexUpset)
amr_complexUpset(cip_bin)
```

### Investigate ciprofloxacin resistance vs quinolone genotype markers, via solo PPV and upset plots
_(Using phenotypes reinterpreted from assay values in the NCBI AST file; note these have been pre-computed to save time for this example)_

```
# import E. coli AST data from NCBI (R and NWT variables have been pre-computed from the raw NCBI AST file using: import_ncbi_ast("testdata/Ecoli_AST_NCBI.tsv", interpret=T, ecoff=T)
pheno <- read_tsv("testdata/Ecoli_AST_NCBI_reinterpreted.tsv.gz") %>%
    mutate(drug_agent=as.ab(drug_agent), spp_pheno=as.mo(spp_pheno), mic=as.mic(mic), disk=as.disk(disk), pheno=as.sir(pheno))

# import AMRfinderplus data
geno <- import_amrfp("testdata/Ecoli_AMRfinderplus.tsv.gz", "Name")

# do positive predictive value (PPV) analysis for quinolone resistance markers vs ciprofloxacin resistance/NWT
soloPPV_cipro <- solo_ppv_analysis(geno, pheno, antibiotic="Ciprofloxacin", drug_class_list=c("Quinolones"), sir_col="Resistance phenotype")

soloPPV_cipro$solo_stats

soloPPV_cipro$combined_plot

# get matrix combining data on ciprofloxacin phenotype (MIC, plus binary R and NWT) and genotype (binary presence/absence for quinolone resistance markers)
cip_bin<- get_binary_matrix(geno, pheno, antibiotic="Ciprofloxacin", drug_class_list=c("Quinolones"), sir_col="Resistance phenotype", keep_assay_values=T, keep_assay_values_from = "mic")

# do upset plot of MIC vs genotype marker combinations (using AMRgen function, not requiring complexUpset)
amr_upset(cip_bin, min_set_size=2, order="mic")

# do upset plot of MIC vs genotype marker combinations (using complexUpset)
amr_complexUpset(cip_bin)
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
```
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
