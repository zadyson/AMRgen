# AMRgen

R Package for Genetic and Phenotypic Resistance Interpretation

## Overview

**AMRgen** is an open-source R package designed to bridge the gap between genotypic and phenotypic antimicrobial resistance (AMR) data. Developed as an extension to the [AMR R package](https://github.com/msberends/AMR), it provides tools to interpret AMR genes, integrate these findings with antimicrobial susceptibility test (AST) data, and calculate genotype-phenotype associations.

This package is developed in collaboration with the ESGEM-AMR Working Group and is tailored for researchers and healthcare professionals tackling AMR globally.

------------------------------------------------------------------------

## Key Features

> NOTE: These features are *intended* to develop in the near future

-   **Genotype-Phenotype Integration**: Links AMR gene presence with phenotypic resistance profiles, enabling deeper insights into resistance mechanisms.
-   **Automated EUCAST MIC Distribution Integration**: Fetch MIC distribution data directly from [EUCAST](https://mic.eucast.org) for seamless comparison with local susceptibility data.
-   **Visualisation**: Generate powerful UpSet plots to identify intersections of AMR gene presence and phenotypic resistance, highlighting multidrug resistance patterns.
-   **NCBI-Compliant Export**: Create NCBI-compliant antibiograms for global data sharing and interoperability with surveillance platforms.
-   **Modular and Extensible**: Leverages the robust foundation of the AMR package, including antibiotic selectors and clinical breakpoint interpretations.

------------------------------------------------------------------------

## Getting Started

To install and explore the package, follow the instructions below:

### Installation

``` r
# Install from GitHub
remotes::install_github("interpretAMR/AMRgen")
```

## Usage Examples

_Will follow shortly_

## Contributions

Contributions are welcome! If you encounter issues or wish to suggest new features, please open an issue or submit a pull request.

## Licence

This package is distributed under the GNU GPL-3.0 Licence. See `LICENSE` for details.
