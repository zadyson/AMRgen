#' E. coli NCBI AST example data
#'
#' A subset of E. coli phenotype data from the NCBI AST browser
#'
#' @format ## `ecoli_ast_raw`
#' A data frame with 10 rows and 17 columns representing unprocessed data from NCBI AST Browser:
#' \describe{
#'   \item{`#BioSample`}{sample identifier}
#'   \item{`Scientific name`}{species identifier}
#'   \item{Antibiotic}{antibiotic name}
#'   \item{`Testing standard`}{interpretation standard (EUCAST or CLSI)}
#'   \item{`Measurement sign`}{measurement sign (>, <, =, etc) relating to MIC measurement}
#'   \item{`MIC (mg/L)`}{minimum inhibitory concentration}
#'   \item{`Disk diffusion (mm)`}{disk diffusion zone}
#'   \item{`Resistance phenotype`}{resistance call (SIR) as submitted}
#'   
#'   ...
#' }
#' @source <https://www.ncbi.nlm.nih.gov/pathogens/ast>
"ecoli_ast_raw"

#' E. coli NCBI AST example data, re-interpreted with AMR package
#'
#' A subset of E. coli phenotype data from the NCBI AST browser
#'
#' @format ## `ecoli_ast`
#' A data frame with 15320 rows and 24 columns representing reinterpreted data from NCBI AST Browser:
#' \describe{
#'   \item{`id`}{sample identifier, imported from BioSample column in raw input}
#'   \item{spp_pheno}{species identifier, interpreted from `Scientific name` using as.mo, used to interpret ecoff and pheno columns}
#'   \item{drug_agent}{antibiotic code, interpreted from `Antibiotic` using as.ab, used to interpret ecoff and pheno columns}
#'   \item{ecoff}{WT/NWT classification, interpreted using as.sir}
#'   \item{pheno}{S/I/R classification, interpreted using as.sir}
#'   \item{mic}{minimum inhibitory concentration, formated using as.mic, used to interpret ecoff and pheno columns}
#'   \item{disk}{disk diffusion zone, formated using as.disk, used to interpret ecoff and pheno columns}
#'   \item{guideline}{interpretation guidelines used to interpret ecoff and pheno columns}
#'   \item{`Scientific name`}{species identifier, from input file}
#'   \item{Antibiotic}{antibiotic name, from input file}
#'   \item{`Testing standard`}{interpretation standard (EUCAST or CLSI), from input file}
#'   \item{`Measurement sign`}{measurement sign (>, <, =, etc) relating to MIC measurement, from input file}
#'   \item{`MIC (mg/L)`}{minimum inhibitory concentration, from input file}
#'   \item{`Disk diffusion (mm)`}{disk diffusion zone, from input file}
#'   \item{`Resistance phenotype`}{resistance call (SIR), from input file}
#'   
#'   ...
#' }
#' @source <https://www.ncbi.nlm.nih.gov/pathogens/ast>
"ecoli_ast"

#' E. coli genotype example data
#'
#' Genotypes called using AMRfinderplus (v3.12.8, DB 2024-01-31.1), sourced from AllTheBacteria project
#'
#' @format ## `ecoli_geno_raw`
#' A data frame with 53591 rows and 10 columns representing genotyping results from AMRfinderplus:
#' \describe{
#'   \item{Name}{sample identifier}
#'   \item{`Gene symbol`}{gene symbol in NCBI refgene}
#'   \item{`Hierarchy node`}{node in NCBI hierarchy}
#'   \item{Class, Subclass}{drug class/es associated with marker (in NCBI refgene)}
#'   \item{`% Coverage of reference sequence`,% Identity to reference sequence`, `Accession of closest sequence` }{sequence match information}
#'   
#'   ...
#' }
#' @source <https://github.com/ncbi/amr/wiki/Interpreting-results>
"ecoli_geno_raw"