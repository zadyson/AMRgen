require(AMR)
require(tidyverse)

#' @title Gene Class and AMR Parsing Functions
#' @description Functions to create a custom "gene" class and parse AMR data.
#' @importFrom AMR as.ab
#' @importFrom tidyverse read_tsv filter add_column separate_longer_delim left_join select mutate
#' @param x A character vector to be converted to a "gene" class.
#' @param input_table A file path to the input table containing AMR data.
#' @param sample_col A character specifying the sample column in the input file.
#' @param amrfp_drugs A data frame containing AMRFinderPlus drug classes and agents. Defaults to `amrfp_drugs_table`.
#' @return For `as.gene`, an object of class "gene". For `import_amrfp`, a data frame with parsed AMR data.
#' @examples
#' \dontrun{
#' # Create a gene object
#' gene <- as.gene(c("gene1", "gene2"))
#' print(gene)
#'
#' # Parse AMR data
#' parsed_data <- import_amrfp("path/to/input_table.tsv", "SampleID")
#' }
#' @export

#' @rdname as.gene
as.gene <- function(x){
  # Check x is a character vector
  if(!is.character(x)) {
    stop("Input must be a character vector.")
  }
  # Create the object with class "gene"
  structure(x, class = c("gene", class(x)))
}

#' @rdname as.gene
#' @export
print.gene <- function(x, ...) {
  cat("An object of class 'gene':\n")
  print(unclass(x), ...)
}

amrfp_drugs_table <- read_tsv("amrfp_drug_classes_agents.tsv")

#' @rdname import_amrfp
#' @export
import_amrfp <- function(input_table, sample_col, amrfp_drugs=amrfp_drugs_table){
  # read in the table
  in_table <- read_tsv(input_table)
  # filter to only include AMR elements
  in_table_filter <- in_table %>% filter(`Element type` == "AMR")
  
  # create the gene class for the gene column
  gene_col <- as.gene(in_table_filter$`Gene symbol`)
  
  # need to do this with hamronized format, note that this will be a combo of gene symbol and then the nt/prot/aa mutation cols if it's not just gene presence/absence
  
  # Find the position of the "Gene symbol" column and insert gene after it (could replace it really)
  position <- which(names(in_table_filter) == "Gene symbol")
  in_table_gene <- add_column(in_table_filter, marker = gene_col, .after = position)
  
  # now split the Subclass column on the "/" to make them one per row, to make adding the ab names easier
  in_table_subclass_split <- in_table_gene %>% separate_longer_delim(Subclass, "/")
  
  # make two new columns - drug_class and drug_agent, where we control the vocab for the AMRFinderPlus Subclass column
  # into something that is comparable with the drugs in the AMR package
  in_table_ab <- in_table_subclass_split %>% left_join(., amrfp_drugs[,c("AFP_Subclass", "drug_agent", "drug_class")], by=c("Subclass"="AFP_Subclass")) %>%
    select(all_of(sample_col), marker, drug_agent, drug_class, everything())

  # convert drug_agent into the "ab" class (will leave NAs as is)
  in_table_ab <- in_table_ab %>% mutate(drug_agent = as.ab(drug_agent))
  
  return(in_table_ab)
}
