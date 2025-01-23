require(AMR)
require(tidyverse)

as.gene <- function(x){
  # Check x is a character vector
  if(!is.character(x)) {
    stop("Input must be a character vector.")
  }
  # Create the object with class "gene"
  structure(x, class = c("gene", class(x)))
}

# Define a print method for the "gene" class
print.gene <- function(x, ...) {
  cat("An object of class 'gene':\n")
  print(unclass(x), ...)
}

amrfp_drugs <- read_tsv("amrfp_drug_classes_agents.tsv")

parse_amrfp <- function(input_table, amrfp_drugs){
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
  in_table_ab <- in_table_subclass_split %>% left_join(., amrfp_drugs, by=c("Subclass"="AFP_Subclass"))
  
  # turn the drug agent column into the ab class with as.ab
  
  return(in_table_ab)
}