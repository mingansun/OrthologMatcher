#' OrthologMatcher: a package to extract the 1-to-1 orthologous genes across multiple species
#'
#' The package can determine the 1-to-1 orthologous genes across multiple species
#' by automaically retriving the Ensembl ortholog annotation data. The generated
#' orthologous matching result can be output as a matrix or detailed lists, which
#' can be used for further evolutionary analysis.
#'
#' @section OrthologMatcher functions:
#'
#' \code{\link{ortholog_match}}
#' Match the 1-to-1 orthologous genes across multiple species
#'
#' \code{\link{extract_ortholog_matrix}}
#' Convert the 1-to-1 orthologous information across multiple species to a matrix
#'
#' \code{\link{list_species}}
#' List all supported species
#' 
#' \code{\link{summarize_ortholog_gene}}
#' Summarize the orthologues for each genetype or chromosome
#' 
#' \code{\link{verify_ortholog_data}}
#' Verify if the orthologue matching data is normal 
#'
#' @name OrthologMatcher-package
NULL