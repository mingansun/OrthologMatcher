################################################################################
#    match 1-to-1 orthologous genes across multiple species
# ------------------------------------------------------------------------------
# list_species
# summarize_ortholog_gene
# verify_ortholog_data
# extract_ortholog_matrix
################################################################################


############################################################################
#' List all supported species
#'
#' Only species with genome available in ENSEMBL database are supported.
#' This function will show the information for all these species, or those filtered
#' by the user specified pattern.
#'
#' @param pattern
#' The pattern to filtering the species list. The pattern will be searched
#' against the common and scientific names, with the information for matched
#' species returned. To be noted, the matching is case insensitive.
#'
#' @return
#' Data frame, with the information for supported species
#'
#' @examples
#' # list all species based on locally saved species list file
#' list_species()
#'
#' # list all species with common or scientific name matching given pattern
#' list_species(pattern = "human")
#' list_species(pattern = "sapiens")
#' list_species(pattern = "macaque")
#'
#' @export
list_species <- function(pattern){
  
  opts <- options(stringsAsFactors = FALSE)
  
  if(missing(pattern)){
    pattern = ""
  }
  
  df_species <- NULL
  # fetch species list
  ensembl <- biomaRt::useMart(biomart = "ensembl")
  df_species <- biomaRt::listDatasets(mart = ensembl)
  # reformat
  colnames(df_species) <- c("Dataset", "Species", "Version")
  df_species$Dataset <- sub("_gene_ensembl$","", df_species$Dataset)
  df_species$Species <- sub("\\s+genes\\s+\\(.*$","", df_species$Species)
  
  # reformat and filter species list
  idx.flt <- sort(
    unique(
      c(
    grep(pattern = pattern, df_species$Dataset, ignore.case = TRUE),
    grep(pattern = pattern, df_species$Species, ignore.case = TRUE)
      )
    )
  )
  df_species.flt <- df_species[idx.flt,]
  
  # return
  return(df_species.flt)
  on.exit(options(opts, add = TRUE))
}

################################################################################
#' Summarize the orthologues for each genetype or chromosome
#'
#' @param x
#' Data frame with ortholog annotations for corresponding species. It is generated
#' using ortholog_match(species_list=c()).
#' @param group
#' genetype or chrom, the number of genes for each group based on specified
#' tags will be listed.
#'
#' @return
#' Data frame with the count of genes for each genetype or chrom in each species.
#'
#' @examples
#' hs2mm2ss.orth <- ortholog_match(species_list=c("human", "mouse","pig"))
#' summarize_ortholog_gene(hs2mm2ss.orth, "genetype")
#'
#' @export
summarize_ortholog_gene <- function(x, group = c("genetype", "chrom")){
  
  # check x
  if(missing(x)){
    stop("x not specified.")
  }
  x <- verify_ortholog_data(x)
  
  # check group
  if(missing(group)){
    stop("group not specified. Should be: genetype or chrom.")
  }
  if(length(group)>1){
    stop("group should be: genetype or chrom.")
  }
  group <- tolower(group)
  if(!group %in% c("genetype", "chrom")){
    stop("group should be: genetype or crom.")
  }
  
  # retrieve and summary information
  species   <- attr(x, "Species")
  gene_grp  <- NULL
  gene_cnt  <- list()
  
  for(i in 1:length(x)) {
    if(group == "genetype") {
      gene_grp      <- unique(c(gene_grp, x[[i]]$GeneType))
      gene_cnt[[i]] <- table(x[[i]]$GeneType)
    }
    if(group == "chrom") {
      gene_grp      <- unique(c(gene_grp, x[[i]]$Chrom))
      gene_cnt[[i]] <- table(x[[i]]$Chrom)
    }
  }
  
  gene_grp <- sort(gene_grp)
  
  mat <- matrix(0, nrow = length(gene_grp), ncol = length(species))
  row.names(mat) <- gene_grp
  colnames(mat)  <- species
  
  for(i in seq_along(species)) {
    mat[match(names(gene_cnt[[i]]), gene_grp),i] = as.vector(gene_cnt[[i]])
  }
  
  return(mat)
}

#' Verify if the orthologue matching data is normal
#'
#' Examine given data to see if its format is normal. Further, if no "Species"
#' attribute is present, try to set as Species_1, Species_2 etc, and return
#' the updated data. The species name is required for some other functions.
#'
#' @param x
#' Data with information for paired 1-to-1 orthologous pairs. It should be the
#' list generated using \code{\link{ortholog_match}}.
#'
#' @return
#' List, unchanged/updated data with the same structre as the input data x.
#'
#' @examples
#' hs2mm.orth <- ortholog_match(c("human", "mouse"))
#' hs2mm.orth <- verify_ortholog_data(hs2mm.orth)
#'
#' @export
verify_ortholog_data <- function (x) {
  
  if(missing(x)) {
    stop("x not specified.")
  }
  if(!is.list(x)) {
    stop("x should be an object generated using ortholog_match.")
  }
  if(length(x) == 0) {
    stop("x is empty.")
  }
  if(length(x) < 2) {
    warning("x has ", length(x), " items. Expect to be >=2.")
  }
  for(i in length(x)) {
    if(nrow(x[[i]]) < 1){
      stop("Element ", i, " of x has 0 row.")
    }
    if(ncol(x[[i]]) != 4) {
      stop("Element ", i, " of x has ", ncol(x[[i]]), " column. Expect 4.")
    }
  }
  # set Species as Species_n if no "Species" attribute
  if(length(attr(x, "Species")) == 0) {
    attr(x,"Species") <- paste0("Species_", 1:length(x))
    warning("x has no Species attribute. Set to Species_1, Species_2 etc.")
  }
  if(length(attr(x, "Species")) != length(x)) {
    attr(x,"Species") <- paste0("Species_", 1:length(x))
    warning("attr(x,\"Species\") has different length to x. Change to Species_1, Species_2 etc.")
  }
  
  return(x)
}


#' Convert the 1-to-1 orthologous information across multiple species to a matrix
#'
#' Examine given data to see if its format is normal. Further, if no "Species"
#' attribute is present, try to set as Species_1, Species_2 etc, and return
#' the updated data. The species name is required for some other functions.
#'
#' @param x
#' Data with information for paired 1-to-1 orthologous pairs. It should be the
#' list generated using \code{\link{ortholog_match}}.
#'
#' @param gene_attr
#' The information to be merged into the matrix. It should be GeneID or GeneName,
#' which are available for each species in the list generated using
#' \code{\link{ortholog_match}}. Default: GeneID.
#'
#'
#' @return
#' Matrix, with rows as GeneID or GeneName, and columns as species. This matrix
#' can be used for subsequent analysis of the 1:1 orthologs across all the involved
#' species. This is different from the list generated using \code{\link{ortholog_match}},
#' which contains rich information of the orthologs for each species, yet are not
#' merged together as a matrix.
#'
#' @examples
#' hs2mm2ss.orth <- ortholog_match(c("human", "mouse", "pig"))
#' hs2mm2ss.orth.GeneID.matrix <- extract_ortholog_matrix(hs2mm2ss.orth, "GeneID")
#'
#' @export
extract_ortholog_matrix <- function (x, gene_attr = "GeneID") {
  
  if(missing(x)) {
    stop("x not specified.")
  }
  x <- verify_ortholog_data(x)
  
  if(!gene_attr %in% c("GeneID", "GeneName")){
    stop("gene_attr should be GeneID or GeneName.")
  }
  
  ## make the matrix
  mat <- matrix(
    rep("", nrow(x[[1]]) * length(x)),
    nrow = nrow(x[[1]])
  )
  # colnames: species
  colnames(mat) <- attr(x, "Species")
  # rownames: the gene_attr for the anchor species
  x.anchor <- x[[attr(x, "Species_anchor")]]
  row.names(mat) <- x.anchor[,colnames(x.anchor) == gene_attr]
  
  ## merge data to the matrix
  for(i in 1:length(x)){
    z <- x[[i]]
    mat[,i] <- z[,colnames(z) == gene_attr]
  }
  
  return(mat)
}
