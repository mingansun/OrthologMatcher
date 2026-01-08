################################################################################
#    match 1-to-1 orthologous genes across multiple species
# ------------------------------------------------------------------------------
# ortholog_match
################################################################################

#' Match the 1-to-1 orthologous genes across multiple species
#'
#' Ortholog matching is performed based on the ortholog information retrieved
#' from ENSEMBL. The R package biomaRt is used to access ENSEMBL.To work properly,
#' make sure biomaRt package is installed and the internet connection works.
#' It may take sometime to download data from ENSEMBL - depending on the
#' network speed. Only species with their genomes availale in ENSEMBL database
#' are supported. Use \code{\link{list_species}} to get the information for
#' all the supported species.
#'
#' @param species_list
#' A vector for all the species to be analyzed. Should be common names such as 
#' "human", "mouse", etc. Use \code{\link{list_species}} to get the names of all
#' supported species.
#' @param species_anchor
#' A species name from the provided species_list. This species will server as the
#' anchor to be compared against by all other species. Default:The first species
#' from the provided species_list.
#' @param host
#' ENSEMBL host (mirror) to connect to. Default: www.ensembl.org. By default, 
#' bioMart uses the latest releases of genome assemblies from the ENSEMBL database,
#' which may not match the genome version of given data sometimes. Earlier releases
#' can be used by specifying an archieved host address, which can be  found from
#' the ENSEMBL website.
#'
#' @return
#' A matrix with the information for 1-to-1 ortholog pairs. The list has corresponding 
#' species data frames (each for one species) with columns GeneID, GeneName, 
#' Chrom,GeneType. Additonal attributes including Species and SpeciesAbbr are 
#' also attached to the return list.
#' 
#' @examples
#' # Get 1:1 orthologs between different pairs of species
#' hs2mm.orth <- ortholog_match(species_list = c("human", "mouse"))
#' hs2mm2ss.orth <- ortholog_match(species_list = c("human", "mouse", "pig"))
#' 
#' #Specify mouse as the anchor species
#' hs2mm2ss.orth <- ortholog_match(
#'   species_list = c("human", "mouse", "pig"), species_anchor = "mouse"
#' )
#'
#' # Choose to use different host (mirrors) for ENSEMBL
#' hs2mm2ss.orth <- ortholog_match(
#'   species_list = c("human", "mouse", "pig"), host = "www.ensembl.org"
#' )
#' hs2mm2ss.orth <- ortholog_match(
#'   species_list = c("human", "mouse", "pig"), host = "useast.ensembl.org"
#' )
#' hs2mm2ss.orth <- ortholog_match(
#'   species_list = c("human", "mouse", "pig"), host = "asia.ensembl.org"
#' )
#'
#' # Choose to use an archived host for ENSEMBL, which contains old genome releases
#' hs2mm2ss.orth <- ortholog_match(
#'   species_list = c("human", "mouse", "pig"), host = "nov2020.archive.ensembl.org"
#' )
#' @export

ortholog_match <- function(species_list, species_anchor, host) {
  
  ## check species_list
  if(missing(species_list)) {
    stop("species_list not specified.")
  }
  if(length(species_list) < 2) {
    stop("At least two species are required for comparison.")
  }
 
  ## connect to ENSEMBL
  if(missing(host)){
    mart = biomaRt::useMart("ensembl", verbose = TRUE)
  }
  else{
    mart = biomaRt::useMart("ensembl", verbose = TRUE, host = host)
  }
  
  # check if all specified species are supported (by comparing against list_species()),
  # and get their abbreviations which will be used by biomaRt
  species.ok   <- list_species()
  species_list <- tolower(species_list)
  species_abbr <- vector("character", length(species_list))
  
  message("Get species name abbreviation to be used by biomaRt ...")
  for(sp in seq_along(species_list)){
    if(! species_list[sp] %in% tolower(species.ok$Species)){
      stop(
        species_list[sp], " is not supported. Use list_species() to list all supported species."
      )
    }
    species_abbr[sp] <- species.ok$Dataset[which(tolower(species.ok$Species) == species_list[sp])]
    message(species_list[sp], " => ", species_abbr[sp])
  }
  
  # biomaRt connection
  conn         <- list()
  species_data <- list() 
  
  message("\nConnect to biomaRt for each species ...")
  for(i in seq_along(species_list)){
    conn[[i]] <- biomaRt::useDataset(
      paste0(species_abbr[i], "_gene_ensembl"), mart, verbose = FALSE
    )
    species_data[[i]] <- biomaRt::getBM(
      attributes = c(
        "ensembl_gene_id",
        "external_gene_name",
        "chromosome_name",
        "gene_biotype"
      ),
      mart    = conn[[i]],
      filters = "ensembl_gene_id",
      values  = ""
    )
  }

  ##check anchor species
  anchor_idx <- NULL
  if(missing(species_anchor)){
    species_anchor <- species_list[1]
    message(
      "\nNo anchor species specified, using the first species as default: ",
      species_anchor
    )
    anchor_idx <- 1
  }
  else{
    if(species_anchor %in% species_list){
      anchor_idx <- which(species_list == species_anchor)
      message("\nUsing user-specified anchor species: ", species_anchor)
    }
    else{
      stop("\nAnchor species is not in the species list.")
    }
  }
     
  anchor   <- species_list[anchor_idx]
  anchor_abbr <- species_abbr[anchor_idx]
  anchor_data <- species_data[[anchor_idx]]
      
  message("anchor species: ", anchor, " (", nrow(anchor_data), " genes)")

  ## Get orthologous genes from the anchor species to all other species
  ortholog_list <- list()
  ortholog_list[[anchor]] <- anchor_data[order(anchor_data$ensembl_gene_id),]
      
  for(i in seq_along(species_list)){
    if(i == anchor_idx){
      next
    }

    target      <- species_list[i]
    target_abbr <- species_abbr[i]
    message(
      "\nSearching for orthologous genes between ", anchor,
      " and ", target, "..."
    )
        
    hom_tag  <- paste0(target_abbr, "_homolog_ensembl_gene")
    message(anchor, " => ", target, ": ", hom_tag)
    
    # pairwise check if the homolog table available for species
    attr_list <- biomaRt::listAttributes(conn[[anchor_idx]])$name
    if(!hom_tag %in% attr_list) {
      stop(paste0(anchor, " has no attribute: ", target))
    }
        
    homolog_info <- biomaRt::getBM(
      attributes = c("ensembl_gene_id","external_gene_name", hom_tag),
      mart    = conn[[anchor_idx]],
      filters = "ensembl_gene_id",
      values  = anchor_data$ensembl_gene_id
    )
   
    colnames(homolog_info) <- c("anchor_gene_id", "anchor_gene_name", "target_gene_id")
    homolog_info <- homolog_info[homolog_info$target_gene_id != "", ]
    message("orthologous genes: ",nrow(homolog_info))
   
    # 1:1 matched genes
    anchor_counts <- table(homolog_info$anchor_gene_id)
    target_counts <- table(homolog_info$target_gene_id)
   
    anchor_unique <- names(anchor_counts)[anchor_counts == 1]
    target_unique <- names(target_counts)[target_counts == 1]
        
    one_to_one <- homolog_info[
      homolog_info$anchor_gene_id %in% anchor_unique & 
      homolog_info$target_gene_id %in% target_unique,
    ]
        
    message("Found ", nrow(one_to_one), " 1:1 orthologous gene pairs.")
        
    # Get complete gene information for the target species
    target_gene_info <- species_data[[i]]
    target_orthologs <- merge(
      one_to_one,
      target_gene_info,
      by.x = "target_gene_id",
      by.y = "ensembl_gene_id",
      all.x = FALSE
    ) 

    colnames(target_orthologs)[colnames(target_orthologs) == "target_gene_id"] <- 
       paste0(target, "_gene_id")
    target_orthologs <- dplyr::select(target_orthologs, 2:3, 1, 4:6)
    ortholog_list[[target]] <- target_orthologs[order(target_orthologs$anchor_gene_id),]
  }
   
  ## merge data from all species
  message(
    "\nFilting orthologous genes across all compared species: ",
    paste(names(ortholog_list), collapse = ", ")
  )
  unified_matrix <- ortholog_list[[anchor]]
  for(sp in names(ortholog_list)){
    if(sp != anchor){
      sp_data <- ortholog_list[[sp]][, c("anchor_gene_id", paste0(sp, "_gene_id"))]
      # merge the common genes across species
      unified_matrix <- merge(
        unified_matrix,
        sp_data, 
        by.x  = "ensembl_gene_id", 
        by.y  = "anchor_gene_id", 
        all.x = FALSE
      )
    }
  }
  
  ## Rename column names and output list in the specified format  
  common_gene_ids <- unified_matrix$ensembl_gene_id
  formatted_orthologs <- list()
  for(sp in names(ortholog_list)){
    sp_data <- ortholog_list[[sp]]
    if(sp == anchor) {
      filtered_data <- sp_data[sp_data$ensembl_gene_id %in% common_gene_ids, ]
      formatted_df <- data.frame(
        GeneID   = filtered_data$ensembl_gene_id,
        GeneName = filtered_data$external_gene_name,
        Chrom    = filtered_data$chromosome_name,
        GeneType = filtered_data$gene_biotype,
        stringsAsFactors = FALSE
      )
    }
    else{
      filtered_data <- sp_data[sp_data$anchor_gene_id %in% common_gene_ids, ]
      formatted_df <- data.frame(
        GeneID   = filtered_data[,c(paste0(sp,"_gene_id"))],
        GeneName = filtered_data$external_gene_name,
        Chrom    = filtered_data$chromosome_name,
        GeneType = filtered_data$gene_biotype,
        stringsAsFactors = FALSE
      )
    }
    formatted_orthologs[[sp]] <- formatted_df
  }

  message("Found ", nrow(formatted_orthologs[[1]]), " 1:1 orthologous genes")  
  attr(formatted_orthologs, "Unified_matrix") <- unified_matrix
  attr(formatted_orthologs, "Species_anchor") <- anchor
  attr(formatted_orthologs, "Species")        <- species_list
  attr(formatted_orthologs, "SpeciesAbbr")    <- species_abbr
  attr(formatted_orthologs, "SpeciesAbbr")    <- species_abbr
  
  attr(formatted_orthologs, "GeneID_matrix")  <- extract_ortholog_matrix(formatted_orthologs, "GeneID")
  attr(formatted_orthologs, "GeneName_matrix")  <- extract_ortholog_matrix(formatted_orthologs, "GeneName")
  return(formatted_orthologs)  
}
