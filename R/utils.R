################################################################################
#    integrate the expression data of ortholog matrix
# ------------------------------------------------------------------------------

#' assign expression data to the ortholog matrix
#'
#' Compare between the ortholog matrix and expression data, hence obtain the
#' expression matrix for the ortholog genes across all compared species.
#' Theoretically, this function can also be used to integrate other types of
#' data apart from the expression files.
#'
#' @param x
#' Data with information for paired 1-to-1 orthologous pairs. It should be the
#' list generated using \code{\link{ortholog_match}}.
#'
#' @param meta_table
#' The metatable with the details for the expression data of all involved
#' species. This table should be of csv format.
#'
#' @param file_dir
#' The directory that containing the expression data files. Default: "."
#'
#' @return
#' Matrix, with rows as GeneID or GeneName, and columns as species. This matrix
#' can be used for subsequent analysis and visualization.
#'
#' @examples
#' # get ortholog matching information
#' species5.orth <- ortholog_match(
#'   species_list = c(
#'   "human", "mouse", "cattle", "pig", "opossum"
#'   )
#' )
#' 
#' # extract ortholog expression matrix for the compared species
#' species5.expr <- assign_expression_data(
#'   x = species5.orth,
#'   meta_table = "expression_metatable.csv",
#'   file_dir   = "inst/extdata"
#' )
#' 
#' # check top lines of the expression matrix
#' head(Species5.expr)
#' 
#' # write the expression matrix to a file
#' write.table(
#'   species5.expr,
#'   file = "Species5_expr.txt",
#'   sep = "\t", quote = FALSE,
#'    row.names = TRUE, col.names = TRUE
#' )
#'  
#' @export
assign_expression_data <- function (x, meta_table, file_dir = ".") {
  
  ## check parameters
  # check ortholog data
  if(missing(x)) {
    stop("x not specified.")
  }
  x <- verify_ortholog_data(x)
  
  # check expression data
  if(missing(meta_table)) {
    stop("meta_table not specified.")
  }
  
  # check file directory
  if(missing(file_dir)) {
    file_dir = "."
    warning("file_dir is missed. Assign to . by default.")
  }
  
  #x <- ortholog_match(
  #  species_list = c(
  #    "human", "mouse", "cattle", "pig", "opossum"
  #  )
  #)
  #saveRDS(x, "inst/extdata/species5_orth.rds")
  
  ## read meta table
  meta <- read.csv(meta_table)
  
  ## make the empty expression matrix
  mat_expr <- matrix(
    rep(NA, nrow(x[[1]]) * nrow(meta)),
    nrow = nrow(x[[1]])
  )
  # colnames: sample ID
  colnames(mat_expr) <- paste0(meta$Species, "_r", meta$RepIndex)
  # rownames: the gene_attr for the anchor species
  x.anchor <- x[[attr(x, "Species_anchor")]]
  row.names(mat_expr) <- paste0(
    x.anchor[,colnames(x.anchor) == "GeneID"], 
    ":",
    x.anchor[,colnames(x.anchor) == "GeneName"]
  )

  ## get geneID matrix, which will be used to extract expression data
  mat_id <- extract_ortholog_matrix(x, "GeneID")
  
  ## merge data to the expression matrix
  for(i in 1:nrow(meta)){
    # set parameters
    this.File       <- paste0(file_dir, "/", meta$File[i])
    this.Species    <- meta$Species[i]
    this.RepIdx     <- meta$RepIndex[i]
    this.IdColumn   <- meta$IdColumn[i]
    this.ExprColumn <- meta$ExprColumn[i]
    this.SampleName <- paste0(this.Species, "_r", this.RepIdx)
    if(!this.Species %in% colnames(mat_id)){
      stop(
        paste0("Error: ", this.Species, " is absent from the ortholog species.")
      )
    }
    this.GeneID     <- mat_id[,colnames(mat_id) == this.Species]
    print(
      paste0("Processing ", this.SampleName, ": ", this.File)
    )
    
    # read data
    this.data <- read.table(this.File, header = TRUE, sep = "\t")

    # assign data to expression matrix
    this.expr <- mapvalues(
      this.GeneID,
      from = this.data[,this.IdColumn],
      to   = this.data[,this.ExprColumn],
      warn_missing = FALSE
    )
    # unassigned data will be converted to NA
    this.expr[grep("ENS", this.expr)] <- NA
    this.expr <- as.numeric(this.expr)
    # assign data  
    mat_expr[,colnames(mat_expr)==this.SampleName] <- this.expr
  }
  
  ## filter the expression matrix to exclude lines with missing values
  mat_expr.missed <- apply(mat_expr, 1, function(i){sum(is.na(i))})
  print(
    paste0(
      "Warn: ",
      sum(mat_expr.missed>0),
      " genes are discarded due to missed data. ",
      sum(mat_expr.missed == 0),
      " genes with full data are retained."
    )
  )
  mat_expr.filtered <- mat_expr[mat_expr.missed==0,]
  
  ## return expression data
  return(mat_expr.filtered)
  
}

#' plot the correlation heatmap based on expression matrix
#'
#' Plot the correlation heatmap based on the expression matrix of the ortholog
#' genes across multiple species.
#'
#' @param x
#' Ortholog expression matrix obtained by using \code{\link{assign_expression_data}}.
#'
#' @param method
#' The method for calculating correlation efficiency. It can be "p" for Pearson,
#' or "s" for Spearman. Default: "p"  
#'
#' @param file_name
#' The file name for the generated figure. Default: "cor_heatmap.svg"
#'
#' @param width
#' The width for the generated figure. Default: 6.
#'
#' @param height
#' The height for the generated figure. Default: 5.
#'
#' @param pointsize
#' The point size for the generated figure. Default: 12
#'
#' @examples
#' # plot the correlation heatmap of SVG format
#' plot_correlation_heatmap(
#'   x = species5.expr,
#'   method     = "s",
#'   file_name  = "cor_heatmap.svg",
#'   width      = 6,
#'   height     = 5,
#'   pointsize  = 12
#' )
#' 
#' @export
plot_correlation_heatmap <- function (
    x, method = "p", file_name = "cor_heatmap.svg",
    width = 6, height = 5, pointsize = 12, ...
) {
  
  ## check parameters
  # check ortholog data
  if(missing(x)) {
    stop("x not specified.")
  }

  ## generate SVG image
  #x <- species5.expr
  #file_name <- "cor_heatmap.svg"
  #method <- "p"
  svg(
    file = file_name,
    width = width, height = height,
    pointsize = pointsize
    )
  pheatmap::pheatmap(
    cor(log10(x+1)),
    method = method,
    ...
    )
  dev.off()
}

#' plot the PCA biplot based on expression matrix
#'
#' Plot the PCA biplot based on the expression matrix of the ortholog
#' genes across multiple species.
#'
#' @param x
#' Ortholog expression matrix obtained by using \code{\link{assign_expression_data}}.
#'

#' @param file_name
#' The file name for the generated figure. Default: "PCA_biplot.svg"
#'
#' @param width
#' The width for the generated figure. Default: 5.
#'
#' @param height
#' The height for the generated figure. Default: 4.5.
#'
#' @param pointsize
#' The point size for the generated figure. Default: 12
#'
#' @examples
#' # plot the PCA biplot of SVG format
#' plot_PCA_biplot(
#'   x = species5.expr,
#'   file_name  = "PCA_biplot.svg",
#'   width      = 5,
#'   height     = 4.5,
#'   point_size = 12
#' )
#' 
#' @export
plot_PCA_biplot <- function (
    x, file_name = "PCA_biplot.svg",
    width = 5, height = 3, pointsize = 12, ...
) {
  
  ## check parameters
  # check ortholog data
  if(missing(x)) {
    stop("x not specified.")
  }
  
  ## convert expr matrix to DESeq2 object
  # prepare data
  require(DESeq2)
  cnt.data <- round(as.matrix(x), 0)
  samples  <- colnames(cnt.data)
  species  <- sub("_r.$", "", samples)
  col.data <- data.frame(
    Species = factor(species),
    Sample = samples
  )
  # make DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cnt.data,
    colData = col.data,
    design = ~ Species
    )
  # normalize data
  dds <- DESeq2::estimateSizeFactors(dds)
  se <- SummarizedExperiment::SummarizedExperiment(
    log10(counts(dds, normalized=TRUE) + 1),
    colData = colData(dds)
  )
  # perform PCA
  pca.data <- DESeq2::plotPCA(
    DESeq2::DESeqTransform(se), intgroup = c("Species"),
    returnData = TRUE, ntop=500
    )
  percent.var <- round(100*attr(pca.data, "percentVar"))
  
  ## generate SVG image
  #width <- 5
  #height <- 3
  #pointsize <- 12
  svg(
    file = file_name, 
    width = width, height = height, pointsize = pointsize
  )
  require(ggplot2)
  ggplot(pca.data, aes(PC1, PC2, color = Species)) +
    theme_bw() + 
    theme(legend.position = "right", legend.direction = "vertical") +
    geom_point(size = 3, alpha = I(0.8)) +
    xlab(paste0("PC1: ", percent.var[1], "% variance")) +
    ylab(paste0("PC2: ", percent.var[2], "% variance"))
  dev.off()
}
