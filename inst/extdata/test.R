library(devtools)

setwd("X:/0-PRJ/0000-Software/2026-OrthologMatcher/OrthologMatcher")


document()
check()
load_all()

## extract matrix

species5.orth.GeneID.matrix <- extract_ortholog_matrix(
  species5.orth, "GeneID"
)
head(species5.orth.GeneID.matrix)

species5.orth.GeneName.matrix <- extract_ortholog_matrix(
  species5.orth, "GeneName"
)
head(species5.orth.GeneName.matrix)
  
  
## assign expression data
species5.orth <- ortholog_match(
 species_list = c(
 "human", "mouse",
 "cattle", "pig", "opossum"
 )
)
saveRDS(x, "inst/extdata/species5_orth.rds")

species5.expr <- assign_expression_data(
  x = species5.orth,
  meta_table = "inst/extdata/expression_metatable.csv",
  file_dir = "inst/extdata/"
  )

head(species5.expr)

## plot expression heatmap
plot_correlation_heatmap(
  x = species5.expr,
  method     = "s",
  file_name  = "cor_heatmap2.svg",
  width      = 6,
  height     = 5,
  pointsize  = 12
 )

## PCA
plot_PCA_biplot(
 x = species5.expr,
 file_name  = "PCA_biplot3.svg",
 width      = 5,
 height     = 4.5,
 point_size = 12
)




