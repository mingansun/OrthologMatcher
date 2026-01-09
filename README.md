# OrthologMatcher: an R package to extract the 1-to-1 orthologous genes across multiple speciess
---
__OrthologMatcher__ is an R package to determine the 1-to-1 orthologous genes
across multiple species by automaically retriving the Ensembl ortholog annotation
data. The generated orthologous matching result can be output as a matrix or detailed
lists, which can be used for further evolutionary analysis.

## How to install OrthologMatcher
__OrthologMatcher__ can be installed using the __install_git__ function from devtools
package. However, it depends on several other R packages, which should be installed first.

__To install OrthologMatcher and its dependencies:__

install.packages("BiocManager")

BiocManager::install("biomaRt")

install.packages("dplyr")

install.packages("devtools")

devtools::install_git("https://github.com/mingansun/OrthologMatcher")

__To load the OrthologMatcher package:__

library(OrthologMatcher)

## How to use OrthologMatcher
You can use __OrthologMatcher__ to obtain the ortholog pairing information across multiple
species with its two major functions:

1) __ortholog_match__: this function determines the 1-to-1 orthologues across all the compared species. It
   return the results as an R object (a list), which not only contains the detailed information (e.g.,
   gene name, gene ID, gene type, genomic position) for each species, but also holds the ortholog
   pairing matrix (for GeneID and GeneName, respectively) for all the compared species as attributes.
   
3) __extract_ortholog_matrix__: this function extracts the ortholog pairing matrix from the R object produced
   by __ortholog__match__. Regarding the returned matrix, the rows represent the genes (GeneID or GeneName,
   depending on the user supplied parameter), and the columns represent the compared species.
   
Please refer to the R manual for more details.
