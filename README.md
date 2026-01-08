# OrthologMatcher: an R package to extract the 1-to-1 orthologous genes across multiple speciess
---
__OrthologMatcher__ is an R package to can determine the 1-to-1 orthologous genes
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
Please refer to the R manual for the details how to use __OrthologMatcher__
to determine and output the interspecies ortholog information by using its
two major function: __ortholog_match__ and __extract_ortholog_matrix__.
