# OrthologMatcher - an R package to match the 1-to-1 orthologous genes across multiple species
---
__OrthologMatcher__ is an R package to match the 1-to-1 orthologous genes across multiple species. It
can automaically retrive, process, and summarize the ortholog annotation data from Ensembl database.
The generated orthologous matching result can be output as a matrix or detailed lists, which can be
used for further evolutionary analysis.

## How to install OrthologMatcher
__OrthologMatcher__ can be installed using the __install_git__ function from devtools
package. However, it depends on several other R packages, which should be 
installed first.

__To install OrthologMatcher and its dependencies:__

install.packages("BiocManager")

BiocManager::install("biomaRt")

install.packages("dplyr")

devtools::install_git("https://github.com/mingansun/OrthologMatcher")

__To load the OrthologMatcher package:__

library(OrthologMatcher)

## How to use OrthologMatcher
Please refer to the Vignette about how to use __OrthologMatcher__ to extract and summarize
the ortholog matching data by using the two main function: __ortholog_match__ and 
__extract__ortholog_matrix__.


