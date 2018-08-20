# lionessR
## Package for single sample network reconstruction in R

LIONESS, or **L**inear **I**nterpolation to **O**btain **N**etwork **E**stimates for **S**ingle **S**amples, can be used to reconstruct single-sample networks (http://arxiv.org/pdf/1505.06440.pdf). This package implements the LIONESS equation in R to reconstruct single-sample networks. The default network reconstruction method we use here is based on Pearson correlation. However, lionessR can run on any network reconstruction algorithms that returns a complete, weighted adjacency matrix. lionessR works for both unipartite and bipartite networks.

The easiest way to install the R package lionessR is via the devtools package from CRAN:
```
install.packages("devtools")
library(devtools)
devtools::install_github("mararie/lionessR")
```
And then load the package using: ```library(lionessR)```.

Please see the vignette for an example of co-expression network analysis using lionessR. For the example, you need to have CRAN packages `igraph` and `reshape2` and Bioconductor package `limma` installed.
