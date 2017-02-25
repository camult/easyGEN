# easyGEN

R package with utility functions to help with data quantitative genetics analysis, create parameter file, run REMLf90 and GIBBS2f90, as well as to summarize the bayesian results from GIBBS2f90 and DMU package.

This package has also an option "ICgraph" that get partially, oriented graph, connected pairs and unshielded colliders, using Complete IC Algorithm adapted to work with variable number of traits (Valente et al., 2010).


# How to Install

To install this package, use devtools:

devtools::install_github("camult/easyGEN")

Some functions require some bioconductor packages in order to work. As this package is not part of bioconductor, devtools will not automatically install it for you (see this thread for more details). You will have to install them manually:

source("https://bioconductor.org/biocLite.R")

biocLite("Rgraphviz")


# Overview

To see the full list of exported functions:

library("tinyutils")

ls(package:easyGEN)
