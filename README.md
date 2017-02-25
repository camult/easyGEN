# easyGEN

R package with utility functions to help with data analysis

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
