% Generated  by roxygen2: do not edit by hand
% Please edit documentation in R/SEMeff.R
\name{SEMeff}
\alias{SEMeff}
\title{SEM Effects}
\usage{
SEMeff(local = setwd(), SEMplot = TRUE)
}
\arguments{
\item{local}{Set here the path of directory where the results of  THRGIBBS* and GIBBS* are. If you are using the command setwd() to change the directory, you don't need to use this parameter.}

\item{SEMplot}{It is a logical name, indicating if the SEM plot should be created.}
}
\value{
Mean and SD for each covariate in the SEM model.
}
\description{
R functions for summarize effects from strucutal equation models.
}
