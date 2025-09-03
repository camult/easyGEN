% Generated  by roxygen2: do not edit by hand23
% Please edit documentation in R/gibbsf90.R
\name{gibbsf90}
\alias{gibbsf90}
\title{A function interfacing gibbs2f90}
\usage{
gibbsf90(formula, phen, ped = NULL, geno = NULL, map = NULL, idName,
  diffVAR = NULL, nMaternal = NULL, weight = NULL, nIter = 1500,
  burnIn = 500, thin = 5, missing = 0, covariate = 0, OPTeff = NULL,
  OPTlist = NULL, intern = TRUE, Gcov = NULL, Rcov = NULL,
  execute = TRUE, PED_DEPTH = 0, covAM = 1, covR = 1, useF = FALSE)
}
\arguments{
\item{formula}{A two-sided linear formula object describing the fixed effects part of the model, 
with the responses on the left of a ~ operator (separated by | if more then one response)
and the terms, separated by + operators, on the right.}

\item{phen}{Name of phenotype file.}

\item{ped}{Name of pedigree file}

\item{geno}{Name of genotype file.}

\item{map}{Name of map file.}

\item{idName}{Identification name of animal's column.}

\item{diffVAR}{List of the trait and their unique effect, if any.}

\item{nMaternal}{Number of traits on maternal effect. By default nMaternal=0.}

\item{weight}{name of the column (numeric) with a vector of weights, may be NULL.  If weights is not NULL, the residual variance of each data-point is set to be proportional to the inverse of the squared-weight.}

\item{nIter, burnIn, thin}{(integer) the number of iterations, burn-in and thinning.}

\item{missing}{a integer specifying the missing value (default 0).}

\item{covariate}{Name of covariate(s) if there is(are) any.}

\item{OPTeff}{OPTIONAL effect, i.e., mat, pe, mpe. It is NULL if it is not given.}

\item{OPTlist}{OPTIONAL effect. It is NULL if it is not given.}

\item{intern}{a logical (not NA) which indicates whether to capture the output of the command as an R character vector.}

\item{Gcov}{Genetic (co)variance. It is NULL if it is not given.}

\item{Rcov}{Residual (co)variance. It is NULL if it is not given.}

\item{execute}{Whether to run the gibbsf90. By default execute=TRUE.}

\item{PED_DEPTH}{a integer specifying the depth of pedigree search. The default is 3, byt all pedigrees are loaded if it is set to 0.}

\item{covAM}{type 0 if covariance between additive and maternal genetic effects must be fixed in zero, 1 otherwise. By default covAM=1.}

\item{covR}{type 0 if residual covariance must be fixed in zero, 1 otherwise. By default covAM=1.}

\item{useF}{logical value indicating whether inbreeding coefficient of this animal should be computed. Default is FALSE.}
}
\value{
gibbsf90 results.
}
\description{
The function uses the model formula language in R to describe the model and from this generates the files needed to do an analysis using gibbs2f90 software!
}
\references{
Misztal I, Tsuruta S.,  Lourenco D., Aguilar I., Legarra A., 
            Vitezica Z (2022). Manual  for BLUPF90   family	of programs. 
            Available at: <https://nce.ads.uga.edu/html/projects/programs/docs/blupf90_all8.pdf>.
}
