% Generated  by roxygen2: do not edit by hand
% Please edit documentation in R/PostGibbs.R
\name{PostGibbs}
\alias{PostGibbs}
\title{Summary of Multi-trait Bayesian Models23}
\usage{
PostGibbs(local = getwd(), burnIn = 0, thinning = 1, line = 1,
  HPD = 0.95, Names = c(), ICgraph = FALSE, Summary = TRUE, covAM = 1,
  DFS = FALSE)
}
\arguments{
\item{local}{Set here the path of directory where the results of 
THRGIBBS* and GIBBS* are. 
If you are using the command setwd() to change the directory, 
you don't need to use this parameter.}

\item{burnIn}{It is the number of iterations used as burn-in.}

\item{thinning}{It is the number of the thinning interval.}

\item{line}{It is the line width relative to the default (default=1). 2 is twice as wide.}

\item{HPD}{It is the Highest Posterior Density (HPD) intervals for the parameters in an MCMC sample.}

\item{Names}{It is a vector with the names of traits.}

\item{ICgraph}{It is a logical name, indicating if the IC Graph should be created.}

\item{Summary}{It is a logical name, indicating if the statistical summary should be run.}

\item{covAM}{type 0 if covariance between additive and maternal genetic effects must be fixed in zero, 1 otherwise. By default covAM=1.}

\item{DFS}{It is a logical name. Set it to TRUE to remove any adjacent cycles if any by using Depth First Search (DFS).}
}
\value{
A PDF file with graphics that summarizes an as.mcmc.obj 
        with a trace of the sampled output, a density estimate for each variable 
        in the chain, and the evolution of the sample quantiles as a function 
        of the number of iterations. 
        A XLS or CSV file with summary statistic for the marginal posteriors of 
        (co)variance components and parameters estimated. 
        So, it will be computed a sets of summary statistics for each variable: 
        mean; mode; median; standard deviation;  time-series standard error 
        based on an estimate of the spectral density at 0; Highest Posterior Density (HPD) 
        intervals (2.5 and 97.5) of the sample distribution. 
        It also possible return a graph 
        in PDF of causal relationship between tratis (Valente et. al., 2010).
}
\description{
This package is easy to use and can be helpful to summarize Bayesian Gibbs results 
when you run BLUPf90 programs (Misztal et al., 2014) and DMU packages (Madsen et al., 2006). 
Note that the random effects that can be read are: genetic (direct and maternal), 
permanent environmental (maternal or animal, not both), and residual. 
Genetic and environmental parameters are estimated as described by Willham (1972). 
Pay attention if you are fitting a model with maternal+direct effects 
because if the maternal effects came after the direct effects, this package will not run. 
So, if you want to use this package, the maternal effects must come before the direct effects. 
This package has also an option "ICgraph" that get partially, oriented graph, connected pairs 
and unshielded colliders, using Complete IC Algorithm adapted to work with variable 
number of traits (Valente et al., 2010)!
}
\examples{

## Not run:
## For an example of the use of a terms object as a formula

## If you are using setwd() function, you might do:
## setwd("/user/yourdir/")
## PostGibbs(burnIn=0, thinning=1) --> burnIn=0 and thinning=1 are defult

## Or, if are not using setwd() function, you should do:
## PostGibbs(local="/user/yourdir/")

## If you would like to use Inductive Causation (IC) algorithm to create
## causal graph, as discribed by Valente et. al., (2010), do it:
## PostGibbs(HPD=.95, Names=c("Var Name 1", "Var Name 2", ..., "Var Name n"), ICgraph=TRUE)
 
## End(Not run)

}
\references{
Misztal I, Tsuruta S.,  Lourenco D., Aguilar I., Legarra A., 
            Vitezica Z (2014). Manual  for BLUPF90   family	of programs. 
            Available at: <http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all1.pdf>.

Madsen P, Sorensen P, Su G, Damgaard LH, Thomsen H, Labouriau R. 
            DMU-a package for analyzing multivariate mixed models. In Proceedings 
            of the 8th World Congress on Genetics Applied to Livestock Production: 
            13-18 August 2006; Belo Horizonte. 2006.

Valente BD, Rosa GJM, de los Campos G, Gianola D, Silva MA (2010). 
            Searching for recursive causal structures in multivariate quantitative 
            genetics mixed models. Genetics.185:633-644.

Willham RL (1972). The role of maternal effects in animal breeding: 
            III. Biometrical aspects of maternal effects in animals. 
            Journal of Animal Science 35:1288-1295.
}
