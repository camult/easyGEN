# `remlf90`: R functions for interfacing remlf90 software

## Description


 The function uses the model formula language in R to describe the model and from this generates the files needed to do an analysis using remlf90 software!


## Usage

```r
remlf90(
  formula,
  phen,
  ped = NULL,
  geno = NULL,
  map = NULL,
  idName,
  diffVAR = NULL,
  nMaternal = NULL,
  weight = NULL,
  Inb = FALSE,
  covariate = 0,
  OPTeff = NULL,
  OPTlist = NULL,
  missing = 0,
  Gcov = NULL,
  Rcov = NULL,
  execute = TRUE,
  PED_DEPTH = 0,
  covAM = 1,
  covR = 1,
  useF = FALSE
)
```


## Arguments

Argument      |Description
------------- |----------------
```formula```     |     A two-sided linear formula object describing the fixed effects part of the model, with the responses on the left of a ~ operator (separated by "|" if more then one response) and the terms, separated by + operators, on the right.
```phen```     |     Name of phenotype file.
```ped```     |     Name of pedigree file
```geno```     |     Name of genotype file.
```map```     |     Name of map file.
```idName```     |     Identification name of animal's column.
```diffVAR```     |     List of the trait and their unique effect, if any.
```nMaternal```     |     Number of traits on maternal effect. By default nMaternal=0.
```weight```     |     name of the column (numeric) with a vector of weights, may be NULL.  If weights is not NULL, the residual variance of each data-point is set to be proportional to the inverse of the squared-weight.
```Inb```     |     Whether to run the inbupgf90 to compute the coefficient of inbreeding. By default Inb=FALSE.
```covariate```     |     Name of covariate(s) if there is(are) any.
```OPTeff```     |     OPTIONAL effect, i.e., mat, pe, mpe. It is NULL if it is not given.
```OPTlist```     |     OPTIONAL effect. It is NULL if it is not given.
```missing```     |     a integer specifying the missing value (default 0).
```Gcov```     |     Genetic (co)variance. It is NULL if it is not given.
```Rcov```     |     Residual (co)variance. It is NULL if it is not given.
```execute```     |     Whether to run the remlf90. By default execute=TRUE.
```PED_DEPTH```     |     a integer specifying the depth of pedigree search. The default is 3, byt all pedigrees are loaded if it is set to 0.
```covAM```     |     type 0 if covariance between additive and maternal genetic effects must be fixed in zero, 1 otherwise. By default covAM=1.
```covR```     |     type 0 if covariance between additive and maternal genetic effects must be fixed in zero, 1 otherwise. By default covAM=1.
```useF```     |     logical value indicating whether inbreeding coefficient of this animal should be computed. Default is FALSE.

## Value


 BLUPf90 results.


## References


 Misztal I, Tsuruta S.,  Lourenco D., Aguilar I., Legarra A.,
 Vitezica Z (2014). Manual  for BLUPF90   familyof programs.
 Available at: <http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all1.pdf>.


# `gibbsf90`: A function interfacing gibbs2f90

## Description


 The function uses the model formula language in R to describe the model and from this generates the files needed to do an analysis using gibbs2f90 software!


## Usage

```r
gibbsf90(
  formula,
  phen,
  ped = NULL,
  geno = NULL,
  map = NULL,
  idName,
  diffVAR = NULL,
  nMaternal = NULL,
  weight = NULL,
  nIter = 1500,
  burnIn = 500,
  thin = 5,
  missing = 0,
  covariate = 0,
  OPTeff = NULL,
  OPTlist = NULL,
  intern = TRUE,
  Gcov = NULL,
  Rcov = NULL,
  execute = TRUE,
  PED_DEPTH = 0,
  covAM = 1,
  covR = 1,
  useF = FALSE
)
```


## Arguments

Argument      |Description
------------- |----------------
```formula```     |     A two-sided linear formula object describing the fixed effects part of the model, with the responses on the left of a ~ operator (separated by | if more then one response) and the terms, separated by + operators, on the right.
```phen```     |     Name of phenotype file.
```ped```     |     Name of pedigree file
```geno```     |     Name of genotype file.
```map```     |     Name of map file.
```idName```     |     Identification name of animal's column.
```diffVAR```     |     List of the trait and their unique effect, if any.
```nMaternal```     |     Number of traits on maternal effect. By default nMaternal=0.
```weight```     |     name of the column (numeric) with a vector of weights, may be NULL.  If weights is not NULL, the residual variance of each data-point is set to be proportional to the inverse of the squared-weight.
```nIter, burnIn, thin```     |     (integer) the number of iterations, burn-in and thinning.
```missing```     |     a integer specifying the missing value (default 0).
```covariate```     |     Name of covariate(s) if there is(are) any.
```OPTeff```     |     OPTIONAL effect, i.e., mat, pe, mpe. It is NULL if it is not given.
```OPTlist```     |     OPTIONAL effect. It is NULL if it is not given.
```intern```     |     a logical (not NA) which indicates whether to capture the output of the command as an R character vector.
```Gcov```     |     Genetic (co)variance. It is NULL if it is not given.
```Rcov```     |     Residual (co)variance. It is NULL if it is not given.
```execute```     |     Whether to run the gibbsf90. By default execute=TRUE.
```PED_DEPTH```     |     a integer specifying the depth of pedigree search. The default is 3, byt all pedigrees are loaded if it is set to 0.
```covAM```     |     type 0 if covariance between additive and maternal genetic effects must be fixed in zero, 1 otherwise. By default covAM=1.
```covR```     |     type 0 if residual covariance must be fixed in zero, 1 otherwise. By default covAM=1.
```useF```     |     logical value indicating whether inbreeding coefficient of this animal should be computed. Default is FALSE.

## Value


 gibbsf90 results.


## References


 Misztal I, Tsuruta S.,  Lourenco D., Aguilar I., Legarra A.,
 Vitezica Z (2014). Manual  for BLUPF90   familyof programs.
 Available at: <http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all1.pdf>.


# `PostGibbs`: Summary of Multi-trait Bayesian Models

## Description


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


## Usage

```r
PostGibbs(
  local = getwd(),
  burnIn = 0,
  thinning = 1,
  line = 1,
  HPD = 0.95,
  Names = c(),
  ICgraph = FALSE,
  Summary = TRUE,
  covAM = 1,
  DFS = FALSE
)
```


## Arguments

Argument      |Description
------------- |----------------
```local```     |     Set here the path of directory where the results of THRGIBBS* and GIBBS* are. If you are using the command setwd() to change the directory, you don't need to use this parameter.
```burnIn```     |     It is the number of iterations used as burn-in.
```thinning```     |     It is the number of the thinning interval.
```line```     |     It is the line width relative to the default (default=1). 2 is twice as wide.
```HPD```     |     It is the Highest Posterior Density (HPD) intervals for the parameters in an MCMC sample.
```Names```     |     It is a vector with the names of traits.
```ICgraph```     |     It is a logical name, indicating if the IC Graph should be created.
```Summary```     |     It is a logical name, indicating if the statistical summary should be run.
```covAM```     |     type 0 if covariance between additive and maternal genetic effects must be fixed in zero, 1 otherwise. By default covAM=1.
```DFS```     |     It is a logical name. Set it to TRUE to remove any adjacent cycles if any by using Depth First Search (DFS).

## Value


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


## References


 Misztal I, Tsuruta S.,  Lourenco D., Aguilar I., Legarra A.,
 Vitezica Z (2014). Manual  for BLUPF90   familyof programs.
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


## Examples

```r 
 
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
 
 ``` 

# `SEMeff`: SEM Effects

## Description


 R functions for summarize effects from strucutal equation models.


## Usage

```r
SEMeff(local = setwd(), SEMplot = TRUE)
```


## Arguments

Argument      |Description
------------- |----------------
```local```     |     Set here the path of directory where the results of  THRGIBBS* and GIBBS* are. If you are using the command setwd() to change the directory, you don't need to use this parameter.
```SEMplot```     |     It is a logical name, indicating if the SEM plot should be created.

## Value


 Mean and SD for each covariate in the SEM model.


