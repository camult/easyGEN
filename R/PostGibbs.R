#' @title Summary of Multi-trait Bayesian Models
#' 
#' @description This package is easy to use and can be helpful to summarize Bayesian Gibbs results 
#' when you run BLUPf90 programs (Misztal et al., 2014) and DMU packages (Madsen et al., 2006). 
#' Note that the random effects that can be read are: genetic (direct and maternal), 
#' permanent environmental (maternal or animal, not both), and residual. 
#' Genetic and environmental parameters are estimated as described by Willham (1972). 
#' Pay attention if you are fitting a model with maternal+direct effects 
#' because if the maternal effects came after the direct effects, this package will not run. 
#' So, if you want to use this package, the maternal effects must come before the direct effects. 
#' This package has also an option "ICgraph" that get partially, oriented graph, connected pairs 
#' and unshielded colliders, using Complete IC Algorithm adapted to work with variable 
#' number of traits (Valente et al., 2010).
#' 
#' @param local Set here the path of directory where the results of 
#'        THRGIBBS* and GIBBS* are. 
#'        If you are using the command setwd() to change the directory, 
#'        you don't need to use this parameter.
#' @param burnIn It is the number of iterations used as burn-in.
#' @param thinning It is the number of the thinning interval.
#' @param line It is the line width relative to the default (default=1). 2 is twice as wide.
#' @param HPD It is the Highest Posterior Density (HPD) intervals for the parameters in an MCMC sample.
#' @param Names It is a vector with the names of traits.
#' @param ICgraph It is a logical name, indicating if the IC Graph should be created.
#' @param Summary It is a logical name, indicating if the statistical summary should be run.
#' @param covAM type 0 if covariance between additive and maternal genetic effects must be fixed in zero, 1 otherwise. By default covAM=1.
#' 
#' 
#' @return A PDF file with graphics that summarizes an as.mcmc.obj 
#'         with a trace of the sampled output, a density estimate for each variable 
#'         in the chain, and the evolution of the sample quantiles as a function 
#'         of the number of iterations. 
#'         A XLS or CSV file with summary statistic for the marginal posteriors of 
#'         (co)variance components and parameters estimated. 
#'         So, it will be computed a sets of summary statistics for each variable: 
#'         mean; mode; median; standard deviation;  time-series standard error 
#'         based on an estimate of the spectral density at 0; Highest Posterior Density (HPD) 
#'         intervals (2.5 and 97.5) of the sample distribution. 
#'         It also possible return a graph 
#'         in PDF of causal relationship between tratis (Valente et. al., 2010).
#'         
#' @references Misztal I, Tsuruta S.,  Lourenco D., Aguilar I., Legarra A., 
#'             Vitezica Z (2014). Manual  for BLUPF90   family	of programs. 
#'             Available at: <http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all1.pdf>.
#'  
#' @references Madsen P, Sorensen P, Su G, Damgaard LH, Thomsen H, Labouriau R. 
#'             DMU-a package for analyzing multivariate mixed models. In Proceedings 
#'             of the 8th World Congress on Genetics Applied to Livestock Production: 
#'             13-18 August 2006; Belo Horizonte. 2006. 
#'  
#' @references Valente BD, Rosa GJM, de los Campos G, Gianola D, Silva MA (2010). 
#'             Searching for recursive causal structures in multivariate quantitative 
#'             genetics mixed models. Genetics.185:633-644.
#'             
#' @references Willham RL (1972). The role of maternal effects in animal breeding: 
#'             III. Biometrical aspects of maternal effects in animals. 
#'             Journal of Animal Science 35:1288-1295.
#'             
#' @examples
#' 
#' ## Not run:
#' ## For an example of the use of a terms object as a formula
#' 
#' ## If you are using setwd() function, you might do:
#' ## setwd("/user/yourdir/")
#' ## PostGibbs(burnIn=0, thinning=1) --> burnIn=0 and thinning=1 are defult
#' 
#' ## Or, if are not using setwd() function, you should do:
#' ## PostGibbs(local="/user/yourdir/")
#' 
#' ## If you would like to use Inductive Causation (IC) algorithm to create
#' ## causal graph, as discribed by Valente et. al., (2010), do it:
#' ## PostGibbs(HPD=.95, Names=c("Var Name 1", "Var Name 2", ..., "Var Name n"), ICgraph=TRUE)
#'  
#' ## End(Not run)
#' 
#' @export PostGibbs
#' 
#' @import Rgraphviz WriteXLS coda Matrix graph gdata
##------------------------------------------------------------------------------------------##
PostGibbs=function(local=getwd(),burnIn=0,thinning=1,line=1,HPD=.95,
                   Names=c(),ICgraph=FALSE,Summary=TRUE, covAM=1){
  ##---------------------------------------------------------------------------------------##
  ## Matrix functions
  ##----------------------------------------------------------------------------------------##
  ginv <- function(X, tol = sqrt(.Machine$double.eps)){
    #
    # based on suggestions of R. M. Heiberger, T. M. Hesterberg and WNV
    #
    if(length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
      stop("'X' must be a numeric or complex matrix")
    if(!is.matrix(X)) X <- as.matrix(X)
    Xsvd <- svd(X)
    if(is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive)) Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if(!any(Positive)) array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop=FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop=FALSE]))
  }
  #------------------------------------------------------------------------------------------#
  # convert vector to symmetric matrix
  #------------------------------------------------------------------------------------------#
  vec2sm = function(vec, diag=TRUE, order=NULL) {
    # dimension of matrix
    n = (sqrt(1+8*length(vec))+1)/2
    if (diag == TRUE) n = n-1
    if ( ceiling(n) != floor(n) )
      stop("Length of vector incompatible with symmetric matrix")
    
    # fill lower triangle of matrix     
    m = matrix(NA, nrow=n, ncol=n)
    lo = lower.tri(m, diag)
    if (is.null(order))
    {
      m[lo] = vec
    }
    else
    {
      # sort vector according to order
      vec.in.order = rep(NA, length(order))
      vec.in.order[order] = vec
      m[lo] = vec.in.order
    }
    
    # symmetrize
    for (i in 1:(n-1))
      for (j in (i+1):n)
        m[i, j] = m[j, i]   
    
    return( m )
  }
  #------------------------------------------------------------------------------------------#
  # upperTriangle and lowerTriangle
  #------------------------------------------------------------------------------------------#
  upperTriangle <- function(x, diag=FALSE, byrow=FALSE){
    if(byrow)
      t(x)[rev(upper.tri(x, diag=diag))]
    else
      x[upper.tri(x, diag=diag)]
  }
  "upperTriangle<-" <- function(x, diag=FALSE, byrow=FALSE, value){
    if(byrow) {
      ret <- t(x)
      ret[rev(upper.tri(x, diag=diag))] <- value
      t(ret)
    }
    else {        
      x[upper.tri(x, diag=diag)] <- value
      x
    }
  }
  lowerTriangle <- function(x, diag=FALSE, byrow=FALSE){
    if(byrow)
      t(x)[rev(lower.tri(x, diag=diag))]
    else
      x[lower.tri(x, diag=diag)]
  }
  "lowerTriangle<-" <- function(x, diag=FALSE, byrow=FALSE, value){
    if(byrow) {
      ret <- t(x)
      ret[rev(lower.tri(x, diag=diag))] <- value
      t(ret)
    }
    else {        
      x[lower.tri(x, diag=diag)] <- value
      x
    }
  }
  ##----------------------------------------------------------------------------------------##
  ## Mode function
  ##----------------------------------------------------------------------------------------##
  estimate_mode <- function(x) {
    d <- density(x)
    d$x[which.max(d$y)]
  }
  ##----------------------------------------------------------------------------------------##
  burnIn=burnIn
  thinning=thinning
  line=line
  HPD=HPD
  nodeNames=Names
  ##----------------------------------------------------------------------------------------##
  cat("\014")
  cat("\n")
  centerText <- function() {
    width <- getOption("width")
    A <- ("                                                       ._____.    \n")
    B <- ("    _/////_            Fernando Brito Lopes           _|_____|_   \n")
    C <- ("   (' o o ')     Animal Scientist (Zootechnician)     (' o o ')   \n")
    D <- ("__ooO_(_)_Ooo_____ Animal Breeding and Genetics _____ooO_(_)_Ooo__\n")
    E <- ("                    e-mail: <camult@gmail.com>                    \n")
    ws <- rep(" ", floor((width - nchar(A))/2))
    cat(ws,A,ws,B,ws,C,ws,D,ws,E,ws,sep = "")
  }
  centerText()
  cat("\n")
  ##----------------------------------------------------------------------------------------##
  gib.samples <- list.files(pattern = "gib.samples")
  if(!identical(gib.samples,character(0))){
    DMUfile <- list.files(pattern = "gib.samples")
    postgibbs <- read.table(DMUfile, col.names=c("iter","effect","row","col","value","none"))
    #---------------------------------------------------------------------------------------#
    iter <- unique(postgibbs$iter)
    ntrait <- length(unique(postgibbs$row))
    nParam <- (((ntrait*(ntrait-1))/2)+ntrait)*(length(unique(postgibbs$effect)))
    #---------------------------------------------------------------------------------------#
    COV <- matrix(0,
                  nrow=length(iter),
                  ncol=3+nParam)
    for(i in 1:length(iter)){
      G <- postgibbs[postgibbs$iter==iter[i]&postgibbs$effect==1, names(postgibbs)%in%c("row","col","value")]
      R <- postgibbs[postgibbs$iter==iter[i]&postgibbs$effect==2, names(postgibbs)%in%c("row","col","value")]
      G <- G[order(G$col),]
      R <- R[order(R$col),]
      COV[i, ] <- c(i, iter[i], nParam, G$value, R$value)
    }
    write.table(COV, "postgibbs_samples", sep=" ", col.names=F, row.names=F, quote=F)
    write.table(as.matrix(c(0,0,0,0,ntrait)), "renf90.par", sep=" ", col.names=F, row.names=F, quote=F)
    #---------------------------------------------------------------------------------------#
    Gmatrix <- vec2sm(c(1:(nParam/2)), diag=TRUE)
    Gmatrix[lower.tri(Gmatrix, diag=FALSE)] <- 0
    Rmatrix <- vec2sm(c(((nParam/2)+1):(nParam)), diag=TRUE)
    Rmatrix[lower.tri(Rmatrix, diag=FALSE)] <- 0
    #---------------------------------------------------------------------------------------#
    write.fwf(as.matrix(" G matrix for effect ="), "postind", sep=" ", colnames=F, rownames=F, quote=F)
    write.fwf(Gmatrix, "postind", sep=" ", colnames=F, rownames=F, quote=F, append=TRUE, justify="right")
    write.fwf(as.matrix(" R matrix"), "postind", sep=" ", colnames=F, rownames=F, quote=F, append=TRUE)
    write.fwf(Rmatrix, "postind", sep=" ", colnames=F, rownames=F, quote=F, append=TRUE, justify="right")
    #---------------------------------------------------------------------------------------#
  }
  ##----------------------------------------------------------------------------------------##
  cat("Reading postgibbs_samples...\n")
  arq=as.matrix(read.table(paste(local,"/postgibbs_samples", sep="")))[,-c(1:3)]
  ##----------------------------------------------------------------------------------------##
  if(burnIn>nrow(arq)) { message("\n...\n") }
  if(burnIn>nrow(arq)) stop("Burn-in is greater than sample size...", call.=F)
  ##----------------------------------------------------------------------------------------##
  indexT<-c(1:((nrow(arq)-burnIn)/thinning))*thinning
  arq <- arq[burnIn:nrow(arq),]
  arq <- arq[indexT,]
  arq = as.matrix(arq)
  r=arq
  #------------------------------------------------------------------------------------------#
  # Determinando Modelo Utilizado
  #------------------------------------------------------------------------------------------#
  cat("\n")
  cat("Performing initial checks...\n")
  #------------------------------------------------------------------------------------------#
  traits=read.table("renf90.par", header=FALSE, skip=4, nrows =1)[1,1]
  mode(traits)="numeric"
  cat(traits, "traits", "were found")
  cat("\n")
  Matrix=as.matrix(strsplit(readLines("postind"), ", 0 "))
  #------------------------------------------------------------------------------------------#
  # Cheking if fixed effects are fitted as random effects
  #------------------------------------------------------------------------------------------#
  if(length(grep("matrix",Matrix))>3){ 
    Matrix=as.matrix(Matrix[-c(1:((traits+1)*traits)),])
  }
  #------------------------------------------------------------------------------------------#
  #
  #
  # The GIBBS Script start here
  #
  #
  #------------------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------------------#
  if(Summary==TRUE){
    #------------------------------------------------------------------------------------------#    
    #------------------------------------------------------------------------------------------#
    # These scripts will run if it were found two matrices
    #------------------------------------------------------------------------------------------#
    if(length(grep("matrix",Matrix))==2){
      #------------------------------------------------------------------------------------------#
      cat("\n")
      cat("Reading matrices from Bayesian analyses...\n")
      #------------------------------------------------------------------------------------------#    
      n1matrix=as.matrix(strsplit(as.character(Matrix[2]), " ")[[1]])
      mode(n1matrix)="numeric"
      n1matrix=c(length(na.omit(n1matrix)))
      mode(n1matrix)="numeric"
      if(n1matrix==1){
        cat("\n")  
        cat("You are summarizing the results from univariate analysis...")
        cat("\n")  
      }
      eff_G=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
      eff_G=as.matrix(diag(eff_G))
      eff_G=nnzero(eff_G)
      if(n1matrix!=1){
        cat("\n")  
        cat("You are summarizing the results from multi-trait analysis...")
        cat("\n")  
        cat("1st (co)variance matrix has size of", eff_G,"x",eff_G)
        cat("\n")  
      }
      #----------------------------------------------------------------------------------------#
      n2matrix=(length(Matrix)-2)-n1matrix
      if(n2matrix!=1){
        cat("2nd (co)variance matrix has size of", n2matrix,"x",n2matrix)
        cat("\n")
      }
      cat("\n")
      #----------------------------------------------------------------------------------------#
      last1M=as.matrix(strsplit(as.character(Matrix[2]), " ")[[1]])
      mode(last1M)="numeric"
      last1M=c(tail(na.omit(last1M), n=1))
      #----------------------------------------------------------------------------------------#
      n1M=as.matrix(strsplit(as.character(Matrix[traits+2]), " ")[[1]])
      mode(n1M)="numeric"
      n1M=rev(na.omit(n1M))[1] # Reading the last value
      #----------------------------------------------------------------------------------------#
      n2M=as.matrix(strsplit(as.character(Matrix[n1matrix+3]), " ")[[1]])
      mode(n2M)="numeric"
      n2M=na.omit(n2M)[2]
      #----------------------------------------------------------------------------------------#
      # Uni=A + R          --> A and R are vectors
      # A = A + M + R      --> M for all traits
      # B = A + M + R      --> M for only for pre weaning traits
      # C = A + R          --> A and R are matrix
      #----------------------------------------------------------------------------------------#
      Uni=(n1matrix==1) ##
      if(covAM==0){
        A=(n1matrix==traits*2 & last1M==0 & n1M!= 0 & n2M!=0) ## A and M for all traits
        B=(n1matrix==traits*2 & last1M==0 & n1M== 0 & n2M!=0) ## Maternal only for the firsts traits
      } else {
        A=(n1matrix==traits*2 & last1M!=0 & n2M!=0) ## A and M for all traits
        B=(n1matrix==traits*2 & last1M==0 & n2M!=0) ## Maternal only for the firsts traits
      }
      C=(n1matrix==traits & n2M!=0)               ## A and R are matrix
      if(is.na(C)){C=FALSE}
      #----------------------------------------------------------------------------------------#
      # D = A + M + R      --> M for all traits
      # E = A + M + R      --> M for only for pre weaning traits
      # F = A + R          --> A matrix and R is diagonal matrix
      #----------------------------------------------------------------------------------------#
      D=(n1matrix==traits*2 & last1M!=0 & n2M==0) ## M for all traits and diag(R)
      E=(n1matrix==traits*2 & last1M==0 & n2M==0) ## Maternal for the firsts traits and diag(R)
      F=(n1matrix==traits & n2M==0)               ## Only additive matrix and diag(R)
      if(is.na(F)){F=FALSE}
      #----------------------------------------------------------------------------------------#
      if(Uni==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        #--------------------------------------------------------------------------------------#
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+2))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dR=as.matrix(diag(mR))
        Ha=matrix(NA, nrow=nrow(arq), ncol=traits)
        for(i in 1:traits){
          aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+(arq[ ,dR[i]]))
          Ha[,i]=aux
          colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
        }
        cat("\nEstimating direct heritability: DONE\n")
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=paste("Vga",seq(1:traits),sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vR)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=cbind(arq, Ha)
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#    
      if(A==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        #--------------------------------------------------------------------------------------#
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+2))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dR=as.matrix(diag(mR))
        Ha=matrix(NA, nrow=nrow(arq), ncol=traits)
        Hm=matrix(NA, nrow=nrow(arq), ncol=traits)
        if(covAM==0){
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[1 ,dG[traits+i]])+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          for(i in 1:traits){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:traits), sep=''))
          }
          effectG=traits
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
        } else {
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[1 ,dG[traits+i]])+
                                (abs((abs((arq[ ,mG[i,traits+i]])))))+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          for(i in 1:traits){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:traits), sep=''))
          }
          effectG=nnzero(dG)
          CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:traits),sep="")))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)      
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        covR=mR[upper.tri(mR)]; covR[covR==0]=NA; covR=c(na.omit(covR));covR=sort(covR)
        matrixR=mR; diag(matrixR)=0; matrixR=which(matrixR!=0, arr.ind=T)        
        if (nrow(matrixR)>1){ matrixR=matrixR[order(matrixR[,1], matrixR[,2]), ] }
        covR=data.frame(id=covR, comp=paste("COV","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=matrix(NA, nrow=nrow(arq), ncol=((traits^2)-traits)/2)
        subCorE=mR
        colnames(arq)=c(1:ncol(arq))
        subCorE[lower.tri(subCorE)]=lowerTriangle(t(subCorE))
        aux=as.matrix(subCorE)
        for(n in 1:nrow(arq)){
          for(k in 1:ncol(arq)){
            for(i in 1:traits){
              for(j in 1:traits){
                if(subCorE[i,j]==as.numeric(colnames(arq)[k])){
                  aux[i,j]=arq[n, as.numeric(colnames(arq)[k])]
                }
              }
            }
          }
          aux=cov2cor(aux)
          aux=aux[lower.tri(aux)]
          CorrE[n,]=aux
          aux=as.matrix(subCorE)
        }
        CorrE=data.frame(id=seq(1:ncol(CorrE)), t(CorrE))
        CorrE=CorrE[!rowSums(CorrE==0)>=1, ]
        corE=data.frame(id=CorrE[,1], comp=paste("COR","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=merge(corE, CorrE, by=c("id", "id"), sort=FALSE)
        CorrE=t(CorrE)
        colnames(CorrE)=CorrE[2,]
        CorrE=as.matrix(CorrE)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vR,covG,covR)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=cbind(arq, CorrG, CorrE, Ha, Hm)
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
      if(B==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        #--------------------------------------------------------------------------------------#
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        post=rowSums(mG==0, 1)[1]
        pre=traits-post
        if(pre!=0){
          m=matrix(0, nrow=traits*2, ncol=traits*2)
          m[1:sum(rowSums(mG)>0),1:sum(rowSums(mG)>0)]=mG[rowSums(mG)>0, colSums(mG)>0]
          mG=m
        }
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+2))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dR=as.matrix(diag(mR))
        #--------------------------------------------------------------------------------------#
        if(covAM==0){
          Hapre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+(arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #--------------------------------------------------------------------------------------#
          Hm=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          effectG=traits
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
        } else {
          Hapre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((abs((arq[ ,mG[i,traits+i]])))))+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+(arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #--------------------------------------------------------------------------------------#
          Hm=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          effectG=nnzero(dG)
          CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[post+1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
        }
        #--------------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:pre),sep="")))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        covR=mR[upper.tri(mR)]; covR[covR==0]=NA; covR=c(na.omit(covR));covR=sort(covR)
        matrixR=mR; diag(matrixR)=0; matrixR=which(matrixR!=0, arr.ind=T)        
        if (nrow(matrixR)>1){ matrixR=matrixR[order(matrixR[,1], matrixR[,2]), ] }
        covR=data.frame(id=covR, comp=paste("COV","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=matrix(NA, nrow=nrow(arq), ncol=((traits^2)-traits)/2)
        subCorE=mR
        colnames(arq)=c(1:ncol(arq))
        subCorE[lower.tri(subCorE)]=lowerTriangle(t(subCorE))
        aux=as.matrix(subCorE)
        for(n in 1:nrow(arq)){
          for(k in 1:ncol(arq)){
            for(i in 1:traits){
              for(j in 1:traits){
                if(subCorE[i,j]==as.numeric(colnames(arq)[k])){
                  aux[i,j]=arq[n, as.numeric(colnames(arq)[k])]
                }
              }
            }
          }
          aux=cov2cor(aux)
          aux=aux[lower.tri(aux)]
          CorrE[n,]=aux
          aux=as.matrix(subCorE)
        }
        CorrE=data.frame(id=seq(1:ncol(CorrE)), t(CorrE))
        CorrE=CorrE[!rowSums(CorrE==0)>=1, ]
        corE=data.frame(id=CorrE[,1], comp=paste("COR","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=merge(corE, CorrE, by=c("id", "id"), sort=FALSE)
        CorrE=t(CorrE)
        colnames(CorrE)=CorrE[2,]
        CorrE=as.matrix(CorrE)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vR,covG,covR)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=data.frame(cbind(arq, CorrG, CorrE, Ha, Hm))
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
      if(C==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        #--------------------------------------------------------------------------------------#
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+2))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dR=as.matrix(diag(mR))
        Ha=matrix(NA, nrow=nrow(arq), ncol=traits)
        for(i in 1:traits){
          aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+(arq[ ,dR[i]]))
          Ha[,i]=aux
          colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
        }
        effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
        options(warn=-1)
        for(i in 1:nrow(arq)){
          aux=vec2sm(arq[i,dG[1]:rev(dG)[1]])
          aux=cov2cor(aux)
          aux=aux[lower.tri(aux)]
          CorrG[i,]=aux
        }
        options(warn=0)      
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=paste("Vga",seq(1:traits),sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))    
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        covR=mR[upper.tri(mR)]; covR[covR==0]=NA; covR=c(na.omit(covR));covR=sort(covR)
        matrixR=mR; diag(matrixR)=0; matrixR=which(matrixR!=0, arr.ind=T)        
        if (nrow(matrixR)>1){ matrixR=matrixR[order(matrixR[,1], matrixR[,2]), ] }
        covR=data.frame(id=covR, comp=paste("COV","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=matrix(NA, nrow=nrow(arq), ncol=((traits^2)-traits)/2)
        subCorE=mR
        colnames(arq)=c(1:ncol(arq))
        subCorE[lower.tri(subCorE)]=lowerTriangle(t(subCorE))
        aux=as.matrix(subCorE)
        for(n in 1:nrow(arq)){
          for(k in 1:ncol(arq)){
            for(i in 1:traits){
              for(j in 1:traits){
                if(subCorE[i,j]==as.numeric(colnames(arq)[k])){
                  aux[i,j]=arq[n, as.numeric(colnames(arq)[k])]
                }
              }
            }
          }
          aux=cov2cor(aux)
          aux=aux[lower.tri(aux)]
          CorrE[n,]=aux
          aux=as.matrix(subCorE)
        }
        CorrE=data.frame(id=seq(1:ncol(CorrE)), t(CorrE))
        CorrE=CorrE[!rowSums(CorrE==0)>=1, ]
        corE=data.frame(id=CorrE[,1], comp=paste("COR","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=merge(corE, CorrE, by=c("id", "id"), sort=FALSE)
        CorrE=t(CorrE)
        colnames(CorrE)=CorrE[2,]
        CorrE=as.matrix(CorrE)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vR,covG,covR)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=cbind(arq, CorrG, CorrE, Ha)
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
      if(D==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+2))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dR=as.matrix(diag(mR))
        Ha=matrix(NA, nrow=nrow(arq), ncol=traits)
        Hm=matrix(NA, nrow=nrow(arq), ncol=traits)
        if(covAM==0){
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[1 ,dG[traits+i]])+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          for(i in 1:traits){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:traits), sep=''))
          }
          effectG=traits
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
        } else {
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          for(i in 1:traits){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:traits), sep=''))
          }
          effectG=nnzero(dG)
          CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:traits),sep="")))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vR,covG)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=cbind(arq, CorrG, Ha, Hm)
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
      #----------------------------------------------------------------------------------------#
      if(E==TRUE){
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        post=rowSums(mG==0, 1)[1]
        if (post==0){pre=traits/2} else {pre=traits-post}
        if(pre!=0){
          m=matrix(0, nrow=traits*2, ncol=traits*2)
          m[1:sum(rowSums(mG)>0),1:sum(rowSums(mG)>0)]=mG[rowSums(mG)>0, colSums(mG)>0]
          mG=m
        }
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+2))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dR=as.matrix(diag(mR))
        #--------------------------------------------------------------------------------------#
        if(covAM==0){
          Hapre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+(arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #--------------------------------------------------------------------------------------#
          Hm=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          effectG=traits
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
        } else {
          Hapre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+(arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #--------------------------------------------------------------------------------------#
          Hm=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG[1:traits+pre, ])[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:pre),sep="")))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vR,covG)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=data.frame(cbind(arq, CorrG, Ha, Hm))
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
      if(F==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+2))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dR=as.matrix(diag(mR))
        Ha=matrix(NA, nrow=nrow(arq), ncol=traits)
        for(i in 1:traits){
          aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+(arq[ ,dR[i]]))
          Ha[,i]=aux
          colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
        }
        effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
        options(warn=-1)
        for(i in 1:nrow(arq)){
          aux=vec2sm(arq[i,dG[1]:rev(dG)[1]])
          aux=cov2cor(aux)
          aux=aux[lower.tri(aux)]
          CorrG[i,]=aux
        }
        options(warn=0)
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=paste("Vga",seq(1:traits),sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vR,covG)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=cbind(arq, CorrG, Ha)
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
    }
    #------------------------------------------------------------------------------------------#
    #------------------------------------------------------------------------------------------#
    # These scripts will run if it were found three matrices
    #------------------------------------------------------------------------------------------#
    if(length(grep("matrix",Matrix))==3){
      #------------------------------------------------------------------------------------------#
      cat("\n")
      cat("Reading matrices from Bayesian analyses...\n")
      #------------------------------------------------------------------------------------------#
      cat("\n")  
      cat("You are summarizing the results from multi-trait analysis...")
      cat("\n")  
      #------------------------------------------------------------------------------------------#
      n1matrix=as.matrix(strsplit(as.character(Matrix[2]), " ")[[1]])
      mode(n1matrix)="numeric"
      n1matrix=c(length(na.omit(n1matrix)))
      mode(n1matrix)="numeric"
      #------------------------------------------------------------------------------------------#
      eff_G=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
      eff_G=as.matrix(diag(eff_G))
      eff_G=nnzero(eff_G)
      #------------------------------------------------------------------------------------------#
      cat("1st (co)variance matrix has size of", eff_G,"x",eff_G)
      cat("\n")
      nMatrix=grep("matrix",Matrix)
      #----------------------------------------------------------------------------------------#
      n2matrix=(nMatrix[3]-1)-(nMatrix[2])
      cat("2nd (co)variance matrix has size of", n2matrix,"x",n2matrix)
      cat("\n")
      #----------------------------------------------------------------------------------------#
      n3matrix=(nMatrix[3]-1)-(nMatrix[2])
      cat("3rd (co)variance matrix has size of", n3matrix,"x",n3matrix)
      cat("\n")
      cat("\n")
      #----------------------------------------------------------------------------------------#
      last1M=as.matrix(strsplit(as.character(Matrix[2]), " ")[[1]])
      mode(last1M)="numeric"
      last1M=c(tail(na.omit(last1M), n=1))
      #----------------------------------------------------------------------------------------#
      n1M=as.matrix(strsplit(as.character(Matrix[traits+2]), " ")[[1]])
      mode(n1M)="numeric"
      n1M=rev(na.omit(n1M))[1] # Reading the last value
      #----------------------------------------------------------------------------------------#
      n2M=as.matrix(strsplit(as.character(Matrix[n1matrix+3]), " ")[[1]])
      mode(n2M)="numeric"
      n2M=na.omit(n2M)[2]
      #----------------------------------------------------------------------------------------#
      last2M=as.matrix(strsplit(as.character(Matrix[n1matrix+3]), " ")[[1]])
      mode(last2M)="numeric"
      last2M=c(tail(na.omit(last2M), n=1))
      #----------------------------------------------------------------------------------------#
      n3M=as.matrix(strsplit(as.character(Matrix[n1matrix+n2matrix+4]), " ")[[1]])
      mode(n3M)="numeric"
      n3M=na.omit(n3M)[2]
      #----------------------------------------------------------------------------------------#
      # A = A + M + MPE + R --> M and MPE for all traits
      # B = A + M + MPE + R --> M for all traits, but MPE only for pre weaning traits
      # Ba= A + M + MPE + R --> M for all traits, but MPE is diagonal
      # C = A + M + MPE + R --> M and MPE only for pre weaning traits
      # D = A + M + MPE + R --> M only for pre weaning traits, but MPE is diagonal
      #----------------------------------------------------------------------------------------#
      if(covAM==0){
        A =(n1matrix==traits*2 & last1M==0 & last2M!=0 & n2M!=0 & n3M!=0)
        B =(n1matrix==traits*2 & last1M==0 & last2M==0 & n2M!=0 & n3M!=0)
        Ba=(n1matrix==traits*2 & last1M==0 & last2M==0 & n2M==0 & n3M!=0)
        C =(n1matrix==traits*2 & last1M==0 & last2M==0 & n1M==0 & n2M!=0 & n3M!=0)
        D =(n1matrix==traits*2 & last1M==0 & last2M==0 & n1M==0 & n2M==0 & n3M!=0)
      } else {
        A =(n1matrix==traits*2 & last1M!=0 & last2M!=0 & n2M!=0 & n3M!=0)
        B =(n1matrix==traits*2 & last1M!=0 & last2M==0 & n2M!=0 & n3M!=0)
        Ba=(n1matrix==traits*2 & last1M!=0 & last2M==0 & n2M==0 & n3M!=0)
        C =(n1matrix==traits*2 & last1M==0 & last2M==0 & n2M!=0 & n3M!=0)
        D =(n1matrix==traits*2 & last1M==0 & last2M==0 & n2M==0 & n3M!=0)
      }
      #----------------------------------------------------------------------------------------#
      # E = A + M + MPE + R --> M and MPE for all traits and diag(R)
      # F = A + M + MPE + R --> M for all traits but MPE only for pre weaning traits and diag(R)
      # Fa= A + M + MPE + R --> M for all traits, but MPE and R are diagonal
      # G = A + M + MPE + R --> M and MPE only for pre weaning traits and diag(R)
      # H = A + M + MPE + R --> M only for pre weaning traits, diag(MPE) and diag(R)
      #----------------------------------------------------------------------------------------#
      if(covAM==0){
        E =(n1matrix==traits*2 & last1M==0 & last2M!=0 & n2M!=0 & n3M==0)
        F =(n1matrix==traits*2 & last1M==0 & last2M==0 & n2M!=0 & n3M==0)
        Fa=(n1matrix==traits*2 & last1M==0 & last2M==0 & n2M==0 & n3M==0)
        G =(n1matrix==traits*2 & last1M==0 & last2M==0 & n1M==0 & n2M!=0 & n3M==0)
        H =(n1matrix==traits*2 & last1M==0 & last2M==0 & n1M==0 & n2M==0 & n3M==0)
      } else {
        E =(n1matrix==traits*2 & last1M!=0 & last2M!=0 & n2M!=0 & n3M==0)
        F =(n1matrix==traits*2 & last1M!=0 & last2M==0 & n2M!=0 & n3M==0)
        Fa=(n1matrix==traits*2 & last1M!=0 & last2M==0 & n2M==0 & n3M==0)
        G =(n1matrix==traits*2 & last1M==0 & last2M==0 & n2M!=0 & n3M==0)
        H =(n1matrix==traits*2 & last1M==0 & last2M==0 & n2M==0 & n3M==0)
      }
      #----------------------------------------------------------------------------------------#
      # R1 = A + PE + R --> All matrix
      # R2 = A + PE + R --> A matrix, PE diagonal and R matrix
      # R3 = A + PE + R --> A matrix, PE and R are diagonals
      #----------------------------------------------------------------------------------------#
      R1 =(n1matrix==traits & n2M!=0 & n3M!=0)
      R2 =(n1matrix==traits & n2M==0 & n3M!=0)
      R3 =(n1matrix==traits & n2M==0 & n3M==0)
      #----------------------------------------------------------------------------------------#
      if(A==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+traits+3))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dPe=as.matrix(diag(mPe))
        dR=as.matrix(diag(mR))
        #--------------------------------------------------------------------------------------#
        Ha=matrix(NA, nrow=nrow(arq), ncol=traits)
        Hm=matrix(NA, nrow=nrow(arq), ncol=traits)
        if(covAM==0){
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          for(i in 1:traits){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:traits), sep=''))
          }
          effectG=traits
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
          #--------------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=traits)
          for(i in 1:traits){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:traits), sep=''))
          }
        } else {
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          for(i in 1:traits){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:traits), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          effectG=nnzero(dG)
          CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
          #--------------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=traits)
          for(i in 1:traits){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (abs((arq[ ,mG[i,traits+i]])))+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:traits), sep=''))
          }
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:traits),sep="")))
        vPe=dPe; vPe[vPe==0]=NA; vPe=na.omit(vPe); vPe=data.frame(id=vPe, comp=paste("Vmpe", seq(1:nrow(vPe)), sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        covPe=mPe[upper.tri(mPe)]; covPe[covPe==0]=NA; covPe=c(na.omit(covPe));covPe=sort(covPe)
        matrixPe=mPe; diag(matrixPe)=0; matrixPe=which(matrixPe!=0,arr.ind=T)        
        if (nrow(matrixPe)>1){ matrixPe=matrixPe[order(matrixPe[,1], matrixPe[,2]), ] }
        covPe=data.frame(id=covPe, comp=paste("COV","mpe",matrixPe[,1],"mpe",matrixPe[,2], sep=""))
        #--------------------------------------------------------------------------------------#
        covR=mR[upper.tri(mR)]; covR[covR==0]=NA; covR=c(na.omit(covR));covR=sort(covR)
        matrixR=mR; diag(matrixR)=0; matrixR=which(matrixR!=0, arr.ind=T)        
        if (nrow(matrixR)>1){ matrixR=matrixR[order(matrixR[,1], matrixR[,2]), ] }
        covR=data.frame(id=covR, comp=paste("COV","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=matrix(NA, nrow=nrow(arq), ncol=((traits^2)-traits)/2)
        subCorE=mR
        colnames(arq)=c(1:ncol(arq))
        subCorE[lower.tri(subCorE)]=lowerTriangle(t(subCorE))
        aux=as.matrix(subCorE)
        for(n in 1:nrow(arq)){
          for(k in 1:ncol(arq)){
            for(i in 1:traits){
              for(j in 1:traits){
                if(subCorE[i,j]==as.numeric(colnames(arq)[k])){
                  aux[i,j]=arq[n, as.numeric(colnames(arq)[k])]
                }
              }
            }
          }
          aux=cov2cor(aux)
          aux=aux[lower.tri(aux)]
          CorrE[n,]=aux
          aux=as.matrix(subCorE)
        }
        CorrE=data.frame(id=seq(1:ncol(CorrE)), t(CorrE))
        CorrE=CorrE[!rowSums(CorrE==0)>=1, ]
        corE=data.frame(id=CorrE[,1], comp=paste("COR","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=merge(corE, CorrE, by=c("id", "id"), sort=FALSE)
        CorrE=t(CorrE)
        colnames(CorrE)=CorrE[2,]
        CorrE=as.matrix(CorrE)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vPe,vR,covG,covPe,covR)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=cbind(arq, CorrG, CorrE, Ha, Hm, c2)
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
      if(B==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        post=rowSums(mG==0, 1)[1]
        if (post==0){pre=traits/2} else {pre=traits-post}
        if(pre!=0){
          m=matrix(0, nrow=traits*2, ncol=traits*2)
          m[1:sum(rowSums(mG)>0),1:sum(rowSums(mG)>0)]=mG[rowSums(mG)>0, colSums(mG)>0]
          mG=m
        }
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        if(pre!=0){
          pe=matrix(0, nrow=traits, ncol=traits)
          pe[1:sum(rowSums(mPe)>0),1:sum(rowSums(mPe)>0)]=mPe[rowSums(mPe)>0, colSums(mPe)>0]
          mPe=pe
        }
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+traits+3))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dPe=as.matrix(diag(mPe))
        dR=as.matrix(diag(mR))
        #--------------------------------------------------------------------------------------#
        Hapre=matrix(NA, nrow=nrow(arq), ncol=pre)
        if(covAM==0){
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          if (post==0){post=pre} else {post=post}
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+
                                    (arq[ ,dG[traits+pre+i]])+
                                    (arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          Hmpre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hmpre[,i]=aux
            colnames(Hmpre)=c(paste('Hmpre', seq(1:pre), sep=''))
          }
          Hmpost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[traits+pre+i]]/(arq[ ,dG[pre+i]]+
                                           (arq[ ,dG[traits+pre+i]])+
                                           (arq[ ,dR[pre+i]]))
            Hmpost[,i]=aux
            colnames(Hmpost)=c(paste('Hmpost', seq(1:post), sep=''))
          }
          Hm=cbind(Hmpre, Hmpost)
          #---------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
          #--------------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        } else {
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          if (post==0){post=pre} else {post=post}
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+
                                    (arq[ ,dG[traits+pre+i]])+
                                    (abs(arq[ ,mG[pre+i,traits+pre+i]]))+
                                    (arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          Hmpre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hmpre[,i]=aux
            colnames(Hmpre)=c(paste('Hmpre', seq(1:pre), sep=''))
          }
          Hmpost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[traits+pre+i]]/(arq[ ,dG[pre+i]]+
                                           (arq[ ,dG[traits+pre+i]])+
                                           (abs(arq[ ,mG[pre+i,traits+pre+i]]))+
                                           (arq[ ,dR[pre+i]]))
            Hmpost[,i]=aux
            colnames(Hmpost)=c(paste('Hmpost', seq(1:post), sep=''))
          }
          Hm=cbind(Hmpre, Hmpost)
          #--------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
          #--------------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (abs((arq[ ,mG[i,traits+i]])))+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:traits),sep="")))
        vPe=dPe; vPe[vPe==0]=NA; vPe=na.omit(vPe); vPe=data.frame(id=vPe, comp=paste("Vmpe", seq(1:nrow(vPe)), sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        covPe=mPe[upper.tri(mPe)]; covPe[covPe==0]=NA; covPe=c(na.omit(covPe));covPe=sort(covPe)
        matrixPe=mPe; diag(matrixPe)=0; matrixPe=which(matrixPe!=0,arr.ind=T)        
        if (nrow(matrixPe)>1){ matrixPe=matrixPe[order(matrixPe[,1], matrixPe[,2]), ] }
        covPe=data.frame(id=covPe, comp=paste("COV","mpe",matrixPe[,1],"mpe",matrixPe[,2], sep=""))
        #--------------------------------------------------------------------------------------#
        covR=mR[upper.tri(mR)]; covR[covR==0]=NA; covR=c(na.omit(covR));covR=sort(covR)
        matrixR=mR; diag(matrixR)=0; matrixR=which(matrixR!=0, arr.ind=T)        
        if (nrow(matrixR)>1){ matrixR=matrixR[order(matrixR[,1], matrixR[,2]), ] }
        covR=data.frame(id=covR, comp=paste("COV","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=matrix(NA, nrow=nrow(arq), ncol=((traits^2)-traits)/2)
        subCorE=mR
        colnames(arq)=c(1:ncol(arq))
        subCorE[lower.tri(subCorE)]=lowerTriangle(t(subCorE))
        aux=as.matrix(subCorE)
        for(n in 1:nrow(arq)){
          for(k in 1:ncol(arq)){
            for(i in 1:traits){
              for(j in 1:traits){
                if(subCorE[i,j]==as.numeric(colnames(arq)[k])){
                  aux[i,j]=arq[n, as.numeric(colnames(arq)[k])]
                }
              }
            }
          }
          aux=cov2cor(aux)
          aux=aux[lower.tri(aux)]
          CorrE[n,]=aux
          aux=as.matrix(subCorE)
        }
        CorrE=data.frame(id=seq(1:ncol(CorrE)), t(CorrE))
        CorrE=CorrE[!rowSums(CorrE==0)>=1, ]
        corE=data.frame(id=CorrE[,1], comp=paste("COR","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=merge(corE, CorrE, by=c("id", "id"), sort=FALSE)
        CorrE=t(CorrE)
        colnames(CorrE)=CorrE[2,]
        CorrE=as.matrix(CorrE)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vPe,vR,covG,covPe,covR)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=data.frame(cbind(arq, CorrG, CorrE, Ha, Hm, c2))
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
      }
      #----------------------------------------------------------------------------------------#
      if(Ba==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        post=rowSums(mG==0, 1)[1]
        if (post==0){pre=traits/2} else {pre=traits-post}
        if(pre!=0){
          pe=matrix(0, nrow=traits, ncol=traits)
          pe[1:sum(rowSums(mPe)>0),1:sum(rowSums(mPe)>0)]=mPe[rowSums(mPe)>0, colSums(mPe)>0]
          mPe=pe
        }
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+traits+3))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dPe=as.matrix(diag(mPe))
        dR=as.matrix(diag(mR))
        #--------------------------------------------------------------------------------------#
        Hapre=matrix(NA, nrow=nrow(arq), ncol=pre)
        if(covAM==0){
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          if (post==0){post=pre} else {post=post}
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+
                                    (arq[ ,dG[traits+pre+i]])+
                                    (arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          Hmpre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hmpre[,i]=aux
            colnames(Hmpre)=c(paste('Hmpre', seq(1:pre), sep=''))
          }
          Hmpost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[traits+pre+i]]/(arq[ ,dG[pre+i]]+
                                           (arq[ ,dG[traits+pre+i]])+
                                           (arq[ ,dR[pre+i]]))
            Hmpost[,i]=aux
            colnames(Hmpost)=c(paste('Hmpost', seq(1:post), sep=''))
          }
          Hm=cbind(Hmpre, Hmpost)
          #---------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
          #--------------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        } else {
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          if (post==0){post=pre} else {post=post}
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+
                                    (arq[ ,dG[traits+pre+i]])+
                                    (abs(arq[ ,mG[pre+i,traits+pre+i]]))+
                                    (arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          Hmpre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hmpre[,i]=aux
            colnames(Hmpre)=c(paste('Hmpre', seq(1:pre), sep=''))
          }
          Hmpost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[traits+pre+i]]/(arq[ ,dG[pre+i]]+
                                           (arq[ ,dG[traits+pre+i]])+
                                           (abs(arq[ ,mG[pre+i,traits+pre+i]]))+
                                           (arq[ ,dR[pre+i]]))
            Hmpost[,i]=aux
            colnames(Hmpost)=c(paste('Hmpost', seq(1:post), sep=''))
          }
          Hm=cbind(Hmpre, Hmpost)
          #--------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
          #--------------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (abs((arq[ ,mG[i,traits+i]])))+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:traits),sep="")))
        vPe=dPe; vPe[vPe==0]=NA; vPe=na.omit(vPe); vPe=data.frame(id=vPe, comp=paste("Vmpe", seq(1:nrow(vPe)), sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        covR=mR[upper.tri(mR)]; covR[covR==0]=NA; covR=c(na.omit(covR));covR=sort(covR)
        matrixR=mR; diag(matrixR)=0; matrixR=which(matrixR!=0, arr.ind=T)        
        if (nrow(matrixR)>1){ matrixR=matrixR[order(matrixR[,1], matrixR[,2]), ] }
        covR=data.frame(id=covR, comp=paste("COV","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=matrix(NA, nrow=nrow(arq), ncol=((traits^2)-traits)/2)
        subCorE=mR
        colnames(arq)=c(1:ncol(arq))
        subCorE[lower.tri(subCorE)]=lowerTriangle(t(subCorE))
        aux=as.matrix(subCorE)
        for(n in 1:nrow(arq)){
          for(k in 1:ncol(arq)){
            for(i in 1:traits){
              for(j in 1:traits){
                if(subCorE[i,j]==as.numeric(colnames(arq)[k])){
                  aux[i,j]=arq[n, as.numeric(colnames(arq)[k])]
                }
              }
            }
          }
          aux=cov2cor(aux)
          aux=aux[lower.tri(aux)]
          CorrE[n,]=aux
          aux=as.matrix(subCorE)
        }
        CorrE=data.frame(id=seq(1:ncol(CorrE)), t(CorrE))
        CorrE=CorrE[!rowSums(CorrE==0)>=1, ]
        corE=data.frame(id=CorrE[,1], comp=paste("COR","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=merge(corE, CorrE, by=c("id", "id"), sort=FALSE)
        CorrE=t(CorrE)
        colnames(CorrE)=CorrE[2,]
        CorrE=as.matrix(CorrE)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vPe,vR,covG,covR)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=data.frame(cbind(arq, CorrG, CorrE, Ha, Hm, c2))
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
      }
      #----------------------------------------------------------------------------------------#
      if(C==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        post=rowSums(mG==0, 1)[1]
        if (post==0){pre=traits/2} else {pre=traits-post}
        if(pre!=0){
          m=matrix(0, nrow=traits*2, ncol=traits*2)
          m[1:sum(rowSums(mG)>0),1:sum(rowSums(mG)>0)]=mG[rowSums(mG)>0, colSums(mG)>0]
          mG=m
        }
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T)
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        if(pre!=0){
          pe=matrix(0, nrow=traits, ncol=traits)
          pe[1:sum(rowSums(mPe)>0),1:sum(rowSums(mPe)>0)]=mPe[rowSums(mPe)>0, colSums(mPe)>0]
          mPe=pe
        }
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+traits+3))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dPe=as.matrix(diag(mPe))
        dR=as.matrix(diag(mR))
        #----------------------------------------------------------------------------------#
        Hapre=matrix(NA, nrow=nrow(arq), ncol=pre)
        if(covAM==0){
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+(arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #----------------------------------------------------------------------------------#
          Hm=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:pre), sep=''))
          }
          #---------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
          #----------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        } else {
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+(arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #----------------------------------------------------------------------------------#
          Hm=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:pre), sep=''))
          }
          #----------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[post+1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          #----------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (abs((arq[ ,mG[i,traits+i]])))+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:pre),sep="")))
        vPe=dPe; vPe[vPe==0]=NA; vPe=na.omit(vPe); vPe=data.frame(id=vPe, comp=paste("Vmpe", seq(1:nrow(vPe)), sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        covPe=mPe[upper.tri(mPe)]; covPe[covPe==0]=NA; covPe=c(na.omit(covPe));covPe=sort(covPe)
        matrixPe=mPe; diag(matrixPe)=0; matrixPe=which(matrixPe!=0,arr.ind=T)        
        if (nrow(matrixPe)>1){ matrixPe=matrixPe[order(matrixPe[,1], matrixPe[,2]), ] }
        covPe=data.frame(id=covPe, comp=paste("COV","mpe",matrixPe[,1],"mpe",matrixPe[,2], sep=""))
        #--------------------------------------------------------------------------------------#
        covR=mR[upper.tri(mR)]; covR[covR==0]=NA; covR=c(na.omit(covR));covR=sort(covR)
        matrixR=mR; diag(matrixR)=0; matrixR=which(matrixR!=0, arr.ind=T)        
        if (nrow(matrixR)>1){ matrixR=matrixR[order(matrixR[,1], matrixR[,2]), ] }
        covR=data.frame(id=covR, comp=paste("COV","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=matrix(NA, nrow=nrow(arq), ncol=((traits^2)-traits)/2)
        subCorE=mR
        colnames(arq)=c(1:ncol(arq))
        subCorE[lower.tri(subCorE)]=lowerTriangle(t(subCorE))
        aux=as.matrix(subCorE)
        for(n in 1:nrow(arq)){
          for(k in 1:ncol(arq)){
            for(i in 1:traits){
              for(j in 1:traits){
                if(subCorE[i,j]==as.numeric(colnames(arq)[k])){
                  aux[i,j]=arq[n, as.numeric(colnames(arq)[k])]
                }
              }
            }
          }
          aux=cov2cor(aux)
          aux=aux[lower.tri(aux)]
          CorrE[n,]=aux
          aux=as.matrix(subCorE)
        }
        CorrE=data.frame(id=seq(1:ncol(CorrE)), t(CorrE))
        CorrE=CorrE[!rowSums(CorrE==0)>=1, ]
        corE=data.frame(id=CorrE[,1], comp=paste("COR","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=merge(corE, CorrE, by=c("id", "id"), sort=FALSE)
        CorrE=t(CorrE)
        colnames(CorrE)=CorrE[2,]
        CorrE=as.matrix(CorrE)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vPe,vR,covG,covPe,covR)
        #--------------------------------------------------------------------------------------#
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=data.frame(cbind(arq, CorrG, CorrE, Ha, Hm, c2))
        arq=as.matrix(arq); mode(arq)="numeric"
        #----------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #----------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
      if(D==TRUE){
        #----------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        post=rowSums(mG==0, 1)[1]
        if (post==0){pre=traits/2} else {pre=traits-post}
        if(pre!=0){
          m=matrix(0, nrow=traits*2, ncol=traits*2)
          m[1:sum(rowSums(mG)>0),1:sum(rowSums(mG)>0)]=mG[rowSums(mG)>0, colSums(mG)>0]
          mG=m
        }
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        if(pre!=0){
          pe=matrix(0, nrow=traits, ncol=traits)
          pe[1:sum(rowSums(mPe)>0),1:sum(rowSums(mPe)>0)]=mPe[rowSums(mPe)>0, colSums(mPe)>0]
          mPe=pe
        }
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+traits+3))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dPe=as.matrix(diag(mPe))
        dR=as.matrix(diag(mR))
        #----------------------------------------------------------------------------------#
        Hapre=matrix(NA, nrow=nrow(arq), ncol=pre)
        if(covAM==0){
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+(arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #----------------------------------------------------------------------------------#
          Hm=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:pre), sep=''))
          }
          #---------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
          #----------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        } else {
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+(arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #----------------------------------------------------------------------------------#
          Hm=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:pre), sep=''))
          }
          #----------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[post+1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          #----------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (abs((arq[ ,mG[i,traits+i]])))+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:pre),sep="")))
        vPe=dPe; vPe[vPe==0]=NA; vPe=na.omit(vPe); vPe=data.frame(id=vPe, comp=paste("Vmpe", seq(1:nrow(vPe)), sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        covR=mR[upper.tri(mR)]; covR[covR==0]=NA; covR=c(na.omit(covR));covR=sort(covR)
        matrixR=mR; diag(matrixR)=0; matrixR=which(matrixR!=0, arr.ind=T)        
        if (nrow(matrixR)>1){ matrixR=matrixR[order(matrixR[,1], matrixR[,2]), ] }
        covR=data.frame(id=covR, comp=paste("COV","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=matrix(NA, nrow=nrow(arq), ncol=((traits^2)-traits)/2)
        subCorE=mR
        colnames(arq)=c(1:ncol(arq))
        subCorE[lower.tri(subCorE)]=lowerTriangle(t(subCorE))
        aux=as.matrix(subCorE)
        for(n in 1:nrow(arq)){
          for(k in 1:ncol(arq)){
            for(i in 1:traits){
              for(j in 1:traits){
                if(subCorE[i,j]==as.numeric(colnames(arq)[k])){
                  aux[i,j]=arq[n, as.numeric(colnames(arq)[k])]
                }
              }
            }
          }
          aux=cov2cor(aux)
          aux=aux[lower.tri(aux)]
          CorrE[n,]=aux
          aux=as.matrix(subCorE)
        }
        CorrE=data.frame(id=seq(1:ncol(CorrE)), t(CorrE))
        CorrE=CorrE[!rowSums(CorrE==0)>=1, ]
        corE=data.frame(id=CorrE[,1], comp=paste("COR","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=merge(corE, CorrE, by=c("id", "id"), sort=FALSE)
        CorrE=t(CorrE)
        colnames(CorrE)=CorrE[2,]
        CorrE=as.matrix(CorrE)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vPe,vR,covG,covR)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=data.frame(cbind(arq, CorrG, CorrE, Ha, Hm, c2))
        arq=as.matrix(arq); mode(arq)="numeric"
        #----------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #----------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
      }
      #----------------------------------------------------------------------------------------#
      if(E==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+traits+3))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dPe=as.matrix(diag(mPe))
        dR=as.matrix(diag(mR))
        Ha=matrix(NA, nrow=nrow(arq), ncol=traits)
        if(covAM==0){
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          Hm=matrix(NA, nrow=nrow(arq), ncol=traits)
          for(i in 1:traits){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:traits), sep=''))
          }
          #---------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
          c2=matrix(NA, nrow=nrow(arq), ncol=traits)
          for(i in 1:traits){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:traits), sep=''))
          }
        } else {
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          Hm=matrix(NA, nrow=nrow(arq), ncol=traits)
          for(i in 1:traits){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hm[,i]=aux
            colnames(Hm)=c(paste('hm', seq(1:traits), sep=''))
          }
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
          c2=matrix(NA, nrow=nrow(arq), ncol=traits)
          for(i in 1:traits){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (abs((arq[ ,mG[i,traits+i]])))+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:traits), sep=''))
          }
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:traits),sep="")))
        vPe=dPe; vPe[vPe==0]=NA; vPe=na.omit(vPe); vPe=data.frame(id=vPe, comp=paste("Vmpe", seq(1:nrow(vPe)), sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        covPe=mPe[upper.tri(mPe)]; covPe[covPe==0]=NA; covPe=c(na.omit(covPe));covPe=sort(covPe)
        matrixPe=mPe; diag(matrixPe)=0; matrixPe=which(matrixPe!=0,arr.ind=T)        
        if (nrow(matrixPe)>1){ matrixPe=matrixPe[order(matrixPe[,1], matrixPe[,2]), ] }
        covPe=data.frame(id=covPe, comp=paste("COV","mpe",matrixPe[,1],"mpe",matrixPe[,2], sep=""))
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vPe,vR,covG,covPe)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=cbind(arq, CorrG, Ha, Hm, c2)
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
      if(F==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        if(pre!=0){
          pe=matrix(0, nrow=traits, ncol=traits)
          pe[1:sum(rowSums(mPe)>0),1:sum(rowSums(mPe)>0)]=mPe[rowSums(mPe)>0, colSums(mPe)>0]
          mPe=pe
        }
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+traits+3))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dPe=as.matrix(diag(mPe))
        dR=as.matrix(diag(mR))
        post=rowSums(mG==0, 1)[1]
        if (post==0){pre=traits/2} else {pre=traits-post}
        #--------------------------------------------------------------------------------------#
        Hapre=matrix(NA, nrow=nrow(arq), ncol=pre)
        if(covAM==0){
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          if (post==0){post=pre} else {post=post}
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+
                                    (arq[ ,dG[traits+pre+i]])+
                                    (arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #--------------------------------------------------------------------------------------#
          Hmpre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hmpre[,i]=aux
            colnames(Hmpre)=c(paste('Hmpre', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          Hmpost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[traits+pre+i]]/(arq[ ,dG[pre+i]]+
                                           (arq[ ,dG[traits+pre+i]])+
                                           (arq[ ,dR[pre+i]]))
            Hmpost[,i]=aux
            colnames(Hmpost)=c(paste('Hmpost', seq(1:post), sep=''))
          }
          Hm=cbind(Hmpre, Hmpost)
          #---------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
          #--------------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        } else {
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          if (post==0){post=pre} else {post=post}
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+
                                    (arq[ ,dG[traits+pre+i]])+
                                    (abs(arq[ ,mG[pre+i,traits+pre+i]]))+
                                    (arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #--------------------------------------------------------------------------------------#
          Hmpre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hmpre[,i]=aux
            colnames(Hmpre)=c(paste('Hmpre', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          Hmpost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[traits+pre+i]]/(arq[ ,dG[pre+i]]+
                                           (arq[ ,dG[traits+pre+i]])+
                                           (abs(arq[ ,mG[pre+i,traits+pre+i]]))+
                                           (arq[ ,dR[pre+i]]))
            Hmpost[,i]=aux
            colnames(Hmpost)=c(paste('Hmpost', seq(1:post), sep=''))
          }
          Hm=cbind(Hmpre, Hmpost)
          #--------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
          #--------------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (abs((arq[ ,mG[i,traits+i]])))+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:traits),sep="")))
        vPe=dPe; vPe[vPe==0]=NA; vPe=na.omit(vPe); vPe=data.frame(id=vPe, comp=paste("Vmpe", seq(1:nrow(vPe)), sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        covPe=mPe[upper.tri(mPe)]; covPe[covPe==0]=NA; covPe=c(na.omit(covPe));covPe=sort(covPe)
        matrixPe=mPe; diag(matrixPe)=0; matrixPe=which(matrixPe!=0,arr.ind=T)        
        if (nrow(matrixPe)>1){ matrixPe=matrixPe[order(matrixPe[,1], matrixPe[,2]), ] }
        covPe=data.frame(id=covPe, comp=paste("COV","mpe",matrixPe[,1],"mpe",matrixPe[,2], sep=""))
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vPe,vR,covG,covPe)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=data.frame(cbind(arq, CorrG, Ha, Hm, c2))
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
      if(Fa==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        if(pre!=0){
          pe=matrix(0, nrow=traits, ncol=traits)
          pe[1:sum(rowSums(mPe)>0),1:sum(rowSums(mPe)>0)]=mPe[rowSums(mPe)>0, colSums(mPe)>0]
          mPe=pe
        }
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+traits+3))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dPe=as.matrix(diag(mPe))
        dR=as.matrix(diag(mR))
        post=rowSums(mG==0, 1)[1]
        if (post==0){pre=traits/2} else {pre=traits-post}
        #--------------------------------------------------------------------------------------#
        Hapre=matrix(NA, nrow=nrow(arq), ncol=pre)
        if(covAM==0){
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          if (post==0){post=pre} else {post=post}
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+
                                    (arq[ ,dG[traits+pre+i]])+
                                    (arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #--------------------------------------------------------------------------------------#
          Hmpre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hmpre[,i]=aux
            colnames(Hmpre)=c(paste('Hmpre', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          Hmpost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[traits+pre+i]]/(arq[ ,dG[pre+i]]+
                                           (arq[ ,dG[traits+pre+i]])+
                                           (arq[ ,dR[pre+i]]))
            Hmpost[,i]=aux
            colnames(Hmpost)=c(paste('Hmpost', seq(1:post), sep=''))
          }
          Hm=cbind(Hmpre, Hmpost)
          #---------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
          #--------------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        } else {
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          if (post==0){post=pre} else {post=post}
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+
                                    (arq[ ,dG[traits+pre+i]])+
                                    (abs(arq[ ,mG[pre+i,traits+pre+i]]))+
                                    (arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #--------------------------------------------------------------------------------------#
          Hmpre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hmpre[,i]=aux
            colnames(Hmpre)=c(paste('Hmpre', seq(1:pre), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          Hmpost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[traits+pre+i]]/(arq[ ,dG[pre+i]]+
                                           (arq[ ,dG[traits+pre+i]])+
                                           (abs(arq[ ,mG[pre+i,traits+pre+i]]))+
                                           (arq[ ,dR[pre+i]]))
            Hmpost[,i]=aux
            colnames(Hmpost)=c(paste('Hmpost', seq(1:post), sep=''))
          }
          Hm=cbind(Hmpre, Hmpost)
          #--------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
          #--------------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (abs((arq[ ,mG[i,traits+i]])))+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:traits),sep="")))
        vPe=dPe; vPe[vPe==0]=NA; vPe=na.omit(vPe); vPe=data.frame(id=vPe, comp=paste("Vmpe", seq(1:nrow(vPe)), sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vPe,vR,covG)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=data.frame(cbind(arq, CorrG, Ha, Hm, c2))
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
      if(G==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        post=rowSums(mG==0, 1)[1]
        if (post==0){pre=traits/2} else {pre=traits-post}
        if(pre!=0){
          m=matrix(0, nrow=traits*2, ncol=traits*2)
          m[1:sum(rowSums(mG)>0),1:sum(rowSums(mG)>0)]=mG[rowSums(mG)>0, colSums(mG)>0]
          mG=m
        }
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        if(pre!=0){
          pe=matrix(0, nrow=traits, ncol=traits)
          pe[1:sum(rowSums(mPe)>0),1:sum(rowSums(mPe)>0)]=mPe[rowSums(mPe)>0, colSums(mPe)>0]
          mPe=pe
        }
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+traits+3))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dPe=as.matrix(diag(mPe))
        dR=as.matrix(diag(mR))
        #----------------------------------------------------------------------------------#
        Hapre=matrix(NA, nrow=nrow(arq), ncol=pre)
        if(covAM==0){
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          #----------------------------------------------------------------------------------#
          if (post==0){post=pre} else {post=post}
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+(arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #----------------------------------------------------------------------------------#
          Hmpre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hmpre[,i]=aux
            colnames(Hmpre)=c(paste('Hmpre', seq(1:pre), sep=''))
          }
          Hm=Hmpre
          #---------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
          #----------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        } else {
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          #----------------------------------------------------------------------------------#
          if (post==0){post=pre} else {post=post}
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+(arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #----------------------------------------------------------------------------------#
          Hmpre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hmpre[,i]=aux
            colnames(Hmpre)=c(paste('Hmpre', seq(1:pre), sep=''))
          }
          Hm=Hmpre
          #----------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG[1:traits+pre, ])[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
          #----------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (abs((arq[ ,mG[i,traits+i]])))+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:pre),sep="")))
        vPe=dPe; vPe[vPe==0]=NA; vPe=na.omit(vPe); vPe=data.frame(id=vPe, comp=paste("Vmpe", seq(1:nrow(vPe)), sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        covPe=mPe[upper.tri(mPe)]; covPe[covPe==0]=NA; covPe=c(na.omit(covPe));covPe=sort(covPe)
        matrixPe=mPe; diag(matrixPe)=0; matrixPe=which(matrixPe!=0,arr.ind=T)        
        if (nrow(matrixPe)>1){ matrixPe=matrixPe[order(matrixPe[,1], matrixPe[,2]), ] }
        covPe=data.frame(id=covPe, comp=paste("COV","mpe",matrixPe[,1],"mpe",matrixPe[,2], sep=""))
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vPe,vR,covG,covPe)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=data.frame(cbind(arq, CorrG, Ha, Hm, c2))
        arq=as.matrix(arq); mode(arq)="numeric"
        #----------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #----------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
      if(H==TRUE){
        #----------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        post=rowSums(mG==0, 1)[1]
        if (post==0){pre=traits/2} else {pre=traits-post}
        if(pre!=0){
          m=matrix(0, nrow=traits*2, ncol=traits*2)
          m[1:sum(rowSums(mG)>0),1:sum(rowSums(mG)>0)]=mG[rowSums(mG)>0, colSums(mG)>0]
          mG=m
        }
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        if(pre!=0){
          pe=matrix(0, nrow=traits, ncol=traits)
          pe[1:sum(rowSums(mPe)>0),1:sum(rowSums(mPe)>0)]=mPe[rowSums(mPe)>0, colSums(mPe)>0]
          mPe=pe
        }
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+traits+3))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dPe=as.matrix(diag(mPe))
        dR=as.matrix(diag(mR))
        #----------------------------------------------------------------------------------#
        Hapre=matrix(NA, nrow=nrow(arq), ncol=pre)
        if(covAM==0){
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          #----------------------------------------------------------------------------------#
          if (post==0){post=pre} else {post=post}
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+(arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #----------------------------------------------------------------------------------#
          Hmpre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hmpre[,i]=aux
            colnames(Hmpre)=c(paste('Hmpre', seq(1:pre), sep=''))
          }
          Hm=Hmpre
          #---------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
          #----------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        } else {
          for(i in 1:pre){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Hapre[,i]=aux
            colnames(Hapre)=c(paste('Hapre', seq(1:pre), sep=''))
          }
          #----------------------------------------------------------------------------------#
          if (post==0){post=pre} else {post=post}
          Hapost=matrix(NA, nrow=nrow(arq), ncol=post)
          for(i in 1:post){
            aux=arq[ ,dG[pre+i]]/(arq[ ,dG[pre+i]]+(arq[ ,dR[pre+i]]))
            Hapost[,i]=aux
            colnames(Hapost)=c(paste('Hapost', seq(1:post), sep=''))
          }
          Ha=cbind(Hapre, Hapost)
          #----------------------------------------------------------------------------------#
          Hmpre=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dG[traits+i]]/(arq[ ,dG[i]]+
                                       (arq[ ,dG[traits+i]])+
                                       (abs((arq[ ,mG[i,traits+i]])))+
                                       (arq[ ,dPe[i]])+
                                       (arq[ ,dR[i]]))
            Hmpre[,i]=aux
            colnames(Hmpre)=c(paste('Hmpre', seq(1:pre), sep=''))
          }
          Hm=Hmpre
          #----------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG[1:traits+pre, ])[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
          #----------------------------------------------------------------------------------#
          c2=matrix(NA, nrow=nrow(arq), ncol=pre)
          for(i in 1:pre){
            aux=arq[ ,dPe[i]]/(arq[ ,dG[i]]+
                                 (arq[ ,dG[traits+i]])+
                                 (abs((arq[ ,mG[i,traits+i]])))+
                                 (arq[ ,dPe[i]])+
                                 (arq[ ,dR[i]]))
            c2[,i]=aux
            colnames(c2)=c(paste('c', seq(1:pre), sep=''))
          }
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=c(paste("Vga",seq(1:traits),sep=""), paste("Vgm",seq(1:pre),sep="")))
        vPe=dPe; vPe[vPe==0]=NA; vPe=na.omit(vPe); vPe=data.frame(id=vPe, comp=paste("Vmpe", seq(1:nrow(vPe)), sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vPe,vR,covG)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=data.frame(cbind(arq, CorrG, Ha, Hm, c2))
        arq=as.matrix(arq); mode(arq)="numeric"
        #----------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #----------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
      if(R1==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+traits+3))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dPe=as.matrix(diag(mPe))
        dR=as.matrix(diag(mR))
        #--------------------------------------------------------------------------------------#
        Ha=matrix(NA, nrow=nrow(arq), ncol=traits)
        if(covAM==0){
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          t=matrix(NA, nrow=nrow(arq), ncol=traits)
          for(i in 1:traits){
            aux=(arq[ ,dG[i]]+arq[ ,dPe[i]])/(arq[ ,dG[i]]+
                                                (arq[ ,dG[traits+i]])+
                                                (arq[ ,dPe[i]])+
                                                (arq[ ,dR[i]]))
            t[,i]=aux
            colnames(t)=c(paste('t', seq(1:traits), sep=''))
          }
          #---------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
        } else {
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          t=matrix(NA, nrow=nrow(arq), ncol=traits)
          for(i in 1:traits){
            aux=(arq[ ,dG[i]]+arq[ ,dPe[i]])/(arq[ ,dG[i]]+
                                                (arq[ ,dG[traits+i]])+
                                                (abs((arq[ ,mG[i,traits+i]])))+
                                                (arq[ ,dPe[i]])+
                                                (arq[ ,dR[i]]))
            t[,i]=aux
            colnames(t)=c(paste('t', seq(1:traits), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=paste("Vga",seq(1:traits),sep=""))
        vPe=dPe; vPe[vPe==0]=NA; vPe=na.omit(vPe); vPe=data.frame(id=vPe, comp=paste("Vpe", seq(1:nrow(vPe)), sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        covPe=mPe[upper.tri(mPe)]; covPe[covPe==0]=NA; covPe=c(na.omit(covPe));covPe=sort(covPe)
        matrixPe=mPe; diag(matrixPe)=0; matrixPe=which(matrixPe!=0,arr.ind=T)        
        if (nrow(matrixPe)>1){ matrixPe=matrixPe[order(matrixPe[,1], matrixPe[,2]), ] }
        covPe=data.frame(id=covPe, comp=paste("COV","pe",matrixPe[,1],"pe",matrixPe[,2], sep=""))
        #--------------------------------------------------------------------------------------#
        covR=mR[upper.tri(mR)]; covR[covR==0]=NA; covR=c(na.omit(covR));covR=sort(covR)
        matrixR=mR; diag(matrixR)=0; matrixR=which(matrixR!=0, arr.ind=T)        
        if (nrow(matrixR)>1){ matrixR=matrixR[order(matrixR[,1], matrixR[,2]), ] }
        covR=data.frame(id=covR, comp=paste("COV","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=matrix(NA, nrow=nrow(arq), ncol=((traits^2)-traits)/2)
        subCorE=mR
        colnames(arq)=c(1:ncol(arq))
        subCorE[lower.tri(subCorE)]=lowerTriangle(t(subCorE))
        aux=as.matrix(subCorE)
        for(n in 1:nrow(arq)){
          for(k in 1:ncol(arq)){
            for(i in 1:traits){
              for(j in 1:traits){
                if(subCorE[i,j]==as.numeric(colnames(arq)[k])){
                  aux[i,j]=arq[n, as.numeric(colnames(arq)[k])]
                }
              }
            }
          }
          aux=cov2cor(aux)
          aux=aux[lower.tri(aux)]
          CorrE[n,]=aux
          aux=as.matrix(subCorE)
        }
        CorrE=data.frame(id=seq(1:ncol(CorrE)), t(CorrE))
        CorrE=CorrE[!rowSums(CorrE==0)>=1, ]
        corE=data.frame(id=CorrE[,1], comp=paste("COR","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=merge(corE, CorrE, by=c("id", "id"), sort=FALSE)
        CorrE=t(CorrE)
        colnames(CorrE)=CorrE[2,]
        CorrE=as.matrix(CorrE)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vPe,vR,covG,covPe,covR)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=cbind(arq, CorrG, CorrE, Ha, t)
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #---------------------------------------------------------------------------------------#
      if(R2==TRUE){
        #----------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+traits+3))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dPe=as.matrix(diag(mPe))
        dR=as.matrix(diag(mR))
        #--------------------------------------------------------------------------------------#
        Ha=matrix(NA, nrow=nrow(arq), ncol=traits)
        if(covAM==0){
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          t=matrix(NA, nrow=nrow(arq), ncol=traits)
          for(i in 1:traits){
            aux=(arq[ ,dG[i]]+arq[ ,dPe[i]])/(arq[ ,dG[i]]+
                                                (arq[ ,dG[traits+i]])+
                                                (arq[ ,dPe[i]])+
                                                (arq[ ,dR[i]]))
            t[,i]=aux
            colnames(t)=c(paste('t', seq(1:traits), sep=''))
          }
          #---------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
        } else {
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          t=matrix(NA, nrow=nrow(arq), ncol=traits)
          for(i in 1:traits){
            aux=(arq[ ,dG[i]]+arq[ ,dPe[i]])/(arq[ ,dG[i]]+
                                                (arq[ ,dG[traits+i]])+
                                                (abs((arq[ ,mG[i,traits+i]])))+
                                                (arq[ ,dPe[i]])+
                                                (arq[ ,dR[i]]))
            t[,i]=aux
            colnames(t)=c(paste('t', seq(1:traits), sep=''))
          }
          #----------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=paste("Vga",seq(1:traits),sep=""))
        vPe=dPe; vPe[vPe==0]=NA; vPe=na.omit(vPe); vPe=data.frame(id=vPe, comp=paste("Vpe", seq(1:nrow(vPe)), sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        covR=mR[upper.tri(mR)]; covR[covR==0]=NA; covR=c(na.omit(covR));covR=sort(covR)
        matrixR=mR; diag(matrixR)=0; matrixR=which(matrixR!=0, arr.ind=T)        
        if (nrow(matrixR)>1){ matrixR=matrixR[order(matrixR[,1], matrixR[,2]), ] }
        covR=data.frame(id=covR, comp=paste("COV","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=matrix(NA, nrow=nrow(arq), ncol=((traits^2)-traits)/2)
        subCorE=mR
        colnames(arq)=c(1:ncol(arq))
        subCorE[lower.tri(subCorE)]=lowerTriangle(t(subCorE))
        aux=as.matrix(subCorE)
        for(n in 1:nrow(arq)){
          for(k in 1:ncol(arq)){
            for(i in 1:traits){
              for(j in 1:traits){
                if(subCorE[i,j]==as.numeric(colnames(arq)[k])){
                  aux[i,j]=arq[n, as.numeric(colnames(arq)[k])]
                }
              }
            }
          }
          aux=cov2cor(aux)
          aux=aux[lower.tri(aux)]
          CorrE[n,]=aux
          aux=as.matrix(subCorE)
        }
        CorrE=data.frame(id=seq(1:ncol(CorrE)), t(CorrE))
        CorrE=CorrE[!rowSums(CorrE==0)>=1, ]
        corE=data.frame(id=CorrE[,1], comp=paste("COR","e",matrixR[,1],"e",matrixR[,2], sep=""))
        CorrE=merge(corE, CorrE, by=c("id", "id"), sort=FALSE)
        CorrE=t(CorrE)
        colnames(CorrE)=CorrE[2,]
        CorrE=as.matrix(CorrE)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vPe,vR,covG,covR)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=data.frame(cbind(arq, CorrG, CorrE, Ha, t))
        arq=as.matrix(arq); mode(arq)="numeric"
        #----------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #----------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #--------------------------------------------------------------------------------------#
      }
      #----------------------------------------------------------------------------------------#
      if(R3==TRUE){
        #--------------------------------------------------------------------------------------#
        cat("Estimating genetic and environmental parameters...\n")
        mG=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[2:(n1matrix+1)]), " "))))), nrow=n1matrix, byrow=T)
        mPe=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[c((n1matrix+3):(n1matrix+traits+2))]), " "))))), nrow=traits, byrow=T) 
        mR=matrix(na.omit(as.numeric(as.matrix(unlist(strsplit(as.character(Matrix[-c(1:(n1matrix+traits+3))]), " "))))), nrow=traits, byrow=T)
        dG=as.matrix(diag(mG))
        dPe=as.matrix(diag(mPe))
        dR=as.matrix(diag(mR))
        #--------------------------------------------------------------------------------------#
        Ha=matrix(NA, nrow=nrow(arq), ncol=traits)
        if(covAM==0){
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          t=matrix(NA, nrow=nrow(arq), ncol=traits)
          for(i in 1:traits){
            aux=(arq[ ,dG[i]]+arq[ ,dPe[i]])/(arq[ ,dG[i]]+
                                                (arq[ ,dG[traits+i]])+
                                                (arq[ ,dPe[i]])+
                                                (arq[ ,dR[i]]))
            t[,i]=aux
            colnames(t)=c(paste('t', seq(1:traits), sep=''))
          }
          #---------------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          CorrG = matrix(NA, nrow=nrow(arq), ncol=2*((effectG^2-effectG)/2))
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:dG[traits]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,1:length(aux)]=aux
            aux=vec2sm(arq[i,dG[traits+1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,(length(aux)+1):ncol(CorrG)]=aux
          }
          options(warn=0)
        } else {
          for(i in 1:traits){
            aux=arq[ ,dG[i]]/(arq[ ,dG[i]]+
                                (arq[ ,dG[traits+i]])+
                                (abs((arq[ ,mG[i,traits+i]])))+
                                (arq[ ,dPe[i]])+
                                (arq[ ,dR[i]]))
            Ha[,i]=aux
            colnames(Ha)=c(paste('ha', seq(1:traits), sep=''))
          }
          #--------------------------------------------------------------------------------------#
          t=matrix(NA, nrow=nrow(arq), ncol=traits)
          for(i in 1:traits){
            aux=(arq[ ,dG[i]]+arq[ ,dPe[i]])/(arq[ ,dG[i]]+
                                                (arq[ ,dG[traits+i]])+
                                                (abs((arq[ ,mG[i,traits+i]])))+
                                                (arq[ ,dPe[i]])+
                                                (arq[ ,dR[i]]))
            t[,i]=aux
            colnames(t)=c(paste('t', seq(1:traits), sep=''))
          }
          #----------------------------------------------------------------------------------#
          effectG=nnzero(dG); CorrG=matrix(NA, nrow=nrow(arq), ncol=((effectG^2-effectG)/2)) 
          options(warn=-1)
          for(i in 1:nrow(arq)){
            aux=vec2sm(arq[i,dG[1]:rev(dG)[1]])
            aux=cov2cor(aux)
            aux=aux[lower.tri(aux)]
            CorrG[i,]=aux
          }
          options(warn=0)
        }
        #----------------------------------------------------------------------------------#
        vG=dG; vG[vG==0]=NA; vG=na.omit(vG); vG=data.frame(id=vG, comp=paste("Vga",seq(1:traits),sep=""))
        vPe=dPe; vPe[vPe==0]=NA; vPe=na.omit(vPe); vPe=data.frame(id=vPe, comp=paste("Vpe", seq(1:nrow(vPe)), sep=""))
        vR=dR; vR[vR==0]=NA; vR=na.omit(vR); vR=data.frame(id=vR, comp=paste("Ve", seq(1:nrow(vR)), sep=""))
        #--------------------------------------------------------------------------------------#
        covG=mG[upper.tri(mG)]; covG[covG==0]=NA; covG=c(na.omit(covG));covG=sort(covG)
        matrixG=mG; diag(matrixG)=0; matrixG=which(matrixG!=0,arr.ind=T)        
        if (nrow(matrixG)>1){ matrixG=matrixG[order(matrixG[,1], matrixG[,2]), ] }
        covG=data.frame(id=covG, comp=paste("COV","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=data.frame(id=seq(1:ncol(CorrG)), t(CorrG))
        corG=data.frame(id=CorrG[,1], comp=paste("COR","g",matrixG[,1],"g",matrixG[,2], sep=""))
        CorrG=merge(corG, CorrG, by=c("id", "id"), sort=FALSE)
        CorrG=t(CorrG)
        colnames(CorrG)=CorrG[2,]
        CorrG=as.matrix(CorrG)[-c(1,2),]
        #--------------------------------------------------------------------------------------#
        componentes=rbind(vG,vPe,vR,covG)
        arq=data.frame(id=seq(1:ncol(arq)), t(arq))
        arq=merge(componentes, arq, by=c("id", "id"), sort=FALSE)
        arq=t(arq)
        colnames(arq)=arq[2,]
        arq=as.matrix(arq)[-c(1,2),]
        mode(arq)="numeric"
        rownames(arq)=NULL
        #--------------------------------------------------------------------------------------#
        arq=data.frame(cbind(arq, CorrG, Ha, t))
        arq=as.matrix(arq); mode(arq)="numeric"
        #--------------------------------------------------------------------------------------#
        cat("\nSummarizing postgibbs results...\n")
        x=data.frame(arq)
        Summary_all=summary(as.mcmc(x))
        summary=data.frame(Summary_all[[1]], Summary_all[[2]])
        summary=data.frame(Mean=summary[,1],
                           data.frame(Mode=apply(x, 2, estimate_mode)),
                           Median=summary[,7],
                           SD=summary[,2],
                           HPD_2.5=summary[,5],
                           HPD_97.5=summary[,9],
                           Naive_SE=summary[,3],
                           Time_series_SE=summary[,4])
        stat=summary
        #--------------------------------------------------------------------------------------#
        cat("\n")
        cat("Creating PDF to save the plots...\n")
        pdf("Postgibbs_Plot.pdf")
        PDF=ncol(arq)
        par(mfrow = c(PDF, 2))  # 3 rows and 2 columns
        ## Make a plot of all the parameters in the dataset and show some Kernel Density
        ## estimates for the marginal posteriors
        plot(as.mcmc(arq), trace=TRUE, density=TRUE, smooth=FALSE, auto.layout=TRUE, lwd=line)
        cat("\nKernel Density estimates for the marginal posteriors: done\n")
        dev.off()
        #----------------------------------------------------------------------------------------#
      }
    }
    #------------------------------------------------------------------------------------------#
    cat("\n")
    # Check For Perl first
    res <- Sys.which("perl")
    if (res!=""){
      message("Perl found.\n")
      cat("As you have perl installed, the summary of statistical results will be saved in Excel format.\n");
      x=data.frame(Parameter=rownames(stat), stat);
      rownames(x)=NULL;
      WriteXLS(x="x",
               ExcelFileName="PostStat.xls",
               AdjWidth=TRUE,
               BoldHeaderRow=TRUE,
               FreezeCol=1,
               FreezeRow=1)
    }else{
      message("Perl was not found on your system. Either check $PATH if installed or please install Perl.\n")
      cat("As you don't have perl installed, the summary of statistical results will be saved in CSV format.\n");
      write.csv(stat, "PostStat.csv", quote=F);
    }
  }
  ##----------------------------------------------------------------------------------------##
  ##----------------------------------------------------------------------------------------##
  ##                                                                                        ##
  ##                                                                                        ##
  ##                                                                                        ##
  ##                      Inductive Causation (IC) algorithm and graph                      ##
  ##                                                                                        ##
  ##                                                                                        ##
  ##                                                                                        ##
  ##----------------------------------------------------------------------------------------##
  ##----------------------------------------------------------------------------------------##
  if(ICgraph==TRUE){
    hpdContent=HPD
    cat("\nReading residual postgibbs samples to use in the Inductive Causation (IC) Algorithm...\n")
    n_R=grep("R matrix", Matrix)
    nR=as.matrix(unlist(strsplit(as.character(Matrix[-c(1:n_R)]), " ")))
    nR[nR==""] <- NA
    nR=na.omit(nR)
    nR=nR[1]
    mode(nR)="numeric"
    r=r[,nR:ncol(r)]
    ##----------------------------------------------------------------------------------------##
    ## Final transposition corrected
    ## Getting partially, oriented graph, connected pairs and unshielded colliders
    ## Complete IC Algorithm 
    ## Adapted to work with variable number of traits
    ## Input: sample of covariance matrix (vech) in 'r'; and optionally: burnin, thinning, 
    ## HPDContent and names of traits 
    ##----------------------------------------------------------------------------------------##
    options(warn=-1)
    n=ncol(vec2sm(as.numeric(r[1,])))
    cat("\nRunning IC algorithm to create Causal Graphic.\n")
    sets=t(combn(c(1:n),2))              
    nPairs=nrow(sets)
    csets=matrix(0,nrow(sets),(n-2))      
    for(t in 1:nrow(sets)){
      c=matrix(TRUE,1,n)                            
      c[sets[t,]]=FALSE
      csets[t,]= c(1:n)[c]                          
    }
    connectedPairs=matrix(FALSE,nrow(sets),1)
    indSets=list()
    cPair=1
    #--------------------------------------------------------------------------------------------#
    # First step
    #--------------------------------------------------------------------------------------------#
    for(p in 1:nrow(sets)){
      time<-proc.time()[3]
      ind=0                                                                       
      setOfSets=list(c())                                       
      for(u in 1:(n-2)){
        setOfSets3=data.frame()
        comb=combn(c(1:(n-2)),u)
        for(v in 1:ncol(comb)){
          setOfSets2=list(csets[p,c(comb[,v])])
          setOfSets3=c(setOfSets3,setOfSets2)
        }
        setOfSets=c(setOfSets,setOfSets3)
        rm(setOfSets3)
      }
      nPCors=length(setOfSets)
      sPCor=1
      pcor = matrix(0,nrow(r),nPCors)                      
      for(j in 1:nPCors){
        set = c(sets[p,],setOfSets[[j]])
        time2<-proc.time()[3]
        for(i in 1:nrow(r)){
          R=vec2sm(r[i,])
          R2 = R[set,set]
          invR2 <- solve(R2)
          pcor[i,j]= -(cov2cor(invR2)[1,2])
        }
        tmp2 <- proc.time()[3]
        cat(paste(c("PCor: ", "of ", "time: "), c(sPCor, nPCors, round(tmp2 - time2, 4))))
        cat("\n")
        sPCor<-sPCor+1
      }
      for(k in 1:ncol(pcor)){
        pcorMCMC<-as.mcmc(pcor[,k])                          
        hpdNum = as.numeric(HPDinterval(pcorMCMC,prob=hpdContent))
        if(hpdNum[1]<0 & hpdNum[2]>0){
          ind=ind+1
          indSets[[length(indSets)+1]]=list(sets[p,],setOfSets[[k]])
        }
      }
      if(ind==0){
        connectedPairs[p]=TRUE
      }
      tmp <- proc.time()[3]
      cat(paste(c("Pair: ", "of ", "time:"), c(cPair, nPairs, round(tmp - time, 4))))
      cat("\n")
      cat(paste("-----------------------------------------------"))
      cat("\n")
      cPair<-cPair+1
    }
    dConnect=sets[connectedPairs,]
    #--------------------------------------------------------------------------------------------#
    #Second Step
    #--------------------------------------------------------------------------------------------#
    unshColl=list()
    if(sum(as.numeric(connectedPairs))>1){
      for(i in 1:(nrow(dConnect)-1)){
        for(j in (i+1):nrow(dConnect)){
          if(TRUE%in%(dConnect[i,]%in%dConnect[j,])){
            colParCand=c(dConnect[i,],dConnect[j,])[c((which(!c(dConnect[i,],dConnect[j,])%in%dConnect[i,])),(which(!c(dConnect[i,],dConnect[j,])%in%dConnect[j,])))]
            colCand=unique(c(dConnect[i,],dConnect[j,])[!c(dConnect[i,],dConnect[j,])%in%colParCand])
            disc=TRUE
            for(k in 1:nrow(dConnect)){
              if((sum((colParCand-dConnect[k,])^2)==0)|(sum((rev(colParCand)-dConnect[k,])^2)==0))
                disc=FALSE
            }
            if(disc==TRUE){
              isUC=TRUE
              for(l in 1:length(indSets)){
                if((sum((colParCand-indSets[[l]][[1]])^2)==0)|(sum((rev(colParCand)-indSets[[l]][[1]])^2)==0)){
                  if(colCand %in% indSets[[l]][[2]]){
                    isUC=FALSE
                  }
                }
              }
              if(isUC==TRUE)
              {
                unshColl[[length(unshColl)+1]]=list(colCand,colParCand)
              }
            }
          }
        }
      }
    }
    if(length(nodeNames)==0){
      nodeNames=as.character(1:n)
    }
    POG<-new("graphNEL", nodes = nodeNames, edgemode = "undirected")
    # Changes in Dec 7 2016 - Fernando Brito
    if(is.vector(dConnect)){
      POG<-addEdge(from=nodeNames[dConnect[1]], to=nodeNames[dConnect[2]], POG, 1)
    }else{
      POG<-addEdge(from=nodeNames[dConnect[,1]], to=nodeNames[dConnect[,2]], POG, 1)
    }
    if(length(unshColl)>0){
      temp1<-c()
      temp2<-c()
      for(i in 1:length(unshColl)){
        temp1<-rbind(temp1,unshColl[[i]][[1]])
        temp2<-cbind(temp2,t(unshColl[[i]][[2]]))
      }  
      dirE<-cbind(t(temp2),rep(temp1,each=2))
    }
    #--------------------------------------------------------------------------------------------#
    # Step 3
    #--------------------------------------------------------------------------------------------#
    pdag<-as(POG, "matrix")
    p<-dim(pdag)[1]
    if("dirE"%in%ls()){
      for(i in 1:nrow(dirE)){
        pdag[dirE[i,1],dirE[i,2]]<-1
        pdag[dirE[i,2],dirE[i,1]]<-0
      }
      old_pdag <- matrix(0, p, p)
      while (!all(old_pdag == pdag)) {
        old_pdag <- pdag
        ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pdag[b, ] == 1 & pdag[, b] == 1) & (pdag[a, ] == 0 & pdag[, a] == 0))
          if (length(indC) > 0) {
            pdag[b, indC] <- 1
            pdag[indC, b] <- 0
          }
        }
        # Madison, Feb. 23, 2017 - Fernando Brito
        # Changed from: ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE)
        #           to: ind <- which((pdag == 0 & t(pdag) == 1), arr.ind = TRUE)
        ind <- which((pdag == 0 & t(pdag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pdag[a, ] == 1 & pdag[, a] == 0) & (pdag[, b] == 1 & pdag[b, ] == 0))
          if (length(indC) > 0) {
            pdag[a, b] <- 1
            pdag[b, a] <- 0
          }
        }
        ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pdag[a, ] == 1 & pdag[, a] == 1) & (pdag[, b] == 1 & pdag[b, ] == 0))
          if (length(indC) >= 2) {
            g2 <- pdag[indC, indC]
            if (length(g2) <= 1) {
              g2 <- 0
            }
            else {
              diag(g2) <- rep(1, length(indC))
            }
            if (any(g2 == 0)) {
              pdag[a, b] <- 1
              pdag[b, a] <- 0
            }
          }
        }
      }
    }
    # graph plot
    POG2 <- as(pdag, "graphNEL")
    dirE2 <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE)
    if(nrow(dirE2)>0){
      dirECh<-as.character(vector(length=nrow(dirE2)))
      for(i in 1:nrow(dirE2)){
        dirECh[i]<-paste(nodeNames[dirE2[i,1]],"~",nodeNames[dirE2[i,2]],sep="")
      }
      edges <- buildEdgeList(POG2)
      nodes <- buildNodeList(POG2)
      for(i in 1:length(edges)){
        edges[[i]]@attrs$arrowhead <-"none"
        edges[[i]]@attrs$dir<-"none"
      }
      if(length(unshColl)>0){
        for(i in 1:length(dirECh)){
          edges[[dirECh[i]]]@attrs$arrowhead <-"open"
          edges[[dirECh[i]]]@attrs$dir<-"forward"
        }
      }
      vv <- agopen(name="name", nodes=nodes, edges=edges, edgeMode=edgemode(POG2))
      png(paste("ICgraph_",HPD*100,".png",sep=""))
      plot(vv, attrs=list(node=list(fillcolor="white"), edge=list(splines="line", arrowsize=1, arrowType="vee")))
      dev.off()
      plot(vv, attrs=list(node=list(fillcolor="white"), edge=list(splines="line", arrowsize=1, arrowType="vee")))
    } else {
      png(paste("ICgraph_",HPD*100,".png",sep=""))
      plot(POG2, attrs=list(node=list(fillcolor="white"), edge=list(splines="line", arrowsize=1, arrowType="vee")))
      dev.off()
      plot(POG2, attrs=list(node=list(fillcolor="white"), edge=list(splines="line", arrowsize=1, arrowType="vee")))
    }
  }
  gib.samples <- list.files(pattern = "gib.samples")
  if(!identical(gib.samples,character(0))){
    RUNfiles <- c("renf90.par","postind","postgibbs_samples")
    for(remove in 1:length(RUNfiles)) file.remove(RUNfiles[remove])
  }
}
