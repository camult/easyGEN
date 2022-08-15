#' @title A  function interfacing gibbs2f90
#' 
#' @description The function uses the model formula language in R to describe the model and from this generates the files needed to do an analysis using gibbs2f90 software!
#' 
#' @param formula A two-sided linear formula object describing the fixed effects part of the model, 
#'        with the responses on the left of a ~ operator (separated by | if more then one response)
#'        and the terms, separated by + operators, on the right.
#' @param phen Name of phenotype file.
#' @param ped Name of pedigree file
#' @param geno Name of genotype file.
#' @param map Name of map file.
#' @param idName Identification name of animal's column.
#' @param diffVAR List of the trait and their unique effect, if any.
#' @param covariate Name of covariate(s) if there is(are) any.
#' @param OPTeff OPTIONAL effect, i.e., mat, pe, mpe. It is NULL if it is not given.
#' @param OPTlist OPTIONAL effect. It is NULL if it is not given.
#' @param Gcov Genetic (co)variance. It is NULL if it is not given.
#' @param Rcov Residual (co)variance. It is NULL if it is not given.
#' @param nMaternal Number of traits on maternal effect. By default nMaternal=0.
#' @param execute Whether to run the gibbsf90. By default execute=TRUE.
#' @param missing a integer specifying the missing value (default 0).
#' @param nIter,burnIn,thin	(integer) the number of iterations, burn-in and thinning.
#' @param intern	a logical (not NA) which indicates whether to capture the output of the command as an R character vector.
#' @param weight name of the column (numeric) with a vector of weights, may be NULL.  If weights is not NULL, the residual variance of each data-point is set to be proportional to the inverse of the squared-weight.
#' @param PED_DEPTH a integer specifying the depth of pedigree search. The default is 3, byt all pedigrees are loaded if it is set to 0. 
#' @param covAM type 0 if covariance between additive and maternal genetic effects must be fixed in zero, 1 otherwise. By default covAM=1.
#' @param covR type 0 if residual covariance must be fixed in zero, 1 otherwise. By default covAM=1.
#' @param useF logical value indicating whether inbreeding coefficient of this animal should be computed. Default is FALSE.
#' 
#' @return gibbsf90 results. 
#'         
#' @references Misztal I, Tsuruta S.,  Lourenco D., Aguilar I., Legarra A., 
#'             Vitezica Z (2014). Manual  for BLUPF90   family	of programs. 
#'             Available at: <http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all1.pdf>.
#'
#' @export gibbsf90
#' 
#' @import gdata Matrix
##------------------------------------------------------------------------------------------##
gibbsf90 <- function(formula, phen, ped=NULL, geno=NULL, map=NULL, idName, 
                     diffVAR=NULL, nMaternal=NULL, weight=NULL, 
                     nIter=1500, burnIn=500, thin = 5,  missing=0,
                     covariate=0, OPTeff=NULL, OPTlist=NULL, intern=TRUE,
                     Gcov=NULL, Rcov=NULL, execute=TRUE, PED_DEPTH=0, 
                     covAM=1, covR=1, useF=FALSE){
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
  cat("Analyzing the given model...\n")
  ##----------------------------------------------------------------------------------------##
  model <- formula
  trait  <- unlist(strsplit(as.character(model[[2]])," "))[!unlist(strsplit(as.character(model[[2]])," "))%in%c("|")]
  nTr <- length(trait)
  effect <- unlist(strsplit(as.character(model[[3]])," "))[!unlist(strsplit(as.character(model[[3]])," "))%in%c("+")]
  nEff <- length(effect)
  if("1"%in%effect){
    file <- phen[, c(idName, trait)]
    colnames(file) <- c(idName, trait)
  }else{
    file <- phen[, c(idName, trait, effect)]
    colnames(file) <- c(idName, trait, effect)
  }
  file <- file[, !duplicated(colnames(file))]
  if(is.null(Rcov)){
    if(nTr>1){
      Rcov <- round(var(file[, trait], use="complete.obs")*0.6, 5)
    }else{
      Rcov <- round(var(file[, trait], use="complete.obs")*0.6, 5)
    }
  }
  if(is.null(Gcov)){
    if(nTr>1){
      if(length(grep("mat", OPTeff))==0){
        Gcov <- round(var(file[, trait], use="complete.obs")*0.4, 5)
      }else{
        if(nMaternal==nTr){
          Gcov <- round(var(file[, rep(trait,2)], use="complete.obs")*0.4, 5)
        }
        if(nMaternal<nTr){
          Gcov <- round(var(file[, rep(trait,2)], use="complete.obs")*0.4, 5)
          Gcov[(nrow(Gcov)-(nTr-nMaternal)+1):nrow(Gcov), ] <- 0
          Gcov[, (nrow(Gcov)-(nTr-nMaternal)+1):nrow(Gcov)] <- 0
        }
      }
    }else{
      if(length(grep("mat", OPTeff))==0){
        Gcov <- round(var(file[, trait], use="complete.obs")*0.4, 5)
      }else{
        if(nMaternal==nTr){
          Gcov <- round(var(file[, rep(trait,2)], use="complete.obs")*0.4, 5)
        }
        if(nMaternal<nTr){
          Gcov <- round(var(file[, rep(trait,2)], use="complete.obs")*0.4, 5)
          Gcov[(nrow(Gcov)-(nTr-nMaternal)+1):nrow(Gcov), ] <- 0
          Gcov[, (nrow(Gcov)-(nTr-nMaternal)+1):nrow(Gcov)] <- 0
        }
      }
    }
    if(covAM == 0){
      Gcov[1:nTr, (nTr+1):nrow(Gcov)] <- 0
      Gcov[(nTr+1):nrow(Gcov), 1:nTr] <- 0
    }
  }
  if(is.null(ped)){
    pedigree <- data.frame(animal=file[, idName], sire=0, dam=0)
    cat("   writing fake pedigree file...\n")
    write.fwf(pedigree,  "pedigree.dat",  rownames=F, colnames=F, quote=F, justify="right")
  }else{
    pedigree <- ped
    cat("   writing pedigree file...\n")
    write.fwf(pedigree,  "pedigree.dat",  rownames=F, colnames=F, quote=F, justify="right")
  }
  if(is.null(covariate)){
    covariate=NULL
  }else{
    covariate=covariate
  }
  covar <- unique(unlist(sapply(covariate, grep, effect)))
  if(length(covar)==0) covar <- NULL
  dimCovar <- length(covar)
  listCovariate <- list()
  if(dimCovar>0){
    for(i in 1:dimCovar) listCovariate[i] <- covariate[i]
    listCovariate <- c(do.call(cbind, listCovariate))
    listCovariate <- effect[c(covar)]
  }
  if(is.null(covar)){
    if("1"%in%effect){
      ClassEffects <- NULL
    }
  }
  if(is.null(covar)){
    if(!"1"%in%effect){
      ClassEffects <- effect
    }
  }else{
    ClassEffects <- effect[-c(covar)]
  }
  dimClassEffect <- length(ClassEffects)
  listClassEffects <- list()
  if(dimClassEffect>0){
    for(i in 1:dimClassEffect) listClassEffects[i] <- ClassEffects[i]
    listClassEffects <- c(do.call(cbind, listClassEffects))
  }
  if(!is.null(weight)){
    keep_phen <- c(idName, trait, if(dimClassEffect>0){listClassEffects}, if(dimCovar>0){listCovariate}, weight)
    wt <- grep(weight, keep_phen)
    keep_phen <- keep_phen[!duplicated(keep_phen)]
    file <- phen[, keep_phen]
  }else{
    keep_phen <- c(idName, trait, if(dimClassEffect>0){listClassEffects}, if(dimCovar>0){listCovariate})
    keep_phen <- keep_phen[!duplicated(keep_phen)]
  }
  #phen <- data.frame(phen[, colnames(phen)%in%keep_phen])
  #phen <- phen[, match(keep_phen, colnames(phen))]
  file <- file[, keep_phen]
  #------------------------------------------------------------#
  # Added on 04 July, 2018
  #------------------------------------------------------------#
  file <- file[rowSums(!is.na(data.frame(file[,trait]))) > 0, ]
  #------------------------------------------------------------#
  file[is.na(file)] <- missing
  PosTrait <- grep(paste(paste0("^",trait,"$"),collapse="|"), colnames(file))
  cat("   writing phenotype file...\n")
  write.fwf(file, "phen.dat", rownames=F, colnames=F, quote=F, justify="right")
  # (SNP x ID)
  if(!is.null(geno)){
    cat("   writing genotype file...\n")
    if(ncol(geno)>nrow(geno)) geno <- t(geno)
    geno <- geno[,colnames(geno)%in%pedigree[,1]]
    nID <- ncol(geno)
    nSNP <- nrow(geno)
    IDnames <- colnames(geno)
    wGeno <- c(max(nchar(IDnames)), nSNP)
    geno <- t(geno)
    for(gi in 1:nID){
      aux <- data.frame(IDnames[gi], paste(geno[gi,],sep="",collapse=""))
      if(gi==1){
        write.fwf(aux, "geno.dat", rownames=F, colnames=F, quote=F, justify="right", width=wGeno)
      }else{
        write.fwf(aux, "geno.dat", rownames=F, colnames=F, quote=F, justify="right", width=wGeno, append=T)
      }
    }
  }
  if(!is.null(map)){
    # map file (marker, chr, bp)
    map  <- map
    cat("   writting map file...\n")
    write.fwf(map,  "map.dat",  rownames=F, colnames=F, quote=F, justify="right")
  }
  if(dimClassEffect>0){
    PosClassEffect <- c(array(0, dimClassEffect))
    for(i in 1:dimClassEffect){
      PosClassEffect[i] <- c(grep(listClassEffects[i], colnames(file)))
    }
  }
  if(dimCovar>0){
    PosCovariate <- unique(unlist(sapply(listCovariate, grep, colnames(file))))
  }
  EffMat <- matrix(0, nrow=nEff, ncol=nTr)
  rownames(EffMat) <- effect
  colnames(EffMat) <- trait
  EffVec <- paste0(diffVAR, names(diffVAR))
  if(!is.null(diffVAR)){
    Eff_Names <- as.character(diffVAR)
    Eff_All <- ifelse(effect%in%Eff_Names, effect, "NoEff")
    for(i in 1:nrow(EffMat)){
      for(j in 1:ncol(EffMat)){
        EffMat[i,j] <- paste0(Eff_All[i], colnames(EffMat)[j])
      }
    }
    EffMat[grep(paste0(EffVec, collapse="|"), EffMat)] <- 1
    EffMat <- ifelse(EffMat=="1",1,"")
    mode(EffMat) <- "numeric"
    for(i in 1:nrow(EffMat)){
      for(j in 1:ncol(EffMat)){
        EffMat[i,j] <- ifelse(EffMat[i,j]==1, grep(Eff_All[i], keep_phen))
      }
    }
    commonVAR <- effect[-grep(paste(diffVAR,collapse="|"), effect)]
  }else{
    if(!"1"%in%effect){
      commonVAR <- effect
    }else{
      commonVAR <- NULL
    }
  }
  if(length(commonVAR)>0){
    for(i in 1:nrow(EffMat)){
      for(j in 1:length(commonVAR)){
        if(rownames(EffMat)[i]==commonVAR[j]){
          EffMat[i,] <- rep(grep(paste0("^",commonVAR[j],"$"), keep_phen), nTr)
        }
      }
    }
  }
  EffMat <- data.frame(EffMat)
  EffMat$Effec <- ""
  for(i in 1:nrow(EffMat)){
    if(rownames(EffMat)[i]%in%listClassEffects){
      EffMat[i, "Effec"] <- "cross alpha"
    }else if(rownames(EffMat)[i]%in%listCovariate){
      EffMat[i, "Effec"] <- "covariable"
    }
  }
  EffMat[is.na(EffMat)] <- 0
  catMAT <- EffMat[rownames(EffMat)%in%listClassEffects,]
  covMAT <- EffMat[rownames(EffMat)%in%listCovariate,]
  if(nTr>1){
    Rcov <- round(as.matrix(Rcov) + (diag(diag(as.matrix(Rcov)))), 5)
    Gcov <- round(as.matrix(Gcov) + (diag(diag(as.matrix(Gcov)))), 5)
  }
  if(covR==0){
    Rcov <- diag(diag(Rcov))
  }
  
  #---------------------------------------------------------------------------------------#
  # Creating parameter file based on given model
  #---------------------------------------------------------------------------------------#
  cat("\nCreating parameter file...\n")
  sink("renum.par")
  cat("# layout:",colnames(file),"\n")
  cat("DATAFILE\n")
  cat("phen.dat\n")
  cat("TRAITS\n")
  cat(PosTrait,"\n")
  cat("FIELDS_PASSED TO OUTPUT\n\n")
  cat("WEIGHT(S)\n")
  cat(ifelse(length(weight)==0, "", wt),"\n")
  cat("RESIDUAL_VARIANCE\n")
  sink()
  write.table(Rcov,file="renum.par", col.names=F, row.names=F, quote=F, sep=" ", append=T)
  sink("renum.par", append=T)
  if(dimClassEffect>0){
    for(i in 1:dimClassEffect){
      cat("EFFECT\n",as.character(catMAT[i,]),"\n")
    }
  }
  if(dimCovar>0){
    for(i in 1:dimCovar){
      cat("EFFECT\n",as.character(covMAT[i,]),"\n")
    }
  }
  cat("EFFECT\n")
  cat(rep(grep(idName,keep_phen),nTr),"cross alpha\n")
  cat("RANDOM\n")
  cat("animal\n")
  if(!is.null(OPTeff)){
    cat("OPTIONAL\n")
    cat(as.character(OPTeff), "\n")
  }
  cat("FILE\n")
  cat("pedigree.dat\n")
  if(!is.null(geno)){
    cat("SNP_FILE\n")
    cat("geno.dat\n")
  }
  cat("PED_DEPTH\n")
  cat(PED_DEPTH,"\n")
  if(useF){
    cat("INBREEDING\n")
    cat("pedigree\n")
  }
  cat("(CO)VARIANCES\n")
  sink()
  write.table(Gcov,file="renum.par", col.names=F, row.names=F, quote=F, sep=" ", append=T)
  if(!is.null(OPTlist)){
    write.table(do.call(rbind, OPTlist),file="renum.par",col.names=F,row.names=F,quote=F,sep=" ",append=T)
  }
  if(missing!=0){
    cat("OPTION missing", missing,"\n")
  }
  #---------------------------------------------------------------------------------------#
  # Creating parameters files needed to run the BLUPf90 proframs on Windows
  #---------------------------------------------------------------------------------------#
  if(isTRUE(execute)){
    cat("\nPerforming quantitative genetic analysis...\n")
    write.table("renum.par",  file="renum",  row.names=F, col.names=F, quote=F)
    write.table("renumf90.exe < renum",   file="renum.bat",     row.names=F, col.names=F, quote=F)
    sink("renf90")
    cat("renf90.par\n")
    cat(nIter, burnIn,"\n")
    cat(thin,"\n")
    sink()
    write.table("gibbs2f90.exe < renf90",   file="gibbs2f90.bat",   row.names=F, col.names=F, quote=F)
    sink("postgibbs")
    cat("renf90.par\n")
    cat(0,"\n")
    cat(thin,"\n")
    cat(0,"\n")
    sink()
    write.table("postgibbsf90.exe < postgibbs",   file="postgibbsf90.bat",   row.names=F, col.names=F, quote=F)
    #---------------------------------------------------------------------------------------#
    # Running renum
    #---------------------------------------------------------------------------------------#
    cat("   running renumf90 program...\n")
    system("renum.bat", intern=intern)
    #---------------------------------------------------------------------------------------#
    # Running gibbs2f90
    #---------------------------------------------------------------------------------------#
    cat("   running gibbs2f90 program...\n")
    system("gibbs2f90.bat", intern=intern)
    #---------------------------------------------------------------------------------------#
    # Running postgibbsf90
    #---------------------------------------------------------------------------------------#
    cat("   running postgibbsf90 program...\n\n")
    system("postgibbsf90.bat", intern=intern)
    #---------------------------------------------------------------------------------------#
    # Removing files that is no longer needed anymore
    #---------------------------------------------------------------------------------------#
    #cat("Cleaning...\n\n")
    #pattern <- c("rem","ren","phen","ped","map","geno","fsp","gibbs2f90","postgibbsf90")
    #files <- list.files()
    #noRemove <- c("renf90.par","remlf90.log","remlf90.exe","renumf90.exe","Sol_Acc.txt","gibbs2f90.exe","postgibbsf90.exe")
    #Remove <- files[grep(paste(pattern,collapse="|"), files)]
    #Remove <- Remove[-c(grep(paste(noRemove,collapse="|"), Remove))]
    #invisible(file.remove(Remove))
  }
}
