#' @title R  functions for interfacing remlf90+ software
#' 
#' @description The function uses the model formula language in R to describe the model and from this generates the files needed to do an analysis using remlf90 software!
#' 
#' @param formula A two-sided linear formula object describing the fixed effects part of the model, 
#'        with the responses on the left of a ~ operator (separated by "|" if more then one response)
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
#' @param missing a integer specifying the missing value (default 0).
#' @param nMaternal Number of traits on maternal effect. By default nMaternal=0.
#' @param execute Whether to run the remlf90. By default execute=TRUE.
#' @param PED_DEPTH a integer specifying the depth of pedigree search. The default is 3, byt all pedigrees are loaded if it is set to 0. 
#' @param weight name of the column (numeric) with a vector of weights, may be NULL.  If weights is not NULL, the residual variance of each data-point is set to be proportional to the inverse of the squared-weight.
#' @param Inb Whether to run the inbupgf90 to compute the coefficient of inbreeding. By default Inb=FALSE.
#' @param covAM type 0 if covariance between additive and maternal genetic effects must be fixed in zero, 1 otherwise. By default covAM=1.
#' @param covR type 0 if covariance between additive and maternal genetic effects must be fixed in zero, 1 otherwise. By default covAM=1.
#' @param useF logical value indicating whether inbreeding coefficient of this animal should be computed. Default is FALSE.
#' 
#' @return BLUPf90 results. 
#'         
#' @references Misztal I, Tsuruta S.,  Lourenco D., Aguilar I., Legarra A., 
#'             Vitezica Z (2014). Manual  for BLUPF90   family	of programs. 
#'             Available at: <http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all1.pdf>.
#'
#' @export remlf90
#' 
#' @import gdata Matrix
##------------------------------------------------------------------------------------------##
remlf90 <- function(formula, phen, ped=NULL, geno=NULL, map=NULL, idName, 
                    diffVAR=NULL, nMaternal=NULL, weight=NULL, Inb=FALSE,
                    covariate=0, OPTeff=NULL, OPTlist=NULL, missing=0,
                    Gcov=NULL, Rcov=NULL, execute=TRUE, PED_DEPTH=0, 
                    covAM=1, covR=1, useF=FALSE){
  ##----------------------------------------------------------------------------------------##
  cat("\014")
  cat("\n")
  centerText <- function() {
    width <- getOption("width")
    A <- ("                                                       ._____.    \n")
    B <- ("    _//////_            Fernando Brito Lopes           _|_____|_   \n")
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
  }else{
    Rcov <- round(as.matrix(Rcov) + (as.matrix(Rcov)), 5)
    if(covR==0){
      Rcov <- diag(diag(Rcov))
    }
    Gcov <- round(as.matrix(Gcov) + (as.matrix(Gcov)), 5)
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
    write.table("renf90.par", file="renf90", row.names=F, col.names=F, quote=F)
    write.table("renumf90.exe < renum",   file="renum.bat",     row.names=F, col.names=F, quote=F)
    write.table("remlf90.exe < renf90",   file="remlf90.bat",   row.names=F, col.names=F, quote=F)
    #---------------------------------------------------------------------------------------#
    # Running renum
    #---------------------------------------------------------------------------------------#
    cat("   running renumf90 program...\n")
    renum <- data.frame(renum=capture.output(system("renum.bat", intern=TRUE)))
    #---------------------------------------------------------------------------------------#
    # Running remlf90
    #---------------------------------------------------------------------------------------#
    cat("   running remlf90 program...\n\n")
    reml <- data.frame(reml=capture.output(system("remlf90.bat", intern=TRUE)))
    lastRow <- as.character(reml[nrow(reml),])
    if(length(grep("Matrix not positive definite", lastRow))>0){
      warning("\nMatrix not positive definite\n")
    }
    if(file.exists("solutions")){
      renadd <- list.files(pattern="renadd")
      ANIMALpos <- unique(na.omit(as.numeric(unlist(strsplit(unlist(renadd), "[^1-9]+")))))
      Solution <- read.table("solutions", skip = 1, header = F)
      Solution <- Solution[Solution[,2]==ANIMALpos,-c(1,2)]
      if(ncol(Solution)==3){
        colnames(Solution)<-c("animal", "Solution", "se")
      } else {
        colnames(Solution)<-c("animal", "Solution")
      }
      OriginalID <- read.table(renadd, header = F)[,c(1,10)]
      colnames(OriginalID)<-c("animal", "ID")
      Solution <- merge(OriginalID, Solution, by=intersect("animal","animal"))
      write.table(Solution, "Solution.txt", row.names = F, col.names = T, quote = F)
    }
    #---------------------------------------------------------------------------------------#
    if(isTRUE(Inb)){
      #---------------------------------------------------------------------------------------#
      # Reading EBV solution and its s.e.
      #---------------------------------------------------------------------------------------#
      remlf90 <- read.table("remlf90.log", sep="$", colClasses="character")
      Rpos    <- grep("Residual variance", t(remlf90))
      gVAR <- diag(as.matrix(read.table("remlf90.log", skip=3, nrows=nTr)))
      #rVAR <- diag(as.matrix(read.table("remlf90.log", skip=Rpos, nrows=nTr)))
      #renf90 <- read.table("renf90.par", sep="$", colClasses="character")
      renadd <- list.files(pattern="renadd")
      ANIMALpos <- unique(na.omit(as.numeric(unlist(strsplit(unlist(renadd), "[^1-9]+")))))
      Solution <- read.table("solutions", skip = 1, header = F)
      Solution <- Solution[Solution[,2]==ANIMALpos,-c(1,2)]
      if(ncol(Solution)==3){
        cat("Please add: OPTlist=list('OPTION sol se') to compute PEV and Accuracy")
        colnames(Solution)<-c("animal", "Solution", "Se")
        Solution$gVAR <- rep(gVAR, (nrow(Solution)/nTr))
        #---------------------------------------------------------------------------------------#
        # Estimating inbreeding
        #---------------------------------------------------------------------------------------#
        cat("\nComputing Inbreeding...\n")
        inbupgf90 <- data.frame(inbupgf=capture.output(system(paste0("inbupgf90.exe --pedfile ",renadd), intern=TRUE)))
        #---------------------------------------------------------------------------------------#
        Inbreeding <- read.table(paste0(as.character(renadd),".solinb"))
        colnames(Inbreeding) <- c("animal","F")
        Inbreeding$F <-  Inbreeding$F+1
        # Pedigree information and Inbreeding coefficient
        PEDrenum <- read.table(renadd)
        #ped <- PEDrenum[,1:3]
        #colnames(ped) <- c("animal","sire","dam")
        PEDrenum <- PEDrenum[, c(1,ncol(PEDrenum))]
        colnames(PEDrenum) <- c("animal", idName)
        Inbreeding <- merge(PEDrenum, Inbreeding, by=intersect("animal","animal"), all=TRUE)
        Solution <- merge(Solution, Inbreeding, by=intersect("animal","animal"), all=TRUE)
        #---------------------------------------------------------------------------------------#
        # Estimating accuracy
        #---------------------------------------------------------------------------------------#
        # Prediction Error Variance (PEV)
        #---------------------------------------------------------------------------------------#
        cat("   prediction error variance: PEV...\n")
        Solution$PEV  <- Solution$Se^2
        #---------------------------------------------------------------------------------------#
        # Accuracy: sqrt(1-PEV/gVAR)
        # https://gsejournal.biomedcentral.com/articles/10.1186/s12711-016-0188-y
        # http://www-naweb.iaea.org/nafa/news/Blup_Animal_Model.pdf
        #---------------------------------------------------------------------------------------#
        cat("   accuracy: sqrt(1-PEV/gVAR)...\n")
        Solution$Acc  <- suppressWarnings(sqrt(1-(Solution$PEV/Solution$gVAR)))
        #---------------------------------------------------------------------------------------#
        # Accuracy: sqrt(1-PEV/(F*gVAR))
        # http://www.aaabg.org/livestocklibrary/1997/AB97119.pdf
        #---------------------------------------------------------------------------------------#
        cat("   accuracy: sqrt(1-PEV/(F*gVAR))...\n")
        Solution$fAcc <- suppressWarnings(sqrt(1-(Solution$PEV/(Solution$F*Solution$gVAR))))
        #---------------------------------------------------------------------------------------#
        # BIF Accuracy: 1-sqrt(PEV/gVAR)
        # http://beefimprovement.org/content/uploads/2015/08/REVISED-MasterEd-BIF-GuidelinesFinal-08-2015.pdf
        # http://nce.ads.uga.edu/~shogo/html/research/VCE_08.pdf   pag. 36 corrected square to sqrt
        #---------------------------------------------------------------------------------------#
        cat("   BIF accuracy: BIFacc = 1-sqrt(PEV/gVAR)...\n")
        Solution$BIFacc <- suppressWarnings(1-sqrt(Solution$PEV/Solution$gVAR))
        #---------------------------------------------------------------------------------------#
        # Standard Accuracy: sAcc = sqrt(1-(1-BIF)^2)
        # https://www.animalsciencepublications.org/publications/jas/pdfs/92/2/485
        # http://www.beefefficiency.org/info/1297-9686-43-40.pdf
        #---------------------------------------------------------------------------------------#
        cat("   standard accuracy from BIF: sqrt(1-(1-BIF)^2)...\n")
        Solution$sAcc <- suppressWarnings(sqrt(1-((1-Solution$BIFacc)^2)))
        #---------------------------------------------------------------------------------------#
        # The correlation between EPD and true breeding value, rEPD,BV is calculated as:
        # http://beefimprovement.org/content/uploads/2015/08/REVISED-MasterEd-BIF-GuidelinesFinal-08-2015.pdf
        #---------------------------------------------------------------------------------------#
        # rEPD,BV = sqrt((gVAR-PEV)/gVAR)
        #---------------------------------------------------------------------------------------#
        #cat("   true accuracy (r): rEPD,BV = sqrt((gVAR-PEV)/gVAR)...\n")
        #Solution$rEPDxBV <- suppressWarnings(sqrt((Solution$gVAR-Solution$PEV)/(Solution$PEV)))
        #---------------------------------------------------------------------------------------#
        # BIF Accuracy : 1-sqrt(1-rEPD,BV)
        # http://beefimprovement.org/content/uploads/2015/08/REVISED-MasterEd-BIF-GuidelinesFinal-08-2015.pdf
        #---------------------------------------------------------------------------------------#
        #cat("   BIF accuracy 1: rBIF = 1-sqrt(1-rEPD,BV)...\n")
        #Solution$rBIF <- suppressWarnings(1-sqrt(1-Solution$rEPDxBV))
        #Solution$rBIF <- ifelse(Solution$rBIF<0, NA, Solution$rBIF)
        #---------------------------------------------------------------------------------------#
        Solution$Acc   <- ifelse(Solution$Acc   <0, NA, Solution$Acc   )
        Solution$fAcc  <- ifelse(Solution$fAcc  <0, NA, Solution$fAcc  )
        Solution$BIFacc<- ifelse(Solution$BIFacc<0, NA, Solution$BIFacc)
        Solution$sAcc  <- ifelse(Solution$sAcc  <0, NA, Solution$sAcc  )
        #---------------------------------------------------------------------------------------#
        Solution[is.na(Solution)] <- 0
        write.table(Solution, "Sol_Acc.txt", row.names=F, col.names=T, quote=F)
      }
    }
    #---------------------------------------------------------------------------------------#
    # Removing files that is no longer needed anymore
    #---------------------------------------------------------------------------------------#
    cat("\nCleaning...\n\n")
    pattern <- c("rem","ren","phen","ped","map","geno","fsp")
    files <- list.files()
    noRemove <- c("remlf90.log", "remlf90.exe", "renumf90.exe", "Sol_Acc.txt", "renf90.par")
    Remove <- files[grep(paste(pattern,collapse="|"), files)]
    Remove <- Remove[-c(grep(paste(noRemove,collapse="|"), Remove))]
    invisible(file.remove(Remove))
  }
}

