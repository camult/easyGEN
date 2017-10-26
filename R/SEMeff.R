#' @title SEM Effects
#' 
#' @description R functions for summarize effects from strucutal equation models.
#' 
#' @return Mean and SD for each covariate in the SEM model.
#'
#' @export SEMeff
#' 
SEMeff <- function(path=setwd()){
  traits=read.table("renf90.par", header=FALSE, skip=4, nrows =1)[1,1]
  renf90 <- read.table("renf90.par", sep="$", colClasses="character")
  PostMean <- read.table("postmean", sep="$", colClasses="character")
  Rpos    <- grep("RANDOM_RESIDUAL VALUES", t(renf90))
  postR   <- grep("R matrix", t(PostMean))
  Gpos <- grep("VARIANCES", t(renf90))
  postG   <- grep("G matrix for effect", t(PostMean))
  renfSEM <- renf90
  renfSEM[(Gpos+1):(Gpos+traits),] <- PostMean[(postG+1):(postG+traits),]
  renfSEM[(Rpos+1):(Rpos+traits),] <- PostMean[(postR+1):(postR+traits),]
  write.table(renfSEM,  file="renfSEM.par",  row.names=F, col.names=F, quote=F)
  write.table(c("OPTION fixed_var all 1 2"), "renfSEM.par",  row.names=F, col.names=F, quote=F, append=T)
  sink("renfSEM")
  cat("renfSEM.par\n")
  cat(5500, 500,"\n")
  cat(1,"\n")
  sink()
  write.table("gibbs2f90.exe < renfSEM",   file="gibbs2SEM.bat",   row.names=F, col.names=F, quote=F)
  gibbs2SEM <- data.frame(renum=capture.output(system("gibbs2SEM.bat", intern=TRUE)))
  Effects <- read.table("final_solutions", header=F, skip=1)
  Effects <- Effects[Effects[,4]>0,c(2,4,5)]
  colnames(Effects) <- c("Covariate","Mean", "SD")
  return(Effects)
}



