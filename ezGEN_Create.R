library("RcppArmadillo")
library("Rcpp")
library("devtools")
library("roxygen2")

setwd("C:\\Users\\BRITOLOPESF\\Documents\\GitHub\\easyGEN")
roxygenise()


library("Rd2md")
Rd2markdown(rdfile= "C:\\Users\\BRITOLOPESF\\Documents\\GitHub\\easyGEN\\man\\remlf90.Rd",
            outfile="C:\\Users\\BRITOLOPESF\\Documents\\GitHub\\easyGEN\\README.md")
Rd2markdown(rdfile= "C:\\Users\\BRITOLOPESF\\Documents\\GitHub\\easyGEN\\man\\gibbsf90.Rd",
            outfile="C:\\Users\\BRITOLOPESF\\Documents\\GitHub\\easyGEN\\README.md", append=TRUE)
Rd2markdown(rdfile= "C:\\Users\\BRITOLOPESF\\Documents\\GitHub\\easyGEN\\man\\PostGibbs.Rd",
            outfile="C:\\Users\\BRITOLOPESF\\Documents\\GitHub\\easyGEN\\README.md", append=TRUE)
Rd2markdown(rdfile= "C:\\Users\\BRITOLOPESF\\Documents\\GitHub\\easyGEN\\man\\SEMeff.Rd",
            outfile="C:\\Users\\BRITOLOPESF\\Documents\\GitHub\\easyGEN\\README.md", append=TRUE)

setwd("C:\\Users\\BRITOLOPESF\\Documents\\GitHub")
system("R CMD check easyGEN")
system("R CMD build easyGEN")
system("R CMD INSTALL easyGEN_1.01.tar.gz")
