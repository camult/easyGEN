library("RcppArmadillo")
library("Rcpp")
library("devtools")
library("roxygen2")

setwd("C:\\Users\\Camult\\Documents\\GitHub\\easyGEN")
roxygenise()


library("Rd2md")
setwd("C:\\Users\\Camult\\Documents\\GitHub\\ezCV")
Rd2markdown(rdfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\man\\GBLUP_CV.Rd",
            outfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\README.md")
Rd2markdown(rdfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\man\\cvBGBLUP.Rd",
            outfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\README.md", append=TRUE)
Rd2markdown(rdfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\man\\RRBLUP_CV.Rd",
            outfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\README.md", append=TRUE)
Rd2markdown(rdfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\man\\cvBayes.Rd",
            outfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\README.md", append=TRUE)
Rd2markdown(rdfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\man\\MME.Rd",
            outfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\README.md", append=TRUE)
Rd2markdown(rdfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\man\\LOOCV_DG.Rd",
            outfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\README.md", append=TRUE)
Rd2markdown(rdfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\man\\kCV_DG.Rd",
            outfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\README.md", append=TRUE)
Rd2markdown(rdfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\man\\RRBLUP.Rd",
            outfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\README.md", append=TRUE)
Rd2markdown(rdfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\man\\LOOrrDG.Rd",
            outfile="C:\\Users\\Camult\\Documents\\GitHub\\ezCV\\README.md", append=TRUE)

setwd("C:\\Users\\Camult\\Documents\\GitHub")
system("R CMD check ezCV")
system("R CMD build ezCV")
system("R CMD INSTALL ezCV_1.0.tar.gz")
