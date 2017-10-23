args=commandArgs(trailingOnly=TRUE)

source("https://bioconductor.org/biocLite.R")
#biocLite("hicrep", lib=.libPaths()[1]) #old hicrep
biocLite("rhdf5")

install.packages("rmarkdown",lib=.libPaths()[1],repos="http://cran.rstudio.com/")
install.packages("testthat",lib=.libPaths()[1],repos="http://cran.rstudio.com/")
install.packages("reshape2",lib=.libPaths()[1],repos="http://cran.rstudio.com/")
install.packages("pheatmap",lib=.libPaths()[1],repos="http://cran.rstudio.com/")

install.packages("Supplemental_hicrep_1.0.1.tar.gz",dependencies="logical")