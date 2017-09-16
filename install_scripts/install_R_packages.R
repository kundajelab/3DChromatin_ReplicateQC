args=commandArgs(trailingOnly=TRUE)

source("https://bioconductor.org/biocLite.R")
biocLite("hicrep", lib=.libPaths()[1])

print(.libPaths(new))
install.packages("reshape2", lib=.libPaths()[1],repos="http://cran.rstudio.com/")
install.packages("pheatmap", lib=.libPaths()[1],repos="http://cran.rstudio.com/")