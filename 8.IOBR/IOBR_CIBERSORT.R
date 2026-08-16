library(IOBR)
setwd("D:\\sci\\UC+glycolysis\\18.IOBR")
inputFile="normalize.txt" 
exp=read.table(inputFile, header=T, sep="\t", check.names=F,row.names = 1)
res_cibersort<- deconvo_tme(eset = exp,
                            method ="cibersort",
                            arrays = FALSE,
                            perm =1000)
write.table(res_cibersort, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)
