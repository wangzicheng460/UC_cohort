#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("tidyverse")

#install.packages("devtools")
#devtools::install_github("Hy4m/linkET", force = TRUE)


#???ð?
library(limma)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(linkET)

expFile="normalize.txt"       #?????????ļ?
geneFile="modelGene.list.txt"       #?????б??ļ?
immFile="CIBERSORT-Results.txt"     #????ϸ???????Ľ????ļ?
setwd("D:\\sci\\肠道菌群\\14.CIBERSORT")     #???ù???Ŀ¼

#??ȡ?????????ļ?,?????????ļ?????????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#??ȡ?????б??ļ?, ??ȡģ?ͻ????ı???��
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

#, drop=FALSEȥ??????????Ʒ
groucolnames(data)
p=gsub("(datagroup=gsu.*)\\_(.*?)", "\\", colnames(data))
data=data[,group=="treat",drop=F]
data=t(data)

#??ȡ????ϸ???????ļ??????????ݽ???????
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
#?Ա??????ݺ?????ϸ??????ȡ????
sameSample=intersect(row.names(data), row.names(immune))
data=data[sameSample,,drop=F]
immune=immune[sameSample,,drop=F]
immune=immune[,apply(immune,2,sd)>0]
geneLists=list()
for(i in 1:ncol(data)){geneLists[[colnames(data)[i]]]=i}

#??????????ϸ???????Է???
geneCor=data.frame()
for(cell in colnames(immune)){
	if(sd(immune[,cell])==0){next}
	for(gene in colnames(data)){
		x=as.numeric(immune[,cell])
		y=as.numeric(data[,gene])
		corT=cor.test(x, y, method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		geneCor=rbind(geneCor, cbind(spec=gene, env=cell, r=cor, p=pvalue))
	}
}
geneCor$r=as.numeric(geneCor$r)
geneCor$p=as.numeric(geneCor$p)
geneCor$pd=ifelse(geneCor$p<0.05, ifelse(geneCor$r>0, "Postive", "Negative"), "Not")
geneCor$r=abs(geneCor$r)
geneCor=geneCor %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6",">= 0.6")))

#????ͼ??
qcorPlot=qcorrplot(correlate(immune, method="spearman"), type = "lower", diag = FALSE) +
  geom_square() +
  #geom_mark(sep = '\n', size = 1, sig_level = c(0.05, 0.01, 0.001), sig_thres = 0.05, color="black") +
  geom_couple(aes(colour = pd, size = rd), 
              data = geneCor, 
              curvature = nice_curvature()) +
  #????ͼ?ε???ɫ??ͼ????????
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))) +
  scale_size_manual(values = c(0.5, 1.5, 2, 3)) +
  scale_colour_manual(values = c("#1B9E77", "#CCCCCC99", "#D95F02")) +
  guides(size = guide_legend(title = "abs(Cor)",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "pvalue", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Cell-cell cor", order = 3))

#????ͼ??
pdf(file="linkET.pdf", width=9, height=7)
print(qcorPlot)
dev.off()


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

