#install.packages("reshape2")
#install.packages("ggpubr")
#install.packages("corrplot")


#???ð?
library(reshape2)
library(ggpubr)
library(corrplot)

inputFile="CIBERSORT-Results.txt"     #????ϸ???????Ľ????ļ?
setwd("D:\\sci\\UC+glycolysis\\14.CIBERSORT")     #???ù???Ŀ¼

#??ȡ????ϸ???????ļ?
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#????Ʒ???з???(????????ʵ????)
con=grepl("_con", rownames(rt), ignore.case=T)
treat=grepl("_treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData, treatData))

#??????״ͼ
pdf(file="barplot.pdf", width=13, height=7)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,xaxt="n",yaxt="n",ylab="Relative Percent",cex.lab=1.5)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
#??ͼ???б?ע??????????Ʒ,????ɫ??ʾ
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="#008B45FF")
text(a1[conNum]/2,-0.035,"Control",cex=1.0)
#??ͼ???б?עʵ????????Ʒ,?ú?ɫ??ʾ
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="#EE0000FF")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"Treat",cex=1.0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1)
dev.off()

##################????????ͼ##################
#??????ת????ggplot2?????ļ?
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=cbind(as.data.frame(t(data)), Type)
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
#????????ͼ
group=levels(factor(data$Type))
bioCol=c("#5D90BA","#FABFD2","#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Type",
				  xlab="",
				  ylab="Fraction",
				  legend.title="Type",
				  #notch=T, add="point",
				  width=0.8,
				  palette=bioCol)+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#????????ͼ
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()

##################??????ͼ##################
treatData=treatData[,apply(treatData,2,sd)>0]
pdf(file="corHeatmap.pdf", width=12, height=12)
corrplot(corr=cor(treatData, method="spearman"),
         method = "color",        #ͼ?ε?չʾ??ʽ
         order = "hclust",        #????ϸ??????????ʽ
         tl.col="black",          #??????ɫ
         number.cex = 0.8,        #????ϵ????????С
         addCoef.col = "black",   #????ϵ????????ɫ
         col=colorRampPalette(c("blue", "white", "red"))(50),    #??ɫ??????
         )
dev.off()


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio


