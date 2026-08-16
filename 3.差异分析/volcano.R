######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#install.packages("ggplot2")


library(ggplot2)           #???ð?
logFCfilter=1              #logFC????????
adj.P.Val.Filter=0.05      #????????pֵ????????
inputFile="all.txt"        #?????ļ?
setwd("D:\\sci\\UC+glycolysis\\3.差异分析")       #???ù???Ŀ¼

#??ȡ?????ļ?
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
#??????????
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")
rt=cbind(rt, Sig=Sig)

#???ƻ?ɽͼ
p=ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("#00A087B2", "#7E6148B2", "#E64B35B2"))+
    xlim(-5,5)+
    labs(title = "volcano plot")+
    geom_vline(xintercept=c(-logFCfilter,logFCfilter), col="grey", cex=1, linetype=2)+
    geom_hline(yintercept= -log10(adj.P.Val.Filter), col="grey", cex=1, linetype=2)+
    theme(plot.title=element_text(size=16, hjust=0.5, face="bold"))
p=p+theme_bw()

#??????ɽͼ
pdf(file="volcano.pdf", width=6, height=4)
print(p)
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

