library(tidyverse)
library(ggplot2)
library(reshape2)
library(corrplot)
library(xCell)
#devtools::install_github("dviraran/xCell", force = TRUE)
setwd("D:/sci/肠道菌群/18.IOBR")
exp=read.table("normalize.txt", header=T, sep="\t", check.names=F,row.names = 1)
xcell<-xCellAnalysis(exp,rnaseq=TRUE)
write.table(xcell,file="xCell.txt",row.names=T,col.names=T,sep="\t",quote=F)
genelist<-c("PPARG")
goal_exp<-filter(exp,rownames(exp) %in%genelist)
#合并目标基因集表达矩阵和免疫细胞矩阵
combine<-rbind(goal_exp,xcell)
#计算相关系数
comcor<-cor(t(combine))
#计算显著性差异
comp<-cor.mtest(comcor,conf.level=0.95)
pval<-comp$p
#获取目标基因相关性矩阵
goalcor<-select(as.data.frame(comcor),genelist)%>%
  rownames_to_column(var="celltype")
goalcor<-filter(goalcor,!(celltype %in% genelist))
##长宽数据转换
goalcor<-melt(goalcor,id.vars="celltype")
colnames(goalcor)<-c("celltype","Gene","correlation")
#获取目标基因集pvalue矩阵
pval<-select(as.data.frame(pval),genelist)%>%
  rownames_to_column(var="celltype")
pval<-filter(pval,!(celltype %in% genelist))
#长宽数据转换
pval<-melt(pval,id.vars="celltype")
colnames(pval)<-c("celltype","gene","pvalue")#将pvalue和correlation两个文件合并
final<-left_join(goalcor,pval,by=c("celltype"="celltype","Gene"="gene"))
write.table(final,file="xCell.result.txt",sep="\t",quote=F,row.names = F)
#添加一列,来判断pvalue值范围
final$sign<-case_when(final$pvalue<0.05 &final$pvalue>0.01 ~"*",
                      final$pvalue<0.01 &final$pvalue>0.001 ~"**", 
                      final$pvalue<0.001 ~"***",
                      final$pvalue>0.05 ~"")
ggplot(data=final,
       aes(x=Gene,y=celltype))+geom_tile(aes(fill=correlation),colour="white",size=1)+scale_fill_gradient2(low="#2b8cbe",mid="white",high="#e41a1c")+geom_text(aes(label=sign),colour="black")+theme_minimal()+theme(axis.text.x=element_text(angle=45,hjust=1,size=12),      axis.text.y=element_text(size=12),    axis.title.x=element_blank(),    axis.title.y=element_blank(),    axis.ticks.x=element_blank(),    axis.ticks.y=element_blank()) +guides(fill=guide_legend(title="* p<0.05\n\n** p<0.01\n\n*** p<0.001\n\ncorrelation"))
ggsave("correlation.pdf",width=6,height=13)

Epithelial_score <- xcell["Epithelial cells", ]
PPARG_expr <- as.numeric(exp["PPARG", ])

class(Epithelial_score)
# 3. 相关性分析
cor.test(PPARG_expr, Epithelial_score, method="spearman")
library(ggplot2)
library(ggExtra)
library(ggpubr)
# 4. 可视化
# 创建基础散点图
p<-ggplot(data.frame(PPARG = PPARG_expr, Epithelial = Epithelial_score), 
            aes(x = PPARG, y = Epithelial)) +
  geom_point(size = 3, alpha = 0.7, color = "#E64B35") +
  geom_smooth(method = "lm", se = FALSE, color = "#4DBBD5", linewidth = 1.2) +
  labs(title = "PPARG Expression vs Epithelial Cell Infiltration",
       x = "PPARG Expression (log2TPM)",
       y = "Epithelial Cell Abundance Score") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 添加相关系数标注（自动计算R和p值）
p <- p + stat_cor(
  method = "pearson",           # 相关系数类型
  label.x = min(PPARG_expr),     # 标签位置（左下角）
  label.y = max(Epithelial_score),  
  size = 5,                     # 字体大小
  color = "black",              # 字体颜色
  show.legend = FALSE
)

# 添加山脊图（边缘密度分布）
p<-ggMarginal(
  p, 
  type = "density",             # 密度曲线
  fill = "#E64B35",             # 填充色（匹配散点）
  alpha = 0.3,                  # 透明度
  color = NA,                   # 边框颜色（无）
  size = 4                      # 调整山脊图大小
)
ggsave("Epithelial_correlation.pdf",plot = p,height = 6,width = 8)
