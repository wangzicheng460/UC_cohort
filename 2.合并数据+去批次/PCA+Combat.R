BiocManager::install("sva")
# 加载必要的包
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(ggrepel)
library(dplyr)
library(limma)
library(sva)

# 设置工作目录并加载四个数据集
setwd("D:\\sci\\UC+glycolysis\\2.合并数据+去批次")
expr1 <- read.csv("./GSE73661.csv", row.names = 1)  # 假设每个数据集都有基因行名
expr2 <- read.csv("./GSE75214.csv", row.names = 1)
expr3 <- read.csv("./GSE87466.csv", row.names = 1)  # 替换为实际文件名
expr4 <- read.csv("./GSE107499.csv", row.names = 1)  # 替换为实际文件名

# 确保所有数据集使用相同的基因（取交集）
common_genes <- Reduce(intersect, list(rownames(expr1), 
                                       rownames(expr2),
                                       rownames(expr3),
                                       rownames(expr4)))
expr1 <- expr1[common_genes, ]
expr2 <- expr2[common_genes, ]
expr3 <- expr3[common_genes, ]
expr4 <- expr4[common_genes, ]

# 合并数据集并创建批次信息
combined_expr <- cbind(expr1, expr2, expr3, expr4)
batch_info <- factor(c(rep("GSE73661", ncol(expr1)),
                       rep("GSE75214", ncol(expr2)),
                       rep("GSE87466", ncol(expr3)),
                       rep("GSE107499", ncol(expr4))))

# 设置随机种子
set.seed(123)

# 数据预处理
# 1. 过滤低表达基因（在至少10%样本中表达）
keep_genes <- rowSums(combined_expr > 0) >= 0.1 * ncol(combined_expr)
filtered_matrix <- combined_expr[keep_genes, ]

# 2. 对数转换（加1伪计数）
log_matrix <- log2(filtered_matrix + 1)

# ============== 批次效应前PCA ==============
pca_raw <- prcomp(t(log_matrix), center = TRUE, scale. = TRUE)

# 绘制原始PCA
pdf("PCA_Raw.pdf", width = 6, height = 4)
autoplot(pca_raw, 
         data = data.frame(Sample = colnames(log_matrix), Batch = batch_info),
         colour = "Batch", 
         size = 3,
         frame = TRUE,
         frame.type = "norm",
         alpha=0.4) +
  labs(title = "PCA Before Batch Correction",
       subtitle = "Colored by Dataset Source") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# ============== 批次效应校正 ==============
# 使用ComBat进行批次校正
combat_corrected <- ComBat(dat = as.matrix(log_matrix), 
                           batch = batch_info,
                           par.prior = TRUE)

# 将校正后的数据转换回数据框
corrected_df <- as.data.frame(combat_corrected)

expr_unscaled <- 2 ^ corrected_df - 1
# 保存校正后的数据
expr_unscaled %>%
  rownames_to_column(var = "Gene") %>%
  write.csv("Batch_Corrected_Data.csv", row.names = FALSE)



# ============== 批次效应后PCA ==============
pca_corrected <- prcomp(t(combat_corrected), center = TRUE, scale. = TRUE)

# 绘制校正后PCA
pdf("PCA_Corrected.pdf", width = 6, height = 4)
autoplot(pca_corrected, 
         data = data.frame(Sample = colnames(combat_corrected), Batch = batch_info),
         colour = "Batch", 
         size = 3,
         frame = TRUE,
         frame.type = "norm",
         alpha=0.5) +
  labs(title = "PCA After Batch Correction (ComBat)",
       subtitle = "Colored by Dataset Source") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# ============== 方差解释图 ==============
variance <- pca_corrected$sdev^2 / sum(pca_corrected$sdev^2)
variance_df <- data.frame(PC = paste0("PC", 1:length(variance)),
                          Variance = variance,
                          Cumulative = cumsum(variance))

pdf("Variance_Explained.pdf", width = 8, height = 6)
ggplot(variance_df[1:10, ], aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(y = Cumulative, group = 1), color = "red", size = 1) +
  geom_point(aes(y = Cumulative), color = "red", size = 2) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Cumulative Proportion")) +
  labs(title = "Variance Explained by Principal Components",
       x = "Principal Components", 
       y = "Proportion of Variance Explained") +
  theme_minimal()
dev.off()