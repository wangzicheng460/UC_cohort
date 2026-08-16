# 加载必要的包
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(ggrepel)
library(dplyr)
# 示例：加载您的数据
setwd("D:\\fuxian\\UC\\1.5数据合并")
 expr_matrix1 <- read.csv("./merged_data.csv", row.names = 1)
 expr_matrix2 <- read.csv("./GSE87466.csv", row.names = 1)
# 设置随机种子保证结果可重复
set.seed(123)

# 数据预处理
# 1. 过滤低表达基因（例如在少于10%的样本中表达）
keep_genes <- rowSums(expr_matrix1 > 0) >= 0.1 * ncol(expr_matrix1)
filtered_matrix <- expr_matrix1[keep_genes, ]

# 2. 对数转换（加1伪计数避免log(0)）
log_matrixz <- log2(filtered_matrix + 1)

# 3. 中心化和标准化（scale = TRUE表示标准化）
pca_data <- t(log_matrix)  # PCA通常在样本维度进行，所以需要转置

# 执行PCA分析
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

# 查看PCA结果摘要
summary(pca_result)

library(ggfortify)
library(ggplot2)

# 假设 pca_result 是 prcomp() 的结果
# 假设 groups 是长度等于样本数的因子向量（前16实验组，后16对照组）
groups <- factor(c(rep("GSE53306",40 ), rep("GSE87466", 108)))

# 确保 PCA 结果的样本名和 groups 顺序一致
# （如果 expr_matrix 的列名是样本名，且 groups 顺序与之匹配，则无需额外处理）

# 绘制 PCA 得分图 + 椭圆
pdf(file = "PCA.pdf",width = 10,height = 8)
autoplot(pca_result, 
              data = data.frame(Sample = rownames(pca_result$x), Group = groups),
              colour = "Group", 
              size = 3,
              frame = TRUE,
              frame.type = "norm",
              frame.level = 0.95,
              frame.alpha = 0.1) +
  labs(title = "PCA Score Plot",
       color = "Condition") +  # 修改图例标题
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 标题居中，加粗
    plot.margin = margin(10, 10, 10, 10)  # 调整图形边距（可选）
  )

dev.off()

# 可视化2：方差解释比例
variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
variance_df <- data.frame(PC = paste0("PC", 1:length(variance)),
                          Variance = variance,
                          Cumulative = cumsum(variance))

ggplot(variance_df[1:10, ], aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(y = Cumulative, group = 1), color = "red", size = 1) +
  geom_point(aes(y = Cumulative), color = "red", size = 2) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Cumulative Proportion")) +
  labs(title = "Variance Explained by Principal Components",
       x = "Principal Components", y = "Proportion of Variance Explained") +
  theme_minimal()

# 可视化3：PCA载荷图（展示对PC贡献最大的基因）
top_n_genes <- 10  # 展示每个PC中载荷最大的前10个基因

# 提取PC1和PC2的载荷
loadings <- pca_result$rotation[, 1:2]
loadings_df <- as.data.frame(loadings) %>%
  rownames_to_column("Gene") %>%
  arrange(desc(abs(PC1))) %>%
  slice_head(n = top_n_genes)

ggplot(loadings_df, aes(x = PC1, y = PC2, label = Gene)) +
  geom_point(color = "blue") +
  geom_text_repel() +  # 使用ggrepel包避免标签重叠
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "PCA Loading Plot (Top Contributing Genes)",
       x = paste0("PC1 (", round(variance[1]*100, 1), "%)"),
       y = paste0("PC2 (", round(variance[2]*100, 1), "%)")) +
  theme_minimal()

#进行去批次处理
library(limma)
library(sva)
library(dplyr)
#进行的log2转换，也可以进行scale转换
expr_matrix1 <- log2(expr_matrix1+1)
mean_vals <- attr(expr_matrix1, "scaled:center")
sd_vals <- attr(expr_matrix1, "scaled:scale")

# 假设expr_matrix是归一化后的表达矩阵（如log2CPM或vst转换后的数据）
expr_corrected <- ComBat(dat = expr_matrix1, batch = groups)
# 假设log2_data是经过log2转换且加1处理后的数据
expr_unscaled <- 2 ^ expr_corrected - 1
#expr_unscaled <- t(t(expr_corrected) * sd_vals + mean_vals)
expr_unscaled <- as.data.frame(expr_unscaled)
expr_unscaled %>%
  rownames_to_column(var = "Symbol") -> expr_unscaled
write.table(expr_unscaled, "merged.data3.csv", sep = ",", quote = FALSE, row.names = FALSE)
# 重新计算PCA
pca_corrected <- prcomp(t(expr_corrected), scale. = TRUE)
# 按实验分组着色
pdf(file = "PCA2.pdf",width = 8,height = 6)
autoplot(pca_corrected,
                    data = data.frame(Sample = rownames(pca_corrected$x), Group = groups),
                    colour = "Group", size = 3,
                    frame = TRUE, frame.type = "norm", frame.alpha = 0.1) +
  labs(title = "PCA After Batch Correction",
       subtitle = "Batch effect removed using limma::removeBatchEffect") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10))
dev.off()


