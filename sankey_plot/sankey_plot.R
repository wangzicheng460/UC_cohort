setwd("C:\\Users\\wangz\\Desktop\\肠道菌群——完成版\\sankey_plot")
library(ggalluvial)
library(tidyverse)
sankey_data <- tibble(
  microbe = c(
    rep("Enterococcus faecalis", 6),
    rep("Streptococcus salivarius",1)),
  metabolite = c("Agmatine",
                 "Citrulline",
                 "L-Leucic acid",
                 "Lariciresinol",
                 "Levodopa",
                 "Tyramine",
                 "Butyrate"
                 ),
  value = 1)


pdf("sankey_plot.pdf",width = 10,height = 6)
ggplot(sankey_data,
       aes(y = value, axis1 = microbe, axis2 = metabolite)) +
  # 绘制桑基流
  geom_alluvium(aes(fill = microbe),
                width = 1/12, knot.pos = 0.4, alpha = 0.7) +
  # 绘制节点
  geom_stratum(width = 1/12, aes(fill = after_stat(stratum)), color = "black") +
  # 添加节点标签
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3, color = "white", fontface = "bold") +
  # 自定义颜色（匹配原图风格）
  scale_fill_manual(values = c(
    "Enterococcus faecalis" = "#2ECC71",
    "Streptococcus salivarius" = "#9B59B6",
    "Agmatine" = "#A8D5BA",
    "Citrulline" = "#F9D9B1",
    "L-Leucic acid"= "#C5E0F7",
    "Lariciresinol" = "#F4C7C3",
    "Levodopa" = "#D7BDE2",
    "Tyramine" = "#FCE8B2",
    "Butyrate" = "#B0C4D9"
    )) +
  # 设置坐标轴和背景
  scale_x_discrete(limits = c("Microbe", "Metabolite"), expand = c(0.05, 0.05)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )
dev.off()
