library(ggplot2)
library(ggrepel)

setwd("C:/Users/wangz/Desktop/肠道菌群——完成版/3.差异分析")

logFCfilter <- 1
adj.P.Val.Filter <- 0.05
inputFile <- "all.txt"

rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)

Sig <- ifelse((rt$adj.P.Val < adj.P.Val.Filter) & (abs(rt$logFC) > logFCfilter),
              ifelse(rt$logFC > logFCfilter, "Up", "Down"), "Not")
rt <- cbind(rt, Sig = Sig)

labelGenes <- c("LCN2", "PPARG", "SLC2A3")
rt$label <- ifelse(rt$id %in% labelGenes, rt$id, NA)

# low-saturation, colorblind-friendly palette
colUp   <- "#C44E52"  # muted red
colDown <- "#4C72B0"  # muted blue
colNot  <- "#B8B8B8"  # grey

p <- ggplot(rt, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Sig), size = 1.2, alpha = 0.6) +
  scale_color_manual(values = c("Up" = colUp, "Down" = colDown, "Not" = colNot)) +
  geom_hline(yintercept = -log10(adj.P.Val.Filter), linetype = 2, color = "grey45") +
  geom_vline(xintercept = c(-logFCfilter, logFCfilter), linetype = 2, color = "grey45") +
  geom_point(data = subset(rt, id %in% labelGenes),
             aes(x = logFC, y = -log10(adj.P.Val)),
             shape = 21, fill = "white", color = "black", size = 2.8, stroke = 1) +
  geom_text_repel(data = subset(rt, id %in% labelGenes),
                  aes(x = logFC, y = -log10(adj.P.Val), label = label),
                  size = 4, color = "black", fontface = "italic",
                  box.padding = 0.6, point.padding = 0.4,
                  min.segment.length = 0, segment.color = "grey30",
                  max.overlaps = Inf) +
  labs(x = "log2 Fold Change", y = expression(-log[10]~(adjusted~P~value))) +
  xlim(-5, 5) +
  theme_bw(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.key = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92")
  )

ggsave("volcano_annotated.pdf", p, width = 6, height = 5, device = cairo_pdf)
ggsave("volcano_annotated.png", p, width = 6, height = 5, dpi = 300)
