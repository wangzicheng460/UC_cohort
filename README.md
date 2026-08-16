# UC_cohort

溃疡性结肠炎（Ulcerative Colitis, UC）**肠道菌群 — 糖酵解 — PPARG** 多尺度转录组分析

## 项目简介

本研究以溃疡性结肠炎（UC）为疾病模型，整合多个 GEO 黏膜转录组队列，系统分析肠道菌群相关基因、糖酵解通路与 PPARG 在 UC 中的表达特征与诊断价值，并构建机器学习诊断模型。核心结论如下：

- **PPARG** 在 UC 中稳定下调，可作为机制主基因；
- **SLC2A3**（糖酵解分支）在 UC 中稳定上调；
- **LCN2** 是最强的疾病判别基因之一，偏向炎症 / 上皮损伤标志物；
- 构建了由 **8 个稳定基因**组成的 Ridge 诊断模型（锁定外部验证 AUC 0.865–0.957）。

## 目录结构（分析流程）

分析流程按编号组织，从原始数据到下游验证依次进行：

| 目录 | 内容 |
|---|---|
| `1.数据来源/` | 肠道菌群基因、cis-eQTL 与交集基因等输入数据 |
| `2.合并数据+去批次/` | GEO 队列合并、ComBat 去批次（PCA + ComBat） |
| `3.差异分析/` | limma / edgeR 差异表达、火山图、热图 |
| `4.富集分析/` | GO / KEGG 富集分析（clusterProfiler） |
| `5.Venn/` | 差异基因、WGCNA 与糖酵解基因的交集（Venn） |
| `5.WGCNA/` | 加权基因共表达网络分析（WGCNA） |
| `6.ML_reanalysis_20260809/` | 117 个机器学习工作流的候选基因重分析（详见其 README） |
| `6.Venn/` | 肠道菌群基因交集 |
| `7.CIBERSORT/` | CIBERSORT 免疫浸润反卷积与相关性分析 |
| `8.IOBR/` | IOBR 免疫微环境分析（CIBERSORT / MCPcounter / xCell） |
| `sankey_plot/` | 桑基图可视化 |

顶层 `make_fig*.py` / `make_fig*.R` 为论文插图生成脚本（Fig2–Fig6、PPARG 森林图等）。

## 数据来源

- **发现队列（4 个）**：GSE73661、GSE75214、GSE87466、GSE107499，共 414 个样本；
- **外部验证队列（2 个）**：GSE47908、GSE13367；
- **空间转录组**：GSE189184；**单细胞**：GSE214695。

> 原始表达矩阵与大规模中间文件体积过大，未纳入本仓库，请从 GEO 下载后按脚本内的路径组织。

## 依赖环境

### R（建议 4.2+）

`limma`、`edgeR`、`sva`、`clusterProfiler`、`org.Hs.eg.db`、`enrichplot`、`WGCNA`、`glmnet`、`e1071`、`gbm`、`xgboost`、`mboost`、`randomForestSRC`、`pROC`、`IOBR`、`xCell`、`ComplexHeatmap`、`pheatmap`、`ggplot2`、`ggpubr`、`ggrepel`、`ggalluvial`、`ggfortify`、`ggExtra`、`corrplot`、`circlize`、`VennDiagram`、`cowplot`、`patchwork`、`RColorBrewer`、`preprocessCore`、`Seurat`、`harmony`、`CellChat`、`nichenetr`、`TwoSampleMR`、`ieugwasr`、`mr.raps`、`vcfR` 等。

### Python（建议 3.8+）

`numpy`、`pandas`、`matplotlib`、`scipy`、`scikit-learn`、`statsmodels`、`python-docx`、`Pillow`、`pypdfium2`。

## 使用方法

按目录编号顺序执行对应 R 脚本即可复现分析流程；论文插图由顶层 `make_fig*.py` 脚本生成。

## 说明

- 本仓库仅收录代码脚本与必要的小型输入 / 结果表（txt / csv）；图片（pdf / png）、手稿（docx / pptx）、原始测序数据与大规模中间文件已通过 `.gitignore` 排除。
- 机器学习部分的完整结果与解读详见 `6.ML_reanalysis_20260809/README_结果说明.md`。
