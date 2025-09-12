#!/usr/bin/env R
# coding utf-8

###################### Package Loading ###########################
pacman::p_load(
    "Seurat", "Nebulosa", "ggplot2", "future", "reshape2",
    "SingleCellExperiment", "dplyr", "tidyverse", "ggrepel",
    "patchwork", "msigdbr", "GSVA", "RColorBrewer", "ggpubr",
    "ROGUE", "viridis", "magrittr", "data.table",
    "R.utils", "grid", "cowplot", "tidyverse", "dittoSeq",
    "harmony", "scRepertoire", "ggsci", "pheatmap", "ggpie",
    "sscVis", "alakazam", "UpSetR", "CytoTRACE", "ggforce",
    "tidydr", "ggplotify", "doParallel", "readr", "future",
    "monocle", "clusterProfiler", "org.Hs.eg.db", "oposSOM",
    "scrat", "magick", "ggnewscale", "destiny", "ggthemes",
    "diffusionMap", "slingshot", "ggpointdensity", "ggradar",
    "scales", "plotly", "processx", "batchelor", "SeuratWrappers",
    "vegan", "forcats", "Nebulosa", "purrr", "Matrix", "pvclust",
    "ComplexHeatmap", "circlize", "stringr", "tibble", "SeuratDisk",
    "ggsignif", "bseqsc", "fgsea", "ggpmisc", "zoo", "nichenetr"
)
##################################################################


##################################################################
## loading object
setwd("/work/xiaxy/work_2024/wangll_blood_final1/output2/")
data <- readRDS("input_file/input.rds")
mye <- readRDS("input_file/myeloid_sub.rds")
bp <- readRDS("input_file/bp_sub.rds")
nkt <- readRDS("input_file/nkt_sub.rds")
t_data <- readRDS("input_file/tissue_input.rds")

bp_bcr <- readRDS("input_file/bcr_processed_data.rds")
bp01_bcr <- readRDS("input_file/bp01_bcr.rds")
bp06_bcr <- readRDS("bp06_bcr.rds")
bp_tdata <- readRDS("tissue_result/tissue_bp.rds")
bp_pan <- readRDS("../blood_mao/b_cell/b_blood_mao_input.rds")
load("input_file/color.RData")
##################################################################


##################################################################
## Figure 7E & 7F
## 1) 基因清单：仅保留存在于对象中的基因，并保持原顺序
genes_all <- c(
    "FKBP5", "LPCAT1", "PPP3CA", "PPP3CC", "PPP3R1", "NFATC2", "NFATC3",
    "REL", "MAPK1", "MAPK14", "MAP3K2", "MAP3K5",
    "RASGRP1", "RB1CC1", "TAOK3", "TNIK", "CCDC88C", "PLCB1", "PTK2B",
    "ATM", "UVRAG", "WAC", "IFI16", "PIP4K2A",
    "RPTOR", "VMP1", "VPS13A", "VPS13B", "VPS13C", "VPS13D",
    "BCL2", "ATG7", "NLRC3", "NLRC5",
    "DOCK8", "DOCK10", "CXCR4", "FOS", "FOSB",
    "NR4A2", "JUN", "JUNB", "SLC8A1", "ZEB2"
)

## 2) 取 CAMR 子集并只保留目标三类细胞
data_sub <- subset(data, group == "CAMR")
mye_sub <- subset(mye, group == "CAMR")
bp_sub <- subset(bp, group == "CAMR")
nkt_sub <- subset(nkt, group == "CAMR")

cells_keep <- c(colnames(mye_sub), colnames(bp_sub), colnames(nkt_sub))
data_sub <- subset(data_sub, cells = cells_keep)

## 3) 基因表达矩阵：按样本 orig.ident 取平均（slot = "data"）
genes_in_obj <- genes_all[genes_all %in% rownames(data_sub)]
if (length(genes_in_obj) == 0) stop("None of the requested genes are present in data_sub.") # nolint

avg_list <- AverageExpression(
    data_sub,
    assays = "RNA", slot = "data",
    features = genes_in_obj, group.by = "orig.ident", verbose = FALSE
)
mat <- as.matrix(avg_list$RNA) # 行 = 基因, 列 = 样本
# 若有缺失基因，补 NA（可选）
mat_full <- matrix(NA_real_,
    nrow = length(genes_all), ncol = ncol(mat),
    dimnames = list(genes_all, colnames(mat))
)
mat_full[rownames(mat), colnames(mat)] <- mat
mat <- mat_full

## 4) 细胞比例：行=样本，列=亚群；按行归一化
# t1: myeloid
t1_raw <- table(mye_sub$orig.ident, mye_sub$anno2) # 行名=样本, 列=anno2
t1 <- prop.table(t1_raw, 1) %>% as.matrix()

# t2: nkt
t2_raw <- table(nkt_sub$orig.ident, nkt_sub$anno2)
t2 <- prop.table(t2_raw, 1) %>% as.matrix()

# t3: B/P
t3_raw <- table(bp_sub$orig.ident, bp_sub$anno2)
t3 <- prop.table(t3_raw, 1) %>% as.matrix()

## 只保留你关心的亚群（请把下面这些名字改成你对象里真正的列名）
target_myeloid <- c("m06_ZEB2+ monocyte") # 例：你原来取 t1[,6]
target_nkt <- c("tn03_ANK3+CD4+ Tm", "tn10_ZEB2+ T") # 例：你原来取 t2[,c(3,9)]
target_bp <- c("bp01.Pro-BC") # 例：你原来取 t3[,1]

# 防御式按列名选择（缺失列会被静默跳过/或报错，视需求而定）
sel_cols <- function(M, cols) {
    cols_ok <- intersect(cols, colnames(M))
    if (length(cols_ok) == 0) {
        warning("None of target columns found in matrix.")
        return(matrix(numeric(0), nrow = nrow(M), dimnames = list(rownames(M), NULL))) # nolint
    }
    M[, cols_ok, drop = FALSE]
}

t1_sel <- sel_cols(t1, target_myeloid)
t2_sel <- sel_cols(t2, target_nkt)
t3_sel <- sel_cols(t3, target_bp)

fra <- cbind(t1_sel, t2_sel, t3_sel)

## 5) 载入临床，并对齐样本 + 标准化
cl <- read.table("input_file/p_clinical.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE) # nolint
cl <- cl[1:18, c(1:10, 16:19, 21)] # 原选择
# 仅保留与表达/比例都有交集的样本
common_samples <- Reduce(intersect, list(colnames(mat), rownames(fra), rownames(cl))) # nolint
mat <- mat[, common_samples, drop = FALSE]
fra <- fra[common_samples, , drop = FALSE]
cl <- cl[common_samples, , drop = FALSE]

## 6) 标准化（按行 z-score）
gene_scaled <- t(scale(t(mat)))
cell_scaled <- t(scale(fra))
cli_scaled <- t(scale(cl))

## 7) 合并基因+细胞比例用于样本聚类
combined_mat <- rbind(gene_scaled, cell_scaled) # 注意 cell_scaled 先转置成 行=特征，列=样本 # nolint

## 8) 颜色
gene_col <- colorRamp2(c(-1, 0, 1), c("#4575B4", "#FFFFBF", "#D73027"))
clinical_col <- colorRamp2(seq(-1, 1, length = 3), c("#F0F0F0", "#BDBDBD", "#252525")) # nolint

## 9) 基因+细胞比例热图（KM 列分组，并保存列顺序）
ht_gene <- Heatmap(
    combined_mat,
    name = "Gene\nExpression",
    col = gene_col,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    show_row_names = TRUE,
    row_title = paste0(nrow(gene_scaled), " Genes + ", ncol(cell_scaled), " Cells"), # nolint
    clustering_method_columns = "complete",
    row_title_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 8),
    column_km = 2,
    row_km = 2,
    heatmap_legend_param = list(title_position = "leftcenter-rot", legend_height = unit(3, "cm")) # nolint
)

## 10) 临床热图按相同顺序绘制
ht_gene_drawn <- draw(ht_gene)
ordered_columns <- column_order(ht_gene_drawn)

ht_clinical <- Heatmap(
    cli_scaled[, as.numeric(unlist(ordered_columns))],
    name = "Clinical\nFeatures",
    col = clinical_col,
    cluster_columns = FALSE, # 不再聚类
    cluster_rows = TRUE,
    show_row_names = FALSE,
    row_names_side = "left",
    row_title = "10 Clinical",
    column_split = c(rep("group1", 6), rep("group2", 12)),
    row_title_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(
        title_position = "leftcenter-rot",
        legend_height = unit(3, "cm")
    )
)

## 11) 组别与统计（与箱线图一致）
cl$group <- ifelse(rownames(cl) %in% c("CAR26_P", "CAR02_P", "CAR34_P", "CAR18_P", "CAR01_P", "CAR04_P"), "ap1", "au") # nolint

apply(cl[, setdiff(colnames(cl), "group"), drop = FALSE], 2, function(v) {
    wilcox.test(v[cl$group == "au"], v[cl$group == "ap1"])
})

## 12) 箱线图（和你原先一致，但更稳）
p_box <- function(y) {
    ggplot(cl, aes(x = group, y = .data[[y]], color = group)) +
        geom_boxplot(width = 0.6, outlier.colour = "white", outlier.size = 0) +
        scale_color_aaas() +
        theme_classic() +
        theme(axis.text = element_blank(), legend.position = "none", axis.title = element_blank()) # nolint
}

p2 <- p_box("ct")
p3 <- p_box("URPO") + scale_y_log10()
p4 <- p_box("eGFR")
p5 <- p_box("Crea")

final <- p4 + plot_spacer() + p3 + plot_spacer() + p2 + plot_spacer() + p5 +
    plot_layout(nrow = 1, widths = c(1, .1, 1, .1, 1, .1, 1, .1, 1))

## 13) 输出热图 PDF
ht_list <- ht_gene %v% ht_clinical
pdf("Figure7/Figure_7E_hclust.pdf", width = 7, height = 9)
draw(ht_list,
    column_title = "Integrated Multi-omics Analysis",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    padding = unit(c(2, 10, 2, 2), "mm")
) # bottom, left, top, right
dev.off()

## 可选：保存箱线图组合
ggsave("Figure7/Figure_7F_clinical_boxes.pdf", final, width = 12, height = 2.2)
##################################################################


##################################################################
## Figure 7D
## 1) 子集与分母（每样本在 data_sub 的总细胞数）
data_sub <- subset(data, cells = c(colnames(mye), colnames(nkt), colnames(bp)))
denom <- table(data_sub$orig.ident) # 命名向量：names=样本
denom <- as.numeric(denom)
names(denom) <- names(table(data_sub$orig.ident))

## 2) 计数函数：行=样本、列=anno2；再用 data_sub 的分母去除
tab_frac_vs_data_sub <- function(obj, id_col = "orig.ident", type_col = "anno2", denom_vec) { # nolint
    T <- table(obj@meta.data[[id_col]], obj@meta.data[[type_col]]) # 行=样本
    T <- as.matrix(T)
    # 补齐分母中存在而该子集缺失的样本（补0行）
    missing_samples <- setdiff(names(denom_vec), rownames(T))
    if (length(missing_samples) > 0) {
        pad <- matrix(0,
            nrow = length(missing_samples), ncol = ncol(T),
            dimnames = list(missing_samples, colnames(T))
        )
        T <- rbind(T, pad)
    }
    # 只保留在分母中的样本，并按分母顺序排列
    T <- T[names(denom_vec), , drop = FALSE]
    # 用 data_sub 的样本总数做分母（逐行相除）
    frac <- sweep(T, 1, denom_vec, "/")
    return(frac)
}

## 3) 三个子集的“占 data_sub 总细胞数”的比例矩阵
mye_mat <- tab_frac_vs_data_sub(mye, denom_vec = denom)
bp_mat <- tab_frac_vs_data_sub(bp, denom_vec = denom)
nkt_mat <- tab_frac_vs_data_sub(nkt, denom_vec = denom)

## 4) 合并并对齐（行=样本；列=各亚群；缺失样本已补为0）
mt <- cbind(mye_mat, bp_mat, nkt_mat)
colnames(mt) <- make.unique(colnames(mt)) # 避免重名列

## 5) 读取临床，并与 mt 对齐
cl <- read.table("input_file/p_clinical.txt",
    header = TRUE, sep = "\t",
    row.names = 1, check.names = FALSE
)

common <- intersect(rownames(mt), rownames(cl))
stopifnot(length(common) > 0)

mt <- mt[common, , drop = FALSE]
cl <- cl[common, , drop = FALSE]

## 6) 相关性：示例列名请确认存在
stopifnot("m06_ZEB2+ monocyte" %in% colnames(mt))
stopifnot("m09_pDC" %in% colnames(mt))
stopifnot(all(c("g", "PTC", "gs", "cv") %in% colnames(cl)))

d1 <- data.frame(fra = mt[, "m06_ZEB2+ monocyte"], g = cl$g)
p <- ggplot(d1, aes(x = fra, y = g)) +
    geom_point() +
    geom_smooth(method = "lm", linewidth = 0.2, color = "black") +
    stat_cor() +
    theme_classic() +
    theme(
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.title = element_blank(),
        legend.position = "none"
    )

d2 <- data.frame(fra = mt[, "m09_pDC"], ptc = cl$PTC, gs = cl$gs, cv = cl$cv)

mk_scatter <- function(df, y) {
    ggplot(df, aes(x = fra, y = .data[[y]])) +
        geom_point() +
        geom_smooth(method = "lm", linewidth = 0.2, color = "black") +
        stat_cor() +
        theme_classic() +
        theme(
            axis.line = element_line(linewidth = 0.4),
            axis.ticks = element_line(linewidth = 0.4),
            axis.title = element_blank(),
            legend.position = "none"
        )
}

p1 <- mk_scatter(d2, "ptc")
p2 <- mk_scatter(d2, "gs")
p3 <- mk_scatter(d2, "cv")

final <- p + p1 + p2 + p3 + plot_layout(nrow = 2)
final
ggsave("Figure7/Figure_7D_fra_vs_clinical_datasub.pdf", final, width = 8, height = 6) # nolint
##################################################################


##################################################################
## Figure 7A
t_data$anno4 <- factor(t_data$anno4, levels = rev(sort(unique(t_data$anno4))))
t_data_sub1 <- subset(t_data, group == "CAMR")
t_data_sub2 <- subset(t_data, group == "HC")

p1 <- DotPlot(t_data_sub1, features = c("CXCR4", "CXCL12"), scale.max = 80, col.min = -0.5, col.max = 2.5, cols = c("lightgrey", "red"), group.by = "anno4") + # nolint
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
p2 <- DotPlot(t_data_sub2, features = c("CXCR4", "CXCL12"), scale.max = 80, col.min = -0.5, col.max = 2.5, cols = c("lightgrey", "red"), group.by = "anno4") + # nolint
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
final <- p1 + p2 + plot_layout(nrow = 1)
ggsave("Figure7/Figure_7A_DotPlot.pdf", final, width = 9, height = 4.5)
##################################################################


##################################################################
## Figure 7B
## 1) 样本×亚类 计数表
T <- table(t_data$orig.ident, t_data$anno4)

## 2) 用 t_data 的样本总细胞数作分母（行比例）
frac <- prop.table(T, 1) # 行=样本；列=anno4；值=该anno4/该样本在 t_data 的总细胞数

## 3) 转成长表，并合并样本的组别（CAMR/HC）
df <- as.data.frame(as.table(frac))
colnames(df) <- c("orig.ident", "anno4", "Freq")

# 从元数据里拿样本组别（确保每个样本只有一个组别）
grp <- t_data@meta.data |>
    distinct(orig.ident, group)
df <- df |> left_join(grp, by = "orig.ident")

## 4) 只画你关心的两个亚类 —— 用“名字”而不是列号
# 可先查看：colnames(frac)
targets <- c("c06_Fibroblast", "c13_Podocyte") # <- 把实际列名替换进去
df_sub <- df |>
    filter(anno4 %in% targets) |>
    mutate(
        anno4 = factor(anno4, levels = targets),
        group = factor(group, levels = c("CAMR", "HC"))
    )

## 5) 画图（按组并排箱线）
p <- ggplot(df_sub, aes(x = anno4, y = Freq, color = group)) +
    geom_boxplot(outlier.size = 0, outlier.alpha = 0, position = position_dodge(width = 0.7)) + # nolint
    scale_color_manual(values = mycolor$group) +
    scale_y_continuous(limits = c(0, 0.03)) +
    theme_classic() +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3)
    )
ggsave("Figure7/Figure7B_fraction_in_tissue.pdf", p, width = 2, height = 3)
##################################################################


##################################################################
## Figure 7C
## 1) 样本×anno4 计数，并用每样本在 t_data 的总细胞数作分母（行比例）
T <- table(t_data$orig.ident, t_data$anno4) # 行=样本, 列=anno4
frac <- prop.table(T, 1) # 每行除以该样本总数（正是你想要的分母）

## 2) 转成长表，并从元数据按样本并入 group（不要用 rep(...)）
df <- as.data.frame(as.table(frac))
colnames(df) <- c("orig.ident", "anno4", "value")

grp <- t_data@meta.data %>% distinct(orig.ident, group)
df <- df %>% left_join(grp, by = "orig.ident")

## 3) 只保留前 5 个大类（用名字而不是 1:5 列号）
cats <- c("c01_T cell", "c02_NK cell", "c03_B cell", "c04_Plasma cell", "c05_Myeloid") # nolint
df5 <- df %>%
    filter(anno4 %in% cats) %>%
    mutate(
        anno4 = factor(anno4, levels = cats),
        group = factor(group, levels = c("CAMR", "HC"))
    )

## 4) 一个小函数，避免重复代码；按需传 y 轴上限
mk_box <- function(dat, cat, ylim = NULL, add_p = FALSE) {
    p <- ggplot(filter(dat, anno4 == cat), aes(x = group, y = value, color = group)) + # nolint
        geom_boxplot(outlier.size = 0, outlier.alpha = 0, width = 0.6) +
        scale_color_manual(values = mycolor$group) +
        theme_classic() +
        theme(
            axis.title = element_blank(),
            legend.position = "none",
            axis.line = element_line(linewidth = 0.3),
            axis.ticks = element_line(linewidth = 0.3)
        )
    if (!is.null(ylim)) p <- p + scale_y_continuous(limits = ylim)
    if (add_p) p <- p + stat_compare_means()
    p
}

p1 <- mk_box(df5, "c01_T cell", ylim = c(0, 0.6), add_p = FALSE)
p2 <- mk_box(df5, "c02_NK cell", ylim = NULL, add_p = FALSE) # 你原来未设上限
p3 <- mk_box(df5, "c03_B cell", ylim = c(0, 0.05), add_p = TRUE)
p4 <- mk_box(df5, "c04_Plasma cell", ylim = c(0, 0.03), add_p = TRUE)
p5 <- mk_box(df5, "c05_Myeloid", ylim = c(0, 0.25), add_p = TRUE)

final <- p1 + plot_spacer() + p2 + plot_spacer() + p3 + plot_spacer() + p4 + plot_spacer() + p5 + # nolint
    plot_layout(nrow = 1, widths = c(1, 0.15, 1, 0.15, 1, 0.15, 1, 0.15, 1))

ggsave("Figure7/Figure_7C_tissue_immune.pdf", final, width = 15, height = 3)
##################################################################


##################################################################
## Figure 8D
data_exp <- readRDS("input_file/pan_blood_sub13Gene.Rds")
genes <- c(
    "IL6", "TNF", "IL1A", "IL1B", "IL12A",
    "IL12B", "IL18", "JAK1", "JAK2", "JAK3",
    "FKBP1A", "FKBP5", "PPP3CA"
)

kk <- as.character(unique(data_exp$anno))
pvalue <- matrix(NA, ncol = 3 * length(kk), nrow = length(genes)) # nolint
fc <- matrix(NA, ncol = 3 * length(kk), nrow = length(genes)) # nolint

for (i in 1:length(kk)) {
    print(i)
    sub <- subset(data_exp, anno == kk[i])
    for (j in 1:length(genes)) {
        print(j)
        oa <- subset(sub, type == "A_B")@assays$RNA@data[genes[j], ]
        ra <- subset(sub, type == "I_B")@assays$RNA@data[genes[j], ]
        hc <- subset(sub, type == "T_B")@assays$RNA@data[genes[j], ]
        hh <- subset(sub, type == "N_B")@assays$RNA@data[genes[j], ]

        pvalue[j, (3 * (i - 1) + 1)] <- wilcox.test(oa, ra)$p.value
        pvalue[j, (3 * (i - 1) + 2)] <- wilcox.test(oa, hc)$p.value
        pvalue[j, (3 * (i - 1) + 3)] <- wilcox.test(oa, hh)$p.value

        fc[j, (3 * (i - 1) + 1)] <- log(mean(oa) / mean(ra), 2)
        fc[j, (3 * (i - 1) + 2)] <- log(mean(oa) / mean(hc), 2)
        fc[j, (3 * (i - 1) + 3)] <- log(mean(oa) / mean(hh), 2)
    }
}
for (i in 1:nrow(fc)) {
    for (j in 1:ncol(fc)) {
        if (is.nan(fc[i, j]) | is.infinite(fc[i, j])) {
            fc[i, j] <- NA
        }
        if (is.nan(pvalue[i, j])) {
            pvalue[i, j] <- NA
        }
    }
}
rownames(pvalue) <- genes
rownames(fc) <- genes
colnames(pvalue) <- paste(rep(kk, each = 3), rep(1:3, 6), sep = "_")
colnames(fc) <- paste(rep(kk, each = 3), rep(1:3, 6), sep = "_")
mt <- melt(pvalue)
mtt <- melt(fc)
mt$fc <- mtt$value
mt$pvalue <- ifelse(mt$value < 0.05, -log(mt$value, 10), NA)
mt$pvalue <- ifelse(mt$pvalue > 100, 100, mt$pvalue)
mt$fc <- ifelse(mt$fc > 1.6, 1.6, ifelse(mt$fc < -1.6, -1.6, mt$fc))
mt$Var2 <- factor(mt$Var2, levels = levels(mt$Var2)[c(1:9, 16:18, 10:12, 13:15)]) # nolint
mt$Var1 <- factor(mt$Var1, levels = rev(c(
    "JAK2", "FKBP1A", "IL12A", "PPP3CA", "JAK1", "JAK3",
    "FKBP5", "IL1A", "IL1B", "TNF", "IL18", "IL6", "IL12B"
))) # nolint

color <- c(muted("blue"), "white", muted("red"))
col <- colorRampPalette(color)(100)
col1 <- c(col[1:50], rep("#FDFDFD", 5), col[51:100])

p <- ggplot(mt, aes(x = Var2, y = Var1)) +
    geom_point(aes(color = fc, size = pvalue)) +
    theme_bw() +
    scale_color_gradientn(colors = col1) + # nolint
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
p
ggsave("Figure8/Figure_8D_exp.pdf", p, width = 5.5, height = 3.5)

matr <- matrix(NA, nrow = length(unique(data_exp$anno1)), ncol = length(genes)) # nolint

for (j in 1:length(unique(data_exp$anno1))) {
    sub <- subset(data_exp, anno1 == unique(data_exp$anno1)[j])
    for (i in 1:length(genes)) {
        matr[j, i] <- sum(sub@assays$RNA@data[genes[i], ] > 0) / ncol(sub)
    }
}
rownames(matr) <- unique(data_exp$anno1)
colnames(matr) <- genes
matr <- t(matr)
matr <- matr[, c(4, 1, 2, 5, 7, 6, 3, 8)]
pdf("phea.pdf", width = 3, height = 3)
pheatmap(matr,
    cluster_rows = F, cluster_cols = F,
    color = colorRampPalette(brewer.pal(
        n = 9, name = "Reds"
    ))(100), border_color = "black"
)
dev.off()
write.table(matr, "Fig2/fig2a_table.txt", quote = F, sep = "\t")

kk <- as.character(unique(data_exp$anno1))
pvalue <- matrix(NA, ncol = 2 * length(kk), nrow = length(genes)) # nolint
fc <- matrix(NA, ncol = 2 * length(kk), nrow = length(genes)) # nolint

for (i in 1:length(kk)) {
    print(i)
    sub <- subset(data_exp, anno1 == kk[i])
    for (j in 1:length(genes)) {
        print(j)
        oa <- subset(sub, group == "CAMR")@assays$RNA@data[genes[j], ]
        ra <- subset(sub, group == "SAF")@assays$RNA@data[genes[j], ]
        hc <- subset(sub, group == "HC")@assays$RNA@data[genes[j], ]

        pvalue[j, (2 * (i - 1) + 1)] <- wilcox.test(oa, ra)$p.value
        pvalue[j, (2 * (i - 1) + 2)] <- wilcox.test(oa, hc)$p.value

        fc[j, (2 * (i - 1) + 1)] <- log(mean(oa) / mean(ra), 2)
        fc[j, (2 * (i - 1) + 2)] <- log(mean(oa) / mean(hc), 2)
    }
}
for (i in 1:nrow(fc)) {
    for (j in 1:ncol(fc)) {
        if (is.nan(fc[i, j]) | is.infinite(fc[i, j])) {
            fc[i, j] <- NA
        }
        if (is.nan(pvalue[i, j])) {
            pvalue[i, j] <- NA
        }
    }
}
rownames(pvalue) <- genes
rownames(fc) <- genes
colnames(pvalue) <- paste(rep(kk, each = 2), rep(1:2, 5), sep = "_")
colnames(fc) <- paste(rep(kk, each = 2), rep(1:2, 5), sep = "_")
mt <- melt(pvalue)
mtt <- melt(fc)
mt$fc <- mtt$value
mt$pvalue <- ifelse(mt$value < 0.05, -log(mt$value, 10), NA)
mt$pvalue <- ifelse(mt$pvalue > 100, 100, mt$pvalue)
mt$fc <- ifelse(mt$fc > 2, 2, ifelse(mt$fc < -2, -2, mt$fc))
mt$Var2 <- factor(mt$Var2, levels = levels(mt$Var2)[c(5, 6, 1, 2, 3, 4, 7, 8, 9, 10)]) # nolint
mt$Var1 <- factor(mt$Var1, levels = rev(c(
    "JAK2", "FKBP1A", "IL12A", "PPP3CA", "JAK1",
    "JAK3", "FKBP5", "IL1A", "IL1B", "TNF", "IL18", "IL6", "IL12B"
)))
color <- c(muted("blue"), "white", muted("red"))
col <- colorRampPalette(color)(100)
col1 <- c(col[1:50], rep("#FDFDFD", 5), col[51:100])

p <- ggplot(mt, aes(x = Var2, y = Var1)) +
    geom_point(aes(color = fc, size = pvalue)) +
    theme_bw() +
    scale_color_gradientn(colors = col1) + # nolint
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
p
ggsave("Figure8/Figure_8D_exp1.pdf", p, width = 4, height = 3.5)
##################################################################


##################################################################
## Figure 8D
# ---------- Part 1: anno × (A_B vs I_B/T_B/N_B) ----------
data_exp <- readRDS("input_file/pan_blood_sub13Gene.Rds")
genes <- c(
    "IL6", "TNF", "IL1A", "IL1B", "IL12A", "IL12B", "IL18",
    "JAK1", "JAK2", "JAK3", "FKBP1A", "FKBP5", "PPP3CA"
)

# 固定 anno 顺序：优先用 factor 水平，否则按字母序
kk <- if (is.factor(data_exp$anno)) levels(data_exp$anno) else sort(unique(data_exp$anno)) # nolint
comp_names <- c("A_vs_I", "A_vs_T", "A_vs_N")
n_comp <- length(comp_names)

pvalue <- matrix(NA_real_,
    nrow = length(genes), ncol = n_comp * length(kk),
    dimnames = list(genes, as.vector(sapply(kk, function(a) paste(a, comp_names, sep = "_")))) # nolint
)
fc <- matrix(NA_real_,
    nrow = length(genes), ncol = n_comp * length(kk),
    dimnames = list(genes, colnames(pvalue))
)

safe_wilcox_p <- function(x, y) {
    if (length(x) < 2 || length(y) < 2) {
        return(NA_real_)
    }
    out <- tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
    as.numeric(out)
}
fc_from_log1p <- function(x, y, eps = 1e-8) {
    mx <- mean(expm1(x))
    my <- mean(expm1(y))
    if (!is.finite(mx) || !is.finite(my)) {
        return(NA_real_)
    }
    log2((mx + eps) / (my + eps))
}

for (i in seq_along(kk)) {
    sub <- subset(data_exp, anno == kk[i])
    if (ncol(sub) == 0) next
    M <- GetAssayData(sub, assay = "RNA", slot = "data")
    idxA <- which(sub$type == "A_B")
    idxI <- which(sub$type == "I_B")
    idxT <- which(sub$type == "T_B")
    idxN <- which(sub$type == "N_B")

    for (j in seq_along(genes)) {
        g <- genes[j]
        if (!g %in% rownames(M)) {
            pvalue[j, (n_comp * (i - 1) + 1):(n_comp * i)] <- NA
            fc[j, (n_comp * (i - 1) + 1):(n_comp * i)] <- NA
            next
        }
        xA <- as.numeric(M[g, idxA, drop = FALSE])
        xI <- as.numeric(M[g, idxI, drop = FALSE])
        xT <- as.numeric(M[g, idxT, drop = FALSE])
        xN <- as.numeric(M[g, idxN, drop = FALSE])

        pvalue[j, n_comp * (i - 1) + 1] <- safe_wilcox_p(xA, xI)
        pvalue[j, n_comp * (i - 1) + 2] <- safe_wilcox_p(xA, xT)
        pvalue[j, n_comp * (i - 1) + 3] <- safe_wilcox_p(xA, xN)

        fc[j, n_comp * (i - 1) + 1] <- fc_from_log1p(xA, xI)
        fc[j, n_comp * (i - 1) + 2] <- fc_from_log1p(xA, xT)
        fc[j, n_comp * (i - 1) + 3] <- fc_from_log1p(xA, xN)
    }
}

# 熵值/NaN 清理
fc[!is.finite(fc)] <- NA
pvalue[!is.finite(pvalue)] <- NA

mt_p <- melt(pvalue, varnames = c("gene", "anno_comp"), value.name = "p")
mt_fc <- melt(fc, varnames = c("gene", "anno_comp"), value.name = "fc")
mt <- left_join(mt_p, mt_fc, by = c("gene", "anno_comp")) |>
    tidyr::separate(
        col = anno_comp, into = c("anno", "comp"),
        sep = "_", extra = "merge", fill = "right", remove = FALSE
    )

# 变换 & 顺序
mt$neglog10p <- ifelse(is.na(mt$p), NA_real_, pmin(-log10(mt$p), 100))
mt$fc_clip <- pmax(pmin(mt$fc, 1.6), -1.6)
mt$anno <- factor(mt$anno, levels = c("Myeloid", "T Cell", "NK cell", "B Cell", "Plasma", "Plasmablast")) # nolint
mt$comp <- factor(mt$comp, levels = comp_names)
cell_order <- c("Myeloid", "T Cell", "NK", "B Cell", "Plasma", "Plasmablast")
comp_names <- c("A_vs_I", "A_vs_T", "A_vs_N")
lvl <- as.vector(t(outer(cell_order, comp_names, paste, sep = "_")))
mt$anno_comp <- factor(mt$anno_comp, levels = lvl)
mt$gene <- factor(mt$gene, levels = rev(c(
    "JAK2", "FKBP1A", "IL12A", "PPP3CA", "JAK1", "JAK3",
    "FKBP5", "IL1A", "IL1B", "TNF", "IL18", "IL6", "IL12B"
)))

# 颜色
color <- c(scales::muted("blue"), "white", scales::muted("red"))
col <- colorRampPalette(color)(100)
col1 <- c(col[1:50], rep("#FDFDFD", 5), col[51:100])

p1 <- ggplot(mt, aes(x = anno_comp, y = gene)) +
    geom_point(aes(color = fc_clip, size = neglog10p)) +
    theme_bw() +
    scale_color_gradientn(colors = col1) +
    scale_size_continuous(range = c(0.2, 3), limits = c(0, 100)) +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

dir.create("Figure8", showWarnings = FALSE, recursive = TRUE)
ggsave("Figure8/Figure_8D_exp.pdf", p1, width = 5.5, height = 3.5)

# ---------- Part 2: 阳性比例热图（按 anno1） ----------
lev1 <- if (is.factor(data$anno1)) levels(data$anno1) else unique(data$anno1)
matr <- matrix(NA_real_,
    nrow = length(lev1), ncol = length(genes),
    dimnames = list(lev1, genes)
)
for (j in seq_along(lev1)) {
    sub <- subset(data, anno1 == lev1[j])
    if (ncol(sub) == 0) next
    M <- GetAssayData(sub, assay = "RNA", slot = "data")
    matr[j, ] <- vapply(genes, function(g) {
        if (!g %in% rownames(M)) {
            return(NA_real_)
        }
        sum(M[g, ] > 0) / ncol(sub)
    }, numeric(1))
}
matr <- t(matr)

# 如果你确实需要固定列顺序（且列数=8），保留你的顺序；否则按现有列名画
if (ncol(matr) >= 8) {
    ord <- c(4, 1, 2, 5, 7, 6, 3, 8)
    ord <- ord[ord <= ncol(matr)]
    matr <- matr[, ord, drop = FALSE]
}

# ---------- Part 3: anno1 × (CAMR vs SAF/HC) ----------
kk2 <- if (is.factor(data$anno1)) levels(data$anno1) else sort(unique(data$anno1)) # nolint
comp2 <- c("CAMR_vs_SAF", "CAMR_vs_HC")
n_comp2 <- length(comp2)

pvalue2 <- matrix(NA_real_,
    nrow = length(genes), ncol = n_comp2 * length(kk2),
    dimnames = list(genes, as.vector(sapply(kk2, function(a) paste(a, comp2, sep = "_")))) # nolint
)
fc2 <- matrix(NA_real_,
    nrow = length(genes), ncol = n_comp2 * length(kk2),
    dimnames = list(genes, colnames(pvalue2))
)

for (i in seq_along(kk2)) {
    sub <- subset(data, anno1 == kk2[i])
    if (ncol(sub) == 0) next
    M <- GetAssayData(sub, assay = "RNA", slot = "data")
    idxC <- which(sub$group == "CAMR")
    idxS <- which(sub$group == "SAF")
    idxH <- which(sub$group == "HC")

    for (j in seq_along(genes)) {
        g <- genes[j]
        if (!g %in% rownames(M)) {
            pvalue2[j, (n_comp2 * (i - 1) + 1):(n_comp2 * i)] <- NA
            fc2[j, (n_comp2 * (i - 1) + 1):(n_comp2 * i)] <- NA
            next
        }
        xC <- as.numeric(M[g, idxC, drop = FALSE])
        xS <- as.numeric(M[g, idxS, drop = FALSE])
        xH <- as.numeric(M[g, idxH, drop = FALSE])

        pvalue2[j, n_comp2 * (i - 1) + 1] <- safe_wilcox_p(xC, xS)
        pvalue2[j, n_comp2 * (i - 1) + 2] <- safe_wilcox_p(xC, xH)

        fc2[j, n_comp2 * (i - 1) + 1] <- fc_from_log1p(xC, xS)
        fc2[j, n_comp2 * (i - 1) + 2] <- fc_from_log1p(xC, xH)
    }
}

fc2[!is.finite(fc2)] <- NA
pvalue2[!is.finite(pvalue2)] <- NA

mt2_p <- melt(pvalue2, varnames = c("gene", "anno_comp"), value.name = "p")
mt2_fc <- melt(fc2, varnames = c("gene", "anno_comp"), value.name = "fc")
mt2 <- left_join(mt2_p, mt2_fc, by = c("gene", "anno_comp")) |>
    tidyr::separate(
        col = anno_comp, into = c("anno1", "comp"),
        sep = "_", extra = "merge", fill = "right", remove = FALSE
    )
mt2$neglog10p <- ifelse(is.na(mt2$p), NA_real_, pmin(-log10(mt2$p), 100))
mt2$fc_clip <- pmax(pmin(mt2$fc, 2), -2)
mt2$anno1 <- factor(mt2$anno1, levels = kk2)
mt2$comp <- factor(mt2$comp, levels = comp2)
cell_order <- c("Myeloid", "T cell", "NK cell", "B cell", "Plasma")
comp_names <- c("CAMR_vs_SAF", "CAMR_vs_HC")
lvl <- as.vector(t(outer(cell_order, comp_names, paste, sep = "_")))
mt2$anno_comp <- factor(mt2$anno_comp, levels = lvl)
mt2$gene <- factor(mt2$gene, levels = rev(c(
    "JAK2", "FKBP1A", "IL12A", "PPP3CA", "JAK1",
    "JAK3", "FKBP5", "IL1A", "IL1B", "TNF", "IL18", "IL6", "IL12B"
)))

color <- c(scales::muted("blue"), "white", scales::muted("red"))
col <- colorRampPalette(color)(100)
col1 <- c(col[1:50], rep("#FDFDFD", 5), col[51:100])

p2 <- ggplot(mt2, aes(x = anno_comp, y = gene)) +
    geom_point(aes(color = fc_clip, size = neglog10p)) +
    theme_bw() +
    scale_color_gradientn(colors = col1) +
    scale_size_continuous(range = c(0.2, 3), limits = c(0, 100)) +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

dir.create("Figure8", showWarnings = FALSE, recursive = TRUE)
ggsave("Figure8/Figure_8D_exp1.pdf", p2, width = 4, height = 3.5)
##################################################################
