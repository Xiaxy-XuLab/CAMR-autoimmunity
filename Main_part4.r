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
    "ComplexHeatmap", "circlize", "stringr", "tibble", "SeuratDisk"
)
##################################################################


##################################################################
## loading object
setwd("/work/xiaxy/work_2024/wangll_blood_final1/output2/")
data <- readRDS("input_file/input.rds")
nkt <- readRDS("input_file/nkt_sub.rds")
t_data <- readRDS("input_file/tissue_input.rds")
nkt_tdata <- readRDS("tissue_result/tissue_nkt.rds")
nk_pan <- readRDS("../blood_mao/nk_cell/nk_blood_mao_input.rds")
t_pan <- readRDS("../blood_mao/t_cell/t_blood_mao_input.rds")
nkt_cd16 <- readRDS("input_file/cd16nk_new.rds")
load("input_file/color.RData")
##################################################################


##################################################################
## Figure 5B HLA receptors analysis
# 1) 基因集合
receptor <- list(
    "HLA-dependent\ninhibitory receptors"   = c("KIR2DL1", "KIR2DL3", "KIR3DL1", "KIR3DL2", "LILRB1", "LAG3"), # nolint
    "HLA-independent\ninhibitory receptors" = c("PDCD1", "SIGLEC7", "CD300A", "CD96", "IL1RAPL1", "TIGIT", "HAVCR2"), # nolint
    "HLA-dependent\nactivating receptors"   = c("KIR2DL4", "CD160", "KLRC2"),
    "HLA-independent\nactivating receptors" = c("NCR3", "NCR1", "KLRK1", "CRTAM", "FCGR3A") # nolint
)

# 2) 取子集
nkt_data <- subset(nkt, anno2 %in% c("tn16_ZEB2+CD16+ NK", "tn17_CD56+ NK"))

# 3) 一个可复用的“按命名列表批量加模块分数”的函数
#    注意：AddModuleScore 会生成 name + 1 的列名，这里自动改回原名（含换行保持）
fix_module_score_names <- function(obj, gene_sets) {
    md <- obj@meta.data
    for (nm in names(gene_sets)) {
        # AddModuleScore 通常生成 make.names(nm) 再加 "1"
        cand <- paste0(make.names(nm), "1")
        # 也兼容极少数情况下直接是 "nm1" 的可能
        alt <- paste0(nm, "1")
        if (cand %in% colnames(md)) {
            colnames(md)[match(cand, colnames(md))] <- nm
        } else if (alt %in% colnames(md)) {
            colnames(md)[match(alt, colnames(md))] <- nm
        } else if (nm %in% colnames(md)) {
            # 已经是目标名就跳过
            next
        } else {
            warning(
                "找不到模块分数列：", nm,
                " （既没有 ", cand, " 也没有 ", alt, "）"
            )
        }
    }
    obj@meta.data <- md
    obj
}

# 在 add_module_scores/ AddModuleScore 之后调用一次
nkt_data <- fix_module_score_names(nkt_data, receptor)

# 4) 此时你就可以用原来的 select 了
matr_t <- nkt_data@meta.data %>%
    tibble::as_tibble(rownames = "cell") %>%
    dplyr::select(anno2, dplyr::all_of(names(receptor))) %>%
    tidyr::pivot_longer(
        cols = dplyr::all_of(names(receptor)),
        names_to = "variable",
        values_to = "value"
    )

# 5) 画图函数（保持你的样式完全一致）
draw <- function(xx, yy) {
    ggplot(xx, aes(x = anno2, y = value)) +
        geom_violin(aes(fill = anno2), scale = "width", color = "white") +
        scale_fill_manual(values = c("#EBD57C99", "#E93A3B99")) +
        new_scale_fill() +
        scale_y_continuous(limits = yy) +
        geom_boxplot(aes(fill = anno2),
            color = "black", width = 0.1,
            outlier.alpha = 0, outlier.size = 0, linewidth = 0.4
        ) +
        scale_fill_manual(values = c("#EBD57C", "#E93A3B")) +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(color = "white", size = 8),
            axis.title = element_blank(),
            legend.position = "none",
            axis.line = element_line(linewidth = 0.3),
            axis.ticks = element_line(linewidth = 0.3)
        )
}

# 6) 为每个面板配置 y 轴范围（与原代码一致，可复用）
ylims_map <- list(
    "HLA-dependent\ninhibitory receptors"   = c(-0.4, 0.8),
    "HLA-independent\ninhibitory receptors" = c(-0.4, 0.6),
    "HLA-dependent\nactivating receptors"   = c(-1, 1.2),
    "HLA-independent\nactivating receptors" = c(-1, 1.2)
)

# 7) 分面生成四个子图，再用 patchwork 横向拼接（与原来 1 行一致）
vars <- intersect(unique(matr_t$variable), names(ylims_map))

plots <- map(vars, ~ {
    df_i <- filter(matr_t, variable == .x)
    draw(df_i, ylims_map[[.x]])
})

# 0/1/多面板的稳健拼接方案
if (length(plots) == 0) {
    stop("No panels to plot: 请检查 matr_t$variable 与 ylims_map / receptor 是否匹配。")
} else if (length(plots) == 1) {
    final <- plots[[1]] + plot_layout(nrow = 1, guides = "collect")
} else {
    # 用 wrap_plots 更简洁；也可用 base::Reduce(`+`, plots)
    final <- wrap_plots(plots, nrow = 1, guides = "collect")
}

ggsave("Figure5/Figure_5B_NK_hla.pdf", final, width = 7, height = 2)
##################################################################


##################################################################
## Figure 5D trajectory analysis for CD16+ NK
## monocle3 gene analysis
# 1) ident 设置与差异分析（保持原调用）
nkt_cd16@active.ident <- factor(nkt_cd16$group)
marker <- FindAllMarkers(nkt_cd16, only.pos = TRUE, min.pct = 0.25)

# 2) 一次性抓取所需数据（含 ptime 与目标基因）
vars_use <- c(
    "ptime", "JUN", "JUNB", "JUND", "FOS", "FOSB", "CXCR4",
    "FKBP5", "PPP3CA", "DUSP1", "ZEB2"
)
set.seed(100)

mat <- FetchData(nkt_cd16, vars = vars_use) %>%
    as_tibble() %>%
    transmute(
        pse = 1 - ptime,
        jun = JUN,
        junb = JUNB,
        jund = JUND,
        fos = FOS,
        fosb = FOSB,
        cxcr4 = CXCR4,
        fkbp5 = FKBP5,
        ppp3ca = PPP3CA,
        dusp1 = DUSP1,
        zeb2 = ZEB2
    ) %>%
    slice_sample(n = 1000, replace = FALSE)

# 3) 通用的 loess 平滑图（保持你的样式/参数/配色）
smooth_panel <- function(df, yvar) {
    ggplot(df, aes(x = pse, y = .data[[yvar]])) +
        geom_smooth(method = "loess", se = TRUE, level = 0.5, color = "#4C7FB2") + # nolint
        theme_classic() +
        labs(x = "", y = "") +
        theme(
            axis.text         = element_blank(),
            legend.text       = element_text(color = "black", size = 10),
            legend.title      = element_blank(),
            axis.line         = element_line(linewidth = 0.6),
            axis.ticks        = element_line(linewidth = 0.6),
            strip.text        = element_text(color = "black", size = 13),
            strip.background  = element_rect(fill = NA, color = "white")
        )
}

# 4) 三个单基因面板（与你原来的 p1/p2/p4 对应）
p1 <- smooth_panel(mat, "cxcr4")
p2 <- smooth_panel(mat, "ppp3ca")
p4 <- smooth_panel(mat, "zeb2")

# 5) 多基因叠加面板（JUN/JUNB/JUND/FOS/FOSB），保持原样式
mt <- mat %>%
    select(pse, jun, junb, jund, fos, fosb) %>%
    pivot_longer(-pse, names_to = "variable", values_to = "value")

p5 <- ggplot(mt, aes(x = pse, y = value, color = variable)) +
    geom_smooth(method = "loess", se = TRUE, level = 0.5) +
    scale_color_npg() +
    theme_classic() +
    labs(x = "", y = "") +
    theme(
        axis.text         = element_blank(),
        legend.text       = element_text(color = "black", size = 10),
        legend.title      = element_blank(),
        axis.line         = element_line(linewidth = 0.6),
        axis.ticks        = element_line(linewidth = 0.6),
        strip.text        = element_text(color = "black", size = 13),
        strip.background  = element_rect(fill = NA, color = "white")
    )

# 6) 拼图与保存（布局与尺寸保持）
final <- p1 + p2 + p4 + p5 + plot_layout(nrow = 2)
ggsave("Figure5/Figure_5D_marker_nk_pse.pdf", final, width = 5.6, height = 5)
##################################################################


##################################################################
## Figure 5C trajectory analysis for CD16+ NK
# prepare for pyVIA
Seu2Ann <- function(seu_obj = a, filename) {
    SaveH5Seurat(seu_obj, filename = paste0(filename, ".h5Seurat"))
    Convert(paste0(filename, ".h5Seurat"), dest = "h5ad")
}

Seu2Ann(seu_obj = nkt_cd16, filename = "Anndata")
write.csv(data.frame(row.names = 1:ncol(nkt_cd16), cell_id = rownames(nkt_cd16@meta.data), group_id = nkt_cd16$group), "Anndata_id.csv", quote = F) # nolint

# slingshot
# —— 1) 分组因子化（如需改顺序，只改这里）——
grp_levels <- c("HC", "SAF", "CAMR")
nkt_cd16$group <- factor(nkt_cd16$group, levels = grp_levels)

# （可选）保留原 Idents 到元数据，便于追溯
nkt_cd16$idents <- as.character(Idents(nkt_cd16))

# —— 2) 转换为 SCE 并运行 slingshot ——
a.sce <- as.SingleCellExperiment(nkt_cd16)
# 确保 UMAP 已在 nkt_cd16 中；否则需先 RunUMAP，再 as.SingleCellExperiment
sim <- slingshot(a.sce, clusterLabels = "group", reducedDim = "UMAP")

# —— 3) 颜色映射：确保与分组严格对齐 ——
# 你的原色板
base_cols <- c("#749D51", "#4F89BB", "#A65E34")
# 实际分组（去掉 NA），按因子水平排序
present_lvls <- levels(droplevels(colData(sim)$group))
# 若 present_lvls 少于既定 levels，会自动截取对应数量的颜色
col_map <- setNames(
    base_cols[seq_len(min(length(base_cols), length(present_lvls)))],
    present_lvls
)

# —— 4) 绘图（保持你的样式：base plot + lines）——
pdf("Figure5/Figure_5C_trace.pdf", height = 15, width = 15)
on.exit(dev.off(), add = TRUE)

umap_mat <- reducedDim(sim, "UMAP") # 等价于 reducedDims(sim)$UMAP
grp_vec <- as.character(colData(sim)$group)

# 点图（与原参数一致：pch = 16，asp = 1.4）
plot(
    umap_mat,
    col = col_map[grp_vec],
    pch = 16,
    asp = 1.4,
    xlab = "UMAP_1",
    ylab = "UMAP_2"
)
# 轨迹（与原样式一致）
lines(SlingshotDataSet(sim), lwd = 8, type = "lineages", col = "black")
##################################################################


##################################################################
## Figure 5F-5G Fraction & expression
# Figure 5G
DefaultAssay(nkt_cd16) <- "RNA"

# 2) 在 meta.data 上进行变换（不走“在管道里改对象”的花活，简单稳健）
nkt_cd16@meta.data <- nkt_cd16@meta.data %>%
    mutate(
        ptime = as.numeric(ptime),
        ptime = 1 - ptime,
        ss = case_when(
            ptime < 0.3908367 ~ "SI",
            ptime > 0.5604390 ~ "SIII",
            TRUE ~ "SII"
        ),
        ss = factor(ss, levels = c("SI", "SII", "SIII"))
    )

# 3) 同步 Idents
Idents(nkt_cd16) <- nkt_cd16$ss

# 4) 差异分析保持不变
marker <- FindAllMarkers(nkt_cd16, only.pos = TRUE, min.pct = 0.25)

# 5) 小提琴图（等价改写：Seurat::VlnPlot，样式/配色/参数保持）
pdf("Figure5/Figure_5G_ap1_violin.pdf", width = 10, height = 10)
p_vln <- VlnPlot(
    nkt_cd16,
    features = c("JUN", "JUNB", "JUND", "FOS"),
    group.by = "ss",
    split.by = "ss", # 你原来就写了 split.by="ss"
    cols = c("#E2AF6F", "#8E84A7", "#AB6467"),
    pt.size = 0
) +
    theme_classic() +
    theme(
        legend.position   = "none",
        axis.title.x      = element_blank(),
        strip.background  = element_blank(),
        strip.text        = element_blank(),
        axis.title        = element_blank(),
        axis.text         = element_blank(),
        axis.line         = element_line(linewidth = 0.5, color = "black"),
        axis.ticks        = element_line(linewidth = 0.5, color = "black")
    )
print(p_vln)
dev.off()

## -------------------- Figure 5F：分组比例曲线（均值 ± SE） --------------------

# 思路：按 ptime 排序后，固定比例切成 SI/SII/SIII（与你原逻辑一致：40%/20%/剩余）
first <- 0.4
second <- 0.2

mat <- nkt_cd16@meta.data %>%
    as_tibble(rownames = "cell") %>%
    arrange(ptime) %>%
    mutate(
        idx = row_number(),
        n = n(),
        ss = case_when(
            idx <= floor(first * n) ~ "SI",
            idx <= floor((first + second) * n) ~ "SII",
            TRUE ~ "SIII"
        ),
        ss = factor(ss, levels = c("SI", "SII", "SIII"))
    )

# 如果你的分组在 meta.data 里已有（推荐），直接用：
# mat$group 已存在就不用下面这行；否则请确保有 group 列
# mat <- mat %>% mutate(group = nkt_cd16$group[match(cell, rownames(nkt_cd16@meta.data))]) # nolint

# 4) 计算：每个样本（orig.ident）内各阶段的比例，再在组内求均值与 SE
mgg <- mat %>%
    count(orig.ident, group, ss, name = "cells") %>%
    group_by(orig.ident) %>%
    mutate(prop = cells / sum(cells)) %>%
    ungroup() %>%
    group_by(ss, group) %>%
    summarise(
        Proportion = mean(prop, na.rm = TRUE),
        SE         = sd(prop, na.rm = TRUE) / sqrt(dplyr::n()),
        .groups    = "drop"
    ) %>%
    mutate(
        Type     = factor(group, levels = c("HC", "SAF", "CAMR")),
        Category = ss
    ) %>%
    select(Category, Type, Proportion, SE) %>%
    arrange(Category, Type)

# 5) 画图（保持你的风格）
pdf("Figure5/Figure_5F_nk_se.pdf", width = 3.5, height = 2)
ggplot(mgg, aes(x = Category, y = Proportion, color = Type, group = Type)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = mycolor$group) +
    geom_errorbar(aes(ymin = Proportion - SE, ymax = Proportion + SE), width = 0.1) + # nolint
    labs(y = "Proportion", color = "") +
    theme_classic() +
    scale_x_discrete(expand = c(0.05, 0)) +
    theme(
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank(),
        axis.line    = element_line(linewidth = 0.5),
        axis.ticks   = element_line(linewidth = 0.5)
    )
dev.off()
##################################################################


##################################################################
## Figure 5A Heatmap
# 1) 取子集 & 计算 marker
nkt_sub <- subset(nkt, anno3 %in% c("tn16", "tn17"))

marker <- FindAllMarkers(nkt_sub, only.pos = TRUE, min.pct = 0.25) %>%
    as_tibble() %>%
    filter(p_val_adj < 0.05)

# 2) 下采样：可复现
set.seed(2024)
nkt_sub <- subset(nkt_sub, downsample = 1000)

# 3) 目标标签（与你一致），并确保基因存在于对象中
label <- c(
    "FGFBP2", "CX3CR1", "FCGR3A", "GZMH", "GZMB", "CCL4", "ZEB2",
    "PRF1", "GZMA", "GZMM", "CCL5", "GZMK", "SELL", "XCL1", "CD74", "NCAM1" # nolint
)

# 4) 供热图使用的基因向量（去重，并确保在对象里）
genes_use <- marker %>%
    pull(gene) %>%
    unique() %>%
    intersect(rownames(nkt_sub))

# 5) 计算需要标注的行索引（与非聚类行序一致）
#    注意：cluster_row = FALSE → 行顺序就是 genes_use 的顺序
label_idx <- match(label, genes_use) %>% na.omit()

# 6) 出图（样式、配色、参数与你完全一致）
ht <- dittoHeatmap(
    nkt_sub,
    genes = genes_use, # 你准备好的基因向量
    cluster_row = FALSE,
    scaled.to.max = TRUE,
    complex = TRUE, # 关键：返回 Heatmap/HeatmapList
    raster_quality = 6,
    order.by = "anno2",
    heatmap.colors.max.scaled = colorRampPalette(c("#FCFDBFFF", rev(magma(10))))(50), # nolint
    show_colnames = FALSE
)

# 只用 ComplexHeatmap 的注释
label_idx <- match(label, genes_use)
label_idx <- label_idx[!is.na(label_idx)]
if (length(label_idx) > 0) {
    ht <- ht + ComplexHeatmap::rowAnnotation(
        label = ComplexHeatmap::anno_mark(
            at = label_idx,
            labels = label[!is.na(match(label, genes_use))]
        )
    )
}

pdf("Figure5/Figure_5A_nk_deg_heatmap.pdf", width = 5.5, height = 4)
ComplexHeatmap::draw(ht) # 用 ComplexHeatmap 的 draw 来渲染
dev.off()
##################################################################


##################################################################
## Figure 5H SIII vs. SI Deg
# 1) 读入与分段
nkt_sub <- readRDS("input_file/cd16nk_new.rds")
DefaultAssay(nkt_sub) <- "RNA"
nkt_sub$ptime <- 1 - nkt_sub$ptime
nkt_sub$ss <- ifelse(
    nkt_sub$ptime < 0.3908367, "SI",
    ifelse(nkt_sub$ptime > 0.5604390, "SIII", "SII")
)

# 2) 读入 SCENIC 结果 & 名称清洗
sce <- read.table("input_file/cd16_nk_SCENIC.txt", header = TRUE, row.names = 1, sep = "\t") # nolint
colnames(sce) <- gsub("[.]", "-", colnames(sce))
rownames(sce) <- gsub("_extended", "", sapply(strsplit(rownames(sce), " ", fixed = TRUE), "[", 1)) # nolint

# 3) 细胞对齐：两边都限制到共同细胞，并按同一顺序排列（很重要）
common_cells <- intersect(colnames(nkt_sub), colnames(sce))
nkt_sub <- subset(nkt_sub, cells = common_cells)
sce <- sce[, common_cells, drop = FALSE]

# 4) 计算“激活分数-ptime”相关（与 cor.test 默认 Pearson 一致）
#    向量化：对每个 regulon（行）与 ptime 做一次 cor()
ptime_vec <- nkt_sub$ptime %>% as.numeric()
r_activation <- cor(t(as.matrix(sce)), ptime_vec, use = "pairwise.complete.obs", method = "pearson") # nolint
# r_activation: 命名数值向量（names = rownames(sce)）

# 5) 计算“表达-ptime”相关
#    先取表达矩阵（data 槽，log-normalized），对齐基因 & 细胞顺序
expr <- GetAssayData(nkt_sub, assay = "RNA", slot = "data")
genes_in_expr <- intersect(rownames(sce), rownames(expr))
# 对存在于表达矩阵中的基因计算相关
r_expr_part <- cor(t(as.matrix(expr[genes_in_expr, common_cells, drop = FALSE])), # nolint
    ptime_vec,
    use = "pairwise.complete.obs", method = "pearson"
)
# 把没有表达数据的基因补 NA，保证后续合并不丢行
r_expression <- rep(NA_real_, nrow(sce))
names(r_expression) <- rownames(sce)
r_expression[genes_in_expr] <- r_expr_part

# 6) 差异
marker <- FindMarkers(
    nkt_sub,
    ident.1 = "SIII", ident.2 = "SI",
    group.by = "ss", min.pct = 0, logfc.threshold = 0
)
marker$gene <- rownames(marker)
marker <- subset(marker, gene %in% rownames(sce))

# 7) 组装表并做 log2FC 的缩放
mat <- tibble(
    gene = rownames(sce),
    expression = as.numeric(r_expression),
    activation = as.numeric(r_activation)
)

merge_tab <- mat %>%
    inner_join(marker, by = "gene")

ra <- max(merge_tab$avg_log2FC, na.rm = TRUE) / abs(min(merge_tab$avg_log2FC, na.rm = TRUE)) # nolint
merge_tab <- merge_tab %>%
    mutate(
        avg_log2FC = ifelse(avg_log2FC < 0, ra * avg_log2FC, avg_log2FC),
        p_val_plot = pmax(p_val, .Machine$double.xmin)
    ) # 避免 log10(0)

# 8) 作图
p <- ggplot(merge_tab, aes(x = activation, y = expression)) +
    geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    scale_color_gradient2(low = "#2A2C51", mid = "lightgrey", high = "#8D5A50", midpoint = 0) + # nolint
    stat_smooth(method = "lm", se = TRUE, level = 0.5, color = "#4C7FB2") +
    scale_size_continuous(range = c(2, 5)) +
    geom_point(aes(size = -log(p_val_plot, 10), color = avg_log2FC)) +
    theme_classic() +
    theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank()
    )

ggsave("Figure5/Figure_5H_SIII_SI_Deg.pdf", p, width = 3.8, height = 2.8)
##################################################################


##################################################################
## Figure 5I TF analysis
# 1) 读入与分段
nkt_sub <- readRDS("input_file/cd16nk_new.rds")
nkt_sub$ptime <- 1 - nkt_sub$ptime
nkt_sub$ss <- ifelse(nkt_sub$ptime < 0.3908367, "SI",
    ifelse(nkt_sub$ptime > 0.5604390, "SIII", "SII")
)

# 2) 读入 SCENIC 并清洗命名
sce <- read.table("input_file/cd16_nk_SCENIC.txt", header = TRUE, row.names = 1, sep = "\t") # nolint
colnames(sce) <- gsub("[.]", "-", colnames(sce))
rownames(sce) <- gsub("_extended", "", sapply(strsplit(rownames(sce), " ", fixed = TRUE), "[", 1)) # nolint

# 3) 细胞对齐 + 仅保留 SIII（顺序严格一致，避免潜在错位）
common_cells <- intersect(colnames(sce), colnames(nkt_sub))
nkt_sub <- subset(nkt_sub, cells = common_cells)
sce <- sce[, common_cells, drop = FALSE]

nkt_sub <- subset(nkt_sub, ss == "SIII")
sce <- sce[, colnames(nkt_sub), drop = FALSE]

# 4) 计算每个 regulon 的平均活化值并排序 + 排名
mat <- tibble(
    name  = rownames(sce),
    value = rowMeans(sce, na.rm = TRUE)
) %>%
    arrange(value) %>%
    mutate(rank = row_number())

# 5) 画图（保持你的主题与配色；改用 annotate 避免继承映射）
p <- ggplot(mat, aes(x = rank, y = value)) +
    annotate("rect",
        xmin = 30, xmax = 47, ymin = 0, ymax = 0.14,
        fill = "#CBAAAB", color = NA
    ) +
    geom_point(color = "#4C2570") +
    theme_bw() +
    theme(
        panel.grid   = element_blank(),
        axis.title   = element_blank()
    )

ggsave("Figure5/Figure_5I_sce_order.pdf", p, width = 2, height = 2)
##################################################################


##################################################################
## Figure 5K
# 统计并按 group 行归一化
mat <- nk_pan@meta.data %>%
  as_tibble() %>%
  count(group, anno, name = "n") %>%
  group_by(group) %>%
  mutate(value = n / sum(n)) %>%
  ungroup()

# 因子顺序
mat <- mat %>%
  mutate(
    group = factor(group, levels = c("N", "T", "I", "A")),
    anno  = factor(anno,  levels = rev(levels(factor(anno))))
  )

# 画图
p <- ggplot(mat, aes(x = value, y = anno, fill = group)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  theme_classic() +
  scale_fill_manual(values = pal_futurama(alpha = 0.8)(4)[c(3, 4, 1, 2)]) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    axis.line.y = element_blank()
  )
ggsave("Figure5/Figure_5J_pan_nk_percent.pdf", p, width = 4.6, height = 1.9)
##################################################################
