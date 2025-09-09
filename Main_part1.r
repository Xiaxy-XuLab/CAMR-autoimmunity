#!/usr/bin/env R
# coding utf-8

###################### Package Loading ###########################
pacman::p_load(
    "Seurat", "Nebulosa", "ggplot2", "future", "reshape2",
    "SingleCellExperiment", "dplyr", "tidyverse", "ggrepel",
    "patchwork", "msigdbr", "GSVA", "RColorBrewer", "ggpubr",
    "ROGUE", "viridis", "magrittr", "data.table", "tidyr",
    "R.utils", "grid", "cowplot", "tidyverse", "dittoSeq",
    "harmony", "scRepertoire", "ggsci", "pheatmap", "ggpie",
    "sscVis", "alakazam", "UpSetR", "CytoTRACE", "ggforce", "forcats",
    "tidydr", "ggplotify", "vegan", "batchelor", "SeuratWrappers"
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
load("input_file/color.RData")
##################################################################


##################################################################
## Figure 1C Major Cluster
pdf("Figure1/Figure_1c_umap_major.pdf", width = 30, height = 30)
DimPlot(data,
    label = F, group.by = "anno1",
    cols = mycolor$all_type,
    pt.size = 0.1, raster = FALSE, seed = 2
) +
    NoLegend() + NoAxes() +
    labs(title = "")
dev.off()

## Figure 1C Myeloid cells
del_cell <- names(which(mye@reductions$umap@cell.embeddings[, 1] > 7 & mye$anno3 == "m11")) # nolint
mye <- subset(mye, cells = setdiff(rownames(mye@meta.data), del_cell))
pdf("Figure1/Figure_1c_umap_myeloid.pdf", width = 25, height = 25)
DimPlot(mye,
    label = F, group.by = "anno2",
    cols = mycolor$mye_type,
    pt.size = 0.1, raster = FALSE
) +
    NoLegend() + NoAxes() +
    labs(title = "")
dev.off()

## Figure 1C B & Plasma cells
bp@meta.data <- cbind(bp@meta.data, bp@reductions$umap@cell.embeddings)
del_cell <- rownames(bp@meta.data)[intersect(which(bp$UMAP_1 < -5), which(bp$anno2 == "bp01.Pro-BC"))] # nolint
bp <- subset(bp, cells = setdiff(colnames(bp), del_cell))

pdf("Figure1/Figure_1c_umap_bp.pdf", width = 6, height = 6)
DimPlot(bp,
    label = F, group.by = "anno2",
    cols = mycolor$bp_type,
    pt.size = 0.1, raster = FALSE
) +
    NoLegend() + NoAxes() +
    labs(title = "")
dev.off()

## Figure 1C NK & T cells
pdf("Figure1/Figure_1c_umap_nkt.pdf", width = 25, height = 25)
DimPlot(nkt,
    label = F, group.by = "anno2",
    cols = mycolor$nkt_type,
    pt.size = 0.1, raster = FALSE
) +
    NoLegend() + NoAxes() +
    labs(title = "")
dev.off()
##################################################################


##################################################################
## Figure 1D Tissue Major Cluster
pdf("Figure1/Figure_1d_umap_tissue.pdf", width = 30, height = 30)
DimPlot(t_data,
    label = F, group.by = "anno4",
    cols = mycolor$tissue_type,
    pt.size = 0.1, raster = FALSE
) +
    NoLegend() + NoAxes() +
    labs(title = "")
dev.off()
##################################################################


##################################################################
### p11 Ro/e
do.tissueDist <- function(cellInfo.tb = cellInfo.tb, # meta.data
                          meta.cluster = cellInfo.tb$meta.cluster, # 纵坐标，可以是不同的分群  # nolint
                          loc = cellInfo.tb$loc, # 不同的分组，可以是肿瘤，癌旁，正常等
                          out.prefix, # 输出文件的名字
                          pdf.width = 3,
                          pdf.height = 5,
                          verbose = 0 # 如果等于1，则返回对象，可以对对象稍作处理重新画图
) {
    ## input data
    dir.create(dirname(out.prefix), F, T)

    cellInfo.tb <- data.table(cellInfo.tb)
    cellInfo.tb$meta.cluster <- as.character(meta.cluster)

    if (is.factor(loc)) {
        cellInfo.tb$loc <- loc
    } else {
        cellInfo.tb$loc <- as.factor(loc)
    }

    loc.avai.vec <- levels(cellInfo.tb[["loc"]])
    count.dist <- unclass(cellInfo.tb[, table(meta.cluster, loc)])[, loc.avai.vec] # nolint
    freq.dist <- sweep(count.dist, 1, rowSums(count.dist), "/")
    # table(cellInfo.tb[, c("meta.cluster", "loc")])/ as.vector(table(cellInfo.tb$meta.cluster)) # nolint
    freq.dist.bin <- floor(freq.dist * 100 / 10)
    print(freq.dist.bin)

    {
        count.dist.melt.ext.tb <- test.dist.table(count.dist)
        p.dist.tb <- dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "p.value") # nolint
        OR.dist.tb <- dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "OR") # nolint
        OR.dist.mtx <- as.matrix(OR.dist.tb[, -1])
        rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
    }

    sscVis::plotMatrix.simple(round(OR.dist.mtx, 2),
        out.prefix = sprintf("%s.OR.dist", out.prefix),
        show.number = T,
        waterfall.row = T,
        par.warterfall = list(score.alpha = 2, do.norm = T),
        exp.name = expression(italic(OR)),
        z.hi = 3,
        palatte = viridis::viridis(10),
        pdf.width = pdf.width, pdf.height = pdf.height
    )
    if (verbose == 1) {
        return(list(
            "count.dist.melt.ext.tb" = count.dist.melt.ext.tb,
            "p.dist.tb" = p.dist.tb,
            "OR.dist.tb" = OR.dist.tb,
            "OR.dist.mtx" = OR.dist.mtx
        ))
    } else {
        return(OR.dist.mtx)
    }
}

test.dist.table <- function(count.dist, min.rowSum = 0) {
    count.dist <- count.dist[rowSums(count.dist) >= min.rowSum, , drop = F]
    sum.col <- colSums(count.dist)
    sum.row <- rowSums(count.dist)
    count.dist.tb <- as.data.frame(count.dist)
    setDT(count.dist.tb, keep.rownames = T)
    count.dist.melt.tb <- melt(count.dist.tb, id.vars = "rn")
    colnames(count.dist.melt.tb) <- c("rid", "cid", "count")
    count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i) { # nolint
        this.row <- count.dist.melt.tb$rid[i]
        this.col <- count.dist.melt.tb$cid[i]
        this.c <- count.dist.melt.tb$count[i]
        other.col.c <- sum.col[this.col] - this.c
        this.m <- matrix(
            c(
                this.c,
                sum.row[this.row] - this.c,
                other.col.c,
                sum(sum.col) - sum.row[this.row] - other.col.c
            ),
            ncol = 2
        )
        res.test <- fisher.test(this.m)
        data.frame(
            rid = this.row,
            cid = this.col,
            p.value = res.test$p.value,
            OR = res.test$estimate
        )
    }))
    count.dist.melt.ext.tb <- merge(count.dist.melt.tb, count.dist.melt.ext.tb, # nolint
        by = c("rid", "cid")
    )
    count.dist.melt.ext.tb[, adj.p.value := p.adjust(p.value, "BH")]
    return(count.dist.melt.ext.tb)
}

OR.B.list <- do.tissueDist(
    cellInfo.tb = data@meta.data, meta.cluster = data$anno1, loc = data$group, # nolint
    out.prefix = "Figure1/",
    pdf.width = 10, pdf.height = 8, verbose = 1
)

## OR value
or <- round(OR.B.list$OR.dist.mtx, 2)
colnames(or) <- colnames(OR.B.list$OR.dist.mtx)

p11 <- sscVis::plotMatrix.simple(or,
    out.prefix = "Figure1/Figure_1E_Roe_major_cluster.pdf",
    mytitle = " ",
    show.number = T,
    do.clust = F,
    clust.column = F,
    par.heatmap = list(
        row_names_gp = grid::gpar(
            fontsize = 0
        ),
        column_names_gp = grid::gpar(
            fontsize = 0
        )
    ),
    clust.row = F,
    show.dendrogram = T,
    waterfall.row = F,
    returnHT = T,
    par.warterfall = list(score.alpha = 2, do.norm = T),
    exp.name = expression(italic(OR)),
    z.hi = 2,
    palatte = c("#F7EFF3", "#C9CADB", "#7A9BB7", "#3D6896", "#253D58"),
    pdf.width = 3, pdf.height = 5
)
p11 <- as.ggplot(p11) + theme(legend.position = "bottom")
ggsave("Figure1/Figure_1E_Roe_major_cluster.pdf", p11, width = 3, height = 7) # nolint
##################################################################


##################################################################
## Figure 1F Blood myeloid faction statistics
# ---------- 1) 读取元数据（更稳健） ----------
meta_all <- data[[]] |> tibble::as_tibble()

# 检查必需列
required_cols <- c("orig.ident", "anno1", "group")
missing_cols <- setdiff(required_cols, colnames(meta_all))
if (length(missing_cols) > 0) {
    stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", ")) # nolint
}

# ---------- 2) 计算每样本各谱系所占比例（分母=全体细胞） ----------
denom <- meta_all |>
    dplyr::count(.data$orig.ident, name = "n_total")

num <- meta_all |>
    dplyr::count(.data$orig.ident, .data$anno1, name = "n")

fractions1 <- num |>
    dplyr::inner_join(denom, by = "orig.ident") |>
    dplyr::mutate(frac = .data$n / .data$n_total) |>
    dplyr::left_join(
        meta_all |> dplyr::distinct(.data$orig.ident, .data$group),
        by = "orig.ident"
    )
# 如需固定组别顺序，可启用下一行（示例：HC, CAMR）
# fractions1 <- fractions1 |> dplyr::mutate(group = forcats::fct_relevel(factor(group), "HC", "CAMR")) # nolint

# ---------- 3) 小函数：按任意 anno1 谱系快速出图 ----------
plot_fraction_anno1 <- function(lineage = "Myeloid",
                                y_limits = c(0.05, 0.62),
                                palette = mycolor$group,
                                jitter_w = 0.12) {
    fractions1 |>
        dplyr::filter(.data$anno1 == lineage) |>
        ggplot(aes(x = group, y = frac)) +
        geom_boxplot(
            outlier.colour = "white", outlier.size = 0, width = 0.5,
            linetype = "dashed", outlier.alpha = 0, linewidth = 0.2
        ) +
        stat_boxplot(
            aes(ymin = ..lower.., ymax = ..upper..),
            outlier.size = 0,
            outlier.colour = "white", width = 0.5, linewidth = 0.2
        ) +
        stat_boxplot(
            geom = "errorbar", aes(ymin = ..ymax..), width = 0.2,
            linewidth = 0.2
        ) +
        stat_boxplot(
            geom = "errorbar", aes(ymax = ..ymin..), width = 0.2,
            linewidth = 0.2
        ) +
        geom_point(
            aes(fill = group),
            position = position_auto(0.4),
            shape = 21, size = 1, stroke = NA
        ) +
        scale_fill_manual(values = palette) +
        scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0.05, 0.1))) + # nolint
        theme_classic() +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            axis.line = element_line(linewidth = 0.15),
            axis.ticks = element_line(linewidth = 0.15)
        )
}

p_any <- plot_fraction_anno1("Myeloid")
ggsave("Figure1/Figure_1F_blood_Myeloid_fraction.pdf", p_any, width = 1, height = 1.6) # nolint
##################################################################


##################################################################
##  Figure 1H Tissue myeloid fraction statistics
# ---------- 0) 参数 ----------
wanted_lineages <- c("B cell", "T cell", "NK cell", "Plasma", "Myeloid")
group_levels <- c("HC", "CAMR") # 组的显示顺序

# ---------- 1) 读出元数据并做健壮性检查 ----------
meta_all <- t_data[[]] |> tibble::as_tibble()

required_cols <- c("orig.ident", "anno3", "group")
missing_cols <- setdiff(required_cols, colnames(meta_all))
if (length(missing_cols) > 0) {
    stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", ")) # nolint
}

# ---------- 2) 计算每样本各谱系所占比例（分母=全体细胞） ----------
denom <- meta_all |>
    dplyr::count(.data$orig.ident, name = "n_total")

num <- meta_all |>
    dplyr::filter(.data$anno3 %in% wanted_lineages) |>
    dplyr::count(.data$orig.ident, .data$anno3, name = "n")

fractions <- num |>
    dplyr::inner_join(denom, by = "orig.ident") |>
    dplyr::mutate(frac = .data$n / .data$n_total) |>
    dplyr::left_join(
        meta_all |>
            dplyr::distinct(.data$orig.ident, .data$group),
        by = "orig.ident"
    ) |>
    dplyr::mutate(
        group = forcats::fct_relevel(factor(.data$group), group_levels)
    )

# ---------- 3) 可复用小函数：一行切换谱系 ----------
plot_lineage_fraction <- function(lineage = "Myeloid",
                                  y_limits = c(0.01, 0.23),
                                  palette = mycolor$group,
                                  jitter_w = 0.12) {
    fractions |>
        dplyr::filter(.data$anno3 == lineage) |>
        ggplot(aes(x = group, y = frac)) +
        geom_boxplot(
            outlier.colour = "white", outlier.size = 0, width = 0.5,
            linetype = "dashed", outlier.alpha = 0, linewidth = 0.2
        ) +
        stat_boxplot(
            aes(ymin = ..lower.., ymax = ..upper..),
            outlier.size = 0,
            outlier.colour = "white", width = 0.5, linewidth = 0.2
        ) +
        stat_boxplot(
            geom = "errorbar", aes(ymin = ..ymax..), width = 0.2,
            linewidth = 0.2
        ) +
        stat_boxplot(
            geom = "errorbar", aes(ymax = ..ymin..), width = 0.2,
            linewidth = 0.2
        ) +
        geom_point(
            aes(fill = group),
            position = position_auto(0.4),
            shape = 21, size = 1, stroke = NA
        ) +
        scale_fill_manual(values = palette) +
        scale_y_continuous(
            limits = y_limits, expand = expansion(mult = c(0.05, 0.1))
        ) +
        theme_classic() +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            axis.line = element_line(linewidth = 0.15),
            axis.ticks = element_line(linewidth = 0.15)
        )
}

p_any <- plot_lineage_fraction("Myeloid")
ggsave("Figure1/Figure_1H_tissue_Myeloid_fraction.pdf", p_any, width = 1, height = 1.6) #  nolint
##################################################################


##################################################################
## Supplementary Figure 1A-B Expressed gene count & cell count

# -------------------- 公共准备：元数据与分组顺序 --------------------
group_levels <- c("HC", "SAF", "CAMR")

meta <- data[[]] |>
    tibble::as_tibble() |>
    mutate(
        group = str_trim(as.character(group)),
        group = factor(group, levels = group_levels),
        orig.ident = as.character(orig.ident)
    )

# -------------------- p4: nFeature_RNA 箱线图（组内按中位数排序） --------------------
# 1) 计算“样本在各自组内的 nFeature_RNA 中位数”，并按 group_levels + 中位数排序，得到全局 factor 顺序
sample_order_p4 <- meta |>
    group_by(group, orig.ident) |>
    summarize(median_nFeature = median(nFeature_RNA, na.rm = TRUE), .groups = "drop") |> # nolint
    arrange(match(group, group_levels), median_nFeature) |>
    pull(orig.ident)

# 2) 把排序映射回原数据（x 轴使用该因子水平；facet 使用 group）
meta_p4 <- meta |>
    mutate(orig.ident = factor(orig.ident, levels = sample_order_p4))

global_median_nFeature <- median(meta$nFeature_RNA, na.rm = TRUE)

p4 <- ggplot(meta_p4, aes(x = orig.ident, y = nFeature_RNA, color = group)) +
    geom_hline(yintercept = global_median_nFeature, linetype = "dashed", color = "gray") + # nolint
    geom_boxplot(
        width = 0.5, linetype = "dashed", linewidth = 0.2,
        outlier.shape = NA
    ) +
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.size = 0, outlier.colour = NA, width = 0.5) + # nolint
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.2) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.2) +
    facet_grid(~group, scales = "free_x", space = "free") +
    labs(y = "Expressed gene count") +
    scale_y_continuous(limits = c(0, 4000), expand = c(0, 0)) +
    scale_color_manual(values = mycolor$group) +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.3),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.3),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none"
    )

# -------------------- p9: 每样本细胞数柱状图（组内按细胞数排序） --------------------
# 1) 统计每组每样本细胞数，并按 group_levels + 计数排序，得到全局 factor 顺序
counts <- meta |>
    count(group, orig.ident, name = "Freq")

sample_order_p9 <- counts |>
    arrange(match(group, group_levels), Freq) |>
    pull(orig.ident)

# 2) 映射排序；全局中位数作为虚线
counts_p9 <- counts |>
    mutate(orig.ident = factor(orig.ident, levels = sample_order_p9))

global_median_count <- median(counts_p9$Freq, na.rm = TRUE)

p9 <- ggplot(counts_p9, aes(x = orig.ident, y = Freq, color = group, fill = group)) + # nolint
    geom_col(width = 0.7) +
    geom_hline(yintercept = global_median_count, linetype = "dashed", color = "gray") + # nolint
    scale_color_manual(values = mycolor$group) +
    scale_fill_manual(values = mycolor$group) +
    scale_y_continuous(limits = c(0, 19000), expand = c(0, 0)) +
    labs(y = "Cell count") +
    facet_grid(~group, scales = "free_x", space = "free") +
    theme_classic() +
    theme(
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.3),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank()
    )

final <- p4 + plot_spacer() + p9 + plot_layout(nrow = 1, widths = c(1, 0.2, 1))
ggsave("Figure1/Supplementary_Figure1_A-B_gene_cell_count.pdf", final, width = 8, height = 2) # nolint
##################################################################


##################################################################
## Supplementary Figure 1 C&D Dimplot for sample and group
p1 <- DimPlot(data,
    label = F, group.by = "orig.ident",
    cols = mycolor$sample, shuffle = T,
    pt.size = 0.1, raster = FALSE, seed = 3
) +
    NoAxes() +
    labs(title = "")

p2 <- DimPlot(data,
    label = F, group.by = "group",
    cols = mycolor$group, shuffle = T,
    pt.size = 0.1, raster = FALSE, seed = 3
) +
    NoLegend() + NoAxes() +
    labs(title = "")
final <- p1 + plot_spacer() + p2 + plot_layout(nrow = 1, widths = c(1, 1.2, 1))
ggsave("Figure1/Supplementary_Figure1_C_D_umap_sample_group.pdf", final, width = 25, height = 7) # nolint
##################################################################


##################################################################
## Supplementary Figure 1E featureplot marker
marker_gene <- c(
    "CD3D", "IL7R", "KLRD1", "KLRF1", "MS4A1",
    "CD79A", "JCHAIN", "MZB1", "LYZ", "FCN1"
) # nolint
pp <- FeaturePlot(
    data,
    features = marker_gene, cols = c("#F9FAFF", "#102B67"),
    pt.size = 0.0000001, ncol = 5, coord.fixed = T, raster = F,
    min.cutoff = 0, max.cutoff = 4
) & NoAxes() & labs(title = "") # nolint
p8 <- wrap_plots(pp, guides = "collect")
ggsave("Figure1/Supplementary_Figure1_E_featureplot.png", p8, width = 8, height = 4) # nolint
##################################################################
