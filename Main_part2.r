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
    "vegan", "forcats", "Nebulosa", "purrr", "Matrix", "pvclust"
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
mye_tdata <- readRDS("tissue_result/mye_tissue.rds")
mye_pan <- readRDS("../blood_mao/myeloid/mye_blood_mao_input.rds")
load("input_file/color.RData")
##################################################################


##################################################################
## Figure 2A DEG circle plot
# ---- helpers ----
read_deg <- function(path_points, path_bar, drop_cell = "m10_Neutrophil") {
    pts <- readr::read_tsv(path_points, show_col_types = FALSE) |>
        filter(.data$cell != drop_cell) |>
        mutate(cell = forcats::fct_relevel(.data$cell, levels(mye)))
    bar <- readr::read_tsv(path_bar, show_col_types = FALSE) |>
        filter(.data$cell != drop_cell)
    list(points = pts, bar = bar)
}

# 统一的圆环火山样式（保持你原样式不变）
circle_plot <- function(points, bar, angle_vec,
                        fill_pal = mycolor$mye_type,
                        y_lim = c(-8, 4)) {
    # 与原逻辑一致的索引修正
    points <- points |>
        mutate(
            num = if_else(.data$num > 10, .data$num - 1L, .data$num),
            x   = if_else(.data$x > 10, .data$x - 1L, .data$x)
        )

    bar <- bar |>
        mutate(num = if_else(.data$num > 10, .data$num - 1L, .data$num))

    bar_low <- bar |>
        filter(.data$variable == "low") |>
        mutate(angle = angle_vec)

    ggplot(points, aes(x = .data$cell, y = .data$y)) +
        # 极小散点
        geom_point(aes(x = .data$x, y = .data$y, color = .data$sig), size = 0.00001) + # nolint
        scale_color_manual(values = c("#749E6E", "#DDB26A")) +
        guides(color = "none") +
        # 灰色背景环
        geom_bar(
            data = bar, aes(x = .data$num, y = .data$value),
            stat = "identity", position = "stack",
            alpha = 0.3, fill = "#b2b2b2b3"
        ) +
        new_scale_color() +
        # 彩色外环与刻度标签
        geom_tile(
            data = bar_low,
            aes(x = .data$num, y = 0, fill = .data$cell),
            alpha = 1, color = "black", height = 0.8
        ) +
        geom_text(
            data = bar_low,
            aes(x = .data$num, y = 0, label = .data$cell, angle = .data$angle),
            size = 1
        ) +
        scale_fill_manual(values = fill_pal) +
        geom_text_repel(
            data = points |> filter(.data$sig == "up"),
            aes(x = .data$x, y = .data$y, label = .data$label),
            na.rm = TRUE, size = 1,
            point.padding = NA, box.padding = 0,
            segment.size = 0.2, fontface = "bold", force = 5,
            direction = "both", nudge_y = 1.5
        ) +
        geom_text_repel(
            data = points |> filter(.data$sig == "down"),
            aes(x = .data$x, y = .data$y, label = .data$label),
            na.rm = TRUE, size = 1,
            point.padding = NA, box.padding = 0,
            segment.size = 0.2, fontface = "bold", force = 6,
            direction = "both", nudge_y = -1.5
        ) +
        guides(color = "none", fill = "none") +
        coord_polar(start = 0) +
        ylim(y_lim) +
        theme(
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.x.bottom = element_blank()
        )
}

# ---- data in ----
camr_saf <- read_deg(
    "input_file/mye_camr_saf_deg.txt",
    "input_file/bardata_saf_camr_deg.txt"
)
camr_hc <- read_deg(
    "input_file/mye_camr_hc_deg.txt",
    "input_file/bardata_hc_camr_deg.txt"
)

# ---- DEG count（保持原有排序/配色/坐标范围）----
genes_of_interest <- c(
    "FOS", "DDIT4", "GNAS", "FOSB", "FKBP5", "TLR2",
    "SMAP2", "PER1", "JUND", "CXCR4", "JUN", "JDP2",
    "JAK2", "DUSP1", "CEBPD", "ZEB2", "CEBPB"
)

# saf / hc 的计数
cnt_saf <- camr_saf$points |>
    count(rowname, name = "saf") |>
    filter(rowname %in% genes_of_interest)

cnt_hc <- camr_hc$points |>
    count(rowname, name = "nc") |>
    filter(rowname %in% genes_of_interest)

# 合并并按 saf 频次排序
mat <- cnt_saf |>
    full_join(cnt_hc, by = "rowname") |>
    mutate(
        saf = replace_na(saf, 0L),
        nc  = replace_na(nc, 0L)
    ) |>
    filter(rowname %in% genes_of_interest) |>
    arrange(desc(saf)) |>
    mutate(
        gene = forcats::fct_rev(forcats::fct_relevel(rowname, rowname)),
        saf  = -as.numeric(saf),
        nc   = as.numeric(nc)
    )

p_left <- ggplot(mat, aes(x = .data$saf, y = .data$gene)) +
    geom_col(fill = "#A64447", color = "black", width = 0.7, linewidth = 0.3) +
    theme_classic() +
    scale_y_discrete(position = "right") +
    scale_x_continuous(limits = c(-9, 0), breaks = c(-9, -6, -3, 0)) +
    theme(axis.text = element_blank(), axis.title = element_blank())

p_right <- ggplot(mat, aes(x = .data$nc, y = .data$gene)) +
    geom_col(fill = "#6093B5", color = "black", width = 0.7, linewidth = 0.3) +
    theme_classic() +
    theme(axis.text = element_blank(), axis.title = element_blank())

final_degcount <- p_left + plot_spacer() + p_right +
    plot_layout(nrow = 1, widths = c(1, 0.3, 1))

ggsave("Figure2/Figure_2C_deg_count.pdf", final_degcount, width = 3.3, height = 3) # nolint

# ---- Figure 2A/2B：圆环图 ----
myAngle <- seq(-10, -350, length.out = length(levels(mye)))

p_2A <- circle_plot(
    points = camr_saf$points,
    bar = camr_saf$bar,
    angle_vec = myAngle,
    fill_pal = mycolor$mye_type,
    y_lim = c(-8, 4)
)
ggsave("Figure2/Figure_2A_CAMRvsSAF_deg_circle.pdf", p_2A, height = 7, width = 7) # nolint

p_2B <- circle_plot(
    points = camr_hc$points,
    bar = camr_hc$bar,
    angle_vec = myAngle,
    fill_pal = mycolor$mye_type,
    y_lim = c(-8, 4)
)
ggsave("Figure2/Figure_2B_CAMRvsHC_deg_circle.pdf", p_2B, height = 7, width = 7)
##################################################################


##################################################################
## Figure 3A Myeloid subcluster Ro/e
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
    count.dist.melt.ext.tb <- as.data.table(plyr::ldply(seq_len(nrow(count.dist.melt.tb)), function(i) { # nolint
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
    cellInfo.tb = mye@meta.data, meta.cluster = mye$anno2, loc = mye$group, # nolint
    out.prefix = "Figure3/",
    pdf.width = 10, pdf.height = 8, verbose = 1
)

## OR value
or <- round(OR.B.list$OR.dist.mtx, 2)
colnames(or) <- colnames(OR.B.list$OR.dist.mtx)

p11 <- sscVis::plotMatrix.simple(or,
    out.prefix = "Figure3/Figure_3A_Roe_mye_cluster.pdf",
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
ggsave("Figure3/Figure_3A_Roe_mye_cluster.pdf", p11, width = 3, height = 7) # nolint
##################################################################


##################################################################
## Figure 3B Fraction of SLC8A1+ macrophage in tissue data
# ---- 1) 频率矩阵 ----
meta_df <- as_tibble(mye_tdata@meta.data) |>
    select(orig.ident, anno5)

mat <- meta_df |>
    count(orig.ident, anno5, name = "n") |>
    group_by(orig.ident) |>
    mutate(value = n / sum(n)) |>
    ungroup()

# ---- 2) 样本分组 ----
# 以出现顺序获取样本列表，然后生成同长度分组向量并 join 回去
id_order <- meta_df |>
    distinct(orig.ident) |>
    mutate(idx = row_number())

group_map <- id_order |>
    transmute(
        orig.ident,
        group = factor(
            rep(c("CAMR", "HC"), each = 7)[idx], # 与你原来的 rep(...) 等价
            levels = c("HC", "CAMR")
        )
    )

mat <- mat |>
    left_join(group_map, by = "orig.ident")

# ---- 3) 统一的箱线图函数 ----
make_box <- function(df, ymax) {
    ggplot(df, aes(x = group, y = value)) +
        geom_boxplot(
            outlier.colour = "white", outlier.size = 0,
            width = 0.5, linetype = "dashed", outlier.alpha = 0, linewidth = 0.2
        ) +
        stat_boxplot(
            aes(ymin = ..lower.., ymax = ..upper..),
            outlier.size = 0, outlier.colour = "white",
            width = 0.5, linewidth = 0.2
        ) +
        stat_boxplot(
            geom = "errorbar", aes(ymin = ..ymax..),
            width = 0.2, linewidth = 0.2
        ) +
        stat_boxplot(
            geom = "errorbar", aes(ymax = ..ymin..),
            width = 0.2, linewidth = 0.2
        ) +
        geom_point(
            aes(fill = group),
            position = position_auto(0.4),
            shape = 21, size = 1, stroke = NA
        ) +
        scale_fill_manual(values = mycolor$group) +
        stat_compare_means() +
        scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0.05, 0.1))) + # nolint
        theme_classic() +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            axis.line = element_line(linewidth = 0.15),
            axis.ticks = element_line(linewidth = 0.15)
        )
}

# ---- 4) 两个面板 ----
mt1 <- mat |> filter(anno5 == "ZEB2+ Macrophage")
p2 <- make_box(mt1, ymax = 0.4)

mt2 <- mat |> filter(anno5 == "ZEB2+CD163+ Macrophage")
p3 <- make_box(mt2, ymax = 0.2)

final <- p2 + plot_spacer() + p3 + plot_layout(nrow = 1, widths = c(1, 0.1, 1, 0.1, 1)) # nolint
ggsave("Figure3/Figure_3B_SLC8A1_macrophage_tissue_fraction.pdf", final, width = 8, height = 6) # nolint
##################################################################


##################################################################
# ---- data ----
dat <- subset(mye, idents = levels(mye)[7:9])
gene1 <- c("IDO1", "IDO2", "CADM1", "CD59", "BTLA", "THBD", "HAVCR2", "ID2")
gene2 <- c("PHB2", "CD33", "FCGR2B", "LILRB4", "IRF7", "IRF8", "TGFB1", "LILRA4") # nolint

# 统一调色
pal3 <- c("#1078a9CC", "#dd626fCC", "#815AA8CC")

# ---- helpers ----
vln_theme_compact <- function() {
    theme_classic() +
        theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_line(linewidth = 0.5, color = "black"),
            axis.ticks = element_line(linewidth = 0.5, color = "black")
        )
}

make_vln <- function(obj, genes, cols = pal3) {
    VlnPlot(
        obj,
        features = genes,
        stack = TRUE, sort = FALSE, group.by = "anno2",
        flip = TRUE, cols = cols, pt.size = 0, split.by = "anno2"
    ) + vln_theme_compact()
}

# ---- plots ----
p1 <- make_vln(dat, gene1)
p2 <- make_vln(dat, gene2)

final <- p1 + plot_spacer() + p2 + plot_layout(nrow = 1, widths = c(1, 0.3, 1))
ggsave("Figure3/Figure_3H_dc_violin.pdf", final, width = 9.6, height = 4)
##################################################################


##################################################################
## Supplementary Figure 4G FeaturePlot for autograpy-related genes
# ---- 1) 数据子集 + UMAP ----
mye_sub <- subset(mye, idents = levels(mye)[1:6])
mye_sub <- RunUMAP(mye_sub, reduction = "mnn", dims = 1:30)

p_umap <- DimPlot(mye_sub, cols = mycolor$mye_type, label = FALSE) + NoLegend()
ggsave("Figure3/Supplementary_Figure_4G_dimplot.pdf", p_umap, width = 12, height = 12) # nolint

# ---- 2) 密度小提琴底图封装 ----
Edensity <- function(data, features,
                     reduction = "umap",
                     method = "wkde",
                     pt_size = 0.2,
                     ncol = 2,
                     palette = c("lightblue", "lightyellow", "#EF3B36")) {
    plots <- map(features, ~
        Nebulosa::plot_density(
            data,
            reduction = reduction,
            .x,
            method = method,
            size = pt_size
        ) +
            coord_fixed() +
            theme_void() +
            theme(
                plot.title = element_blank(),
                legend.position = "none"
            ) +
            scale_color_gradientn(colours = palette))

    wrap_plots(plots, ncol = ncol)
}

# ---- 3) 生成并保存 ----
p_den <- Edensity(
    data = mye_sub,
    features = c("ATG5", "ATG7", "VMP1", "VPS13C", "VPS8", "DAPK1", "MAPK8", "PRKAG2") # nolint
)
ggsave("Figure3/Supplementary_Figure_4G_autograpy_gene.pdf", p_den, width = 5, height = 8) # nolint
##################################################################


##################################################################
# ---- 并行设置（保持不变）----
plan("multicore", workers = 60)
options(future.globals.maxSize = 100000 * 1024^2)

# ---- 1) Marker 计算与保存（保持逻辑与阈值不变）----
marker <- FindMarkers(
    mye,
    ident.1  = "m06_ZEB2+ monocyte",
    only.pos = TRUE,
    min.pct  = 0.25
)

marker %>%
    rownames_to_column("gene") %>%
    readr::write_tsv("temp/m06_marker.txt")

# ---- 2) 筛选显著基因并做 ID 转换 ----
p_row <- marker %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    filter(p_val_adj < 0.05)

eg <- bitr(p_row$gene,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = "org.Hs.eg.db"
)

p_row1 <- p_row %>%
    left_join(eg, by = c("gene" = "SYMBOL")) %>%
    drop_na(ENTREZID)

# ---- 3) 制作 TERM2GENE：GO:BP + HALLMARK ----
m_t2g <- msigdbr(species = "Homo sapiens") %>%
    filter(gs_subcat == "GO:BP" | gs_cat == "H") %>%
    select(gs_name, entrez_gene)

# （可选）快速查看 GOBP 出现的类别
# msigdbr(species = "Homo sapiens") %>%
#   filter(str_detect(gs_name, "GOBP")) %>%
#   count(gs_cat)

# ---- 4) 富集分析（与原 enricher 等价）----
em <- enricher(p_row1$ENTREZID, TERM2GENE = as.data.frame(m_t2g))

# enrichResult 取结果表更明确：em@result
x <- em@result %>%
    as_tibble() %>%
    mutate(Description = stringr::str_to_lower(Description))

readr::write_csv(x, "temp/m06_hallmarker_profile.csv")

# ---- 5) 选定通路并构造绘图数据 ----
select_ids <- c(
    "GOBP_REGULATION_OF_GTPASE_ACTIVITY",
    "GOBP_REGULATION_OF_SMALL_GTPASE_MEDIATED_SIGNAL_TRANSDUCTION",
    "GOBP_HISTONE_MODIFICATION",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_MITOTIC_SPINDLE",
    "GOBP_PEPTIDYL_SERINE_MODIFICATION",
    "HALLMARK_UV_RESPONSE_DN",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "GOBP_MYELOID_CELL_DIFFERENTIATION",
    "GOBP_REGULATION_OF_RAS_PROTEIN_SIGNAL_TRANSDUCTION",
    "GOBP_MACROAUTOPHAGY",
    "GOBP_REGULATION_OF_WNT_SIGNALING_PATHWAY",
    "GOBP_REGULATION_OF_AUTOPHAGY",
    "GOBP_POSITIVE_REGULATION_OF_CELL_ADHESION",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "GOBP_DEPHOSPHORYLATION",
    "GOBP_FC_RECEPTOR_SIGNALING_PATHWAY",
    "GOBP_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
    "GOBP_JNK_CASCADE"
)

x_sub <- x %>%
    filter(ID %in% select_ids) %>%
    # 维持“按当前行顺序再反转”的因子顺序
    mutate(Description = factor(Description, levels = rev(Description))) %>%
    # group
    mutate(group = sub("_.+$", "", ID)) %>%
    # ratio
    tidyr::separate_wider_delim(GeneRatio, "/",
        names = c("a", "b"), cols_remove = FALSE
    ) %>%
    mutate(
        a = as.numeric(a),
        b = as.numeric(b),
        ratio = a / b
    )

# ---- 6) 绘图 ----
p <- ggplot(x_sub, aes(x = -log10(qvalue), y = Description)) +
    geom_bar(aes(color = ratio, fill = ratio),
        width = 0.11, stat = "identity"
    ) +
    geom_point(aes(shape = group, color = ratio), size = 6) +
    scale_color_viridis() +
    scale_fill_viridis() +
    scale_x_continuous(position = "top") +
    scale_shape_manual(values = c(19, 18)) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.border     = element_blank(),
        axis.ticks       = element_blank(),
        axis.text.y      = element_blank(),
        axis.title       = element_blank(),
        legend.position  = "none"
    )

ggsave("Figure3/Figure_3G_SLC8A1_GO.pdf", p, width = 2.4, height = 5) # nolint

## Figure 3E SLC8A1+ monocyte highly expressed genes
lab <- c("SLC8A1", "TLR2", "VPS13B", "NFKBIZ", "ATG7", "ZEB2", "JAK2", "NOTCH2", "NLRP3", "NFKB1") # nolint

marker_tbl <- marker %>%
    rownames_to_column("gene") %>%
    mutate(
        p_val_adj = if_else(p_val_adj == 0, 9.888883e-303, p_val_adj)
    ) %>%
    arrange(avg_log2FC) %>%
    mutate(
        rank  = row_number(),
        label = if_else(gene %in% lab, gene, NA_character_)
    )

p <- ggplot(marker_tbl, aes(x = rank, y = avg_log2FC)) +
    geom_rect(xmin = 797, xmax = 1100, ymin = 0, ymax = 1, fill = "#CBAAAB") +
    geom_point(aes(color = -log10(p_val_adj)), size = 0.3) +
    scale_color_viridis() +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.title = element_blank()
    )

ggsave("Figure3/Figure_3E_highly_expressed_gene.pdf", p, width = 3.8, height = 2) # nolint
##################################################################


##################################################################
## Figure 3D Fraction of myeloid in big blood data
# 1) 计算各 group 内的比例
meta_df <- as_tibble(mye_pan@meta.data) %>%
    select(group, anno) %>%
    mutate(
        group = factor(group, levels = c("N", "T", "I", "A")),
        anno  = factor(anno, levels = anno_levels)
    )

mat <- meta_df %>%
    count(group, anno, name = "n") %>%
    group_by(group) %>%
    mutate(value = n / sum(n)) %>%
    ungroup() %>%
    mutate(anno = fct_rev(anno))

# 2) 作图
p <- ggplot(mat, aes(x = value, y = anno, fill = group)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    theme_classic() +
    scale_fill_manual(values = pal_futurama(alpha = 0.8)(4)[c(3, 4, 1, 2)]) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(axis.line.y = element_blank())

ggsave("Figure3/Figure_3D_pan_myeloid_percent.pdf", p, width = 4.6, height = 3)
##################################################################


##################################################################
## Supplementary Figure 4F pvcluster tree
# ---- 1) 子集 ----
mye_sub <- subset(mye, idents = levels(mye)[6])
mye_sub$batch <- "batch2"
mye_sub$anno <- mye_sub$anno2

# ---- 2) 下采样并合并（保持参数不变）----
mye_pan@active.ident <- factor(mye_pan$anno)
mye_pan_sub <- subset(mye_pan, downsample = 1000)
mye_pan_sub$batch <- "batch1"

merge_data <- merge(mye_sub, mye_pan_sub)

merge_data <- NormalizeData(
    merge_data,
    normalization.method = "LogNormalize",
    scale.factor = 10000
)
merge_data <- FindVariableFeatures(
    merge_data,
    selection.method = "vst",
    nfeatures = 2000
)
all.genes <- rownames(merge_data)
merge_data <- ScaleData(merge_data, features = all.genes)
merge_data <- RunPCA(merge_data, features = VariableFeatures(object = merge_data)) # nolint
merge_data <- RunFastMNN(object.list = SplitObject(merge_data, split.by = "batch")) # nolint

# 取 MNN 重构矩阵
mat <- GetAssayData(merge_data, assay = "mnn.reconstructed", slot = "data")

# ---- 3) 计算前 2000 行在各 anno 分组的平均值（向量化，避免 for 循环）----
mat_sub <- mat[seq_len(2000), , drop = FALSE]

g <- factor(merge_data$anno)
# 构造组指示矩阵（稀疏）
G <- sparse.model.matrix(~ g - 1) # 列名如 gLevel
colnames(G) <- levels(g)

# 组内均值 = mat_sub %*% G / 每组细胞数
group_sizes <- colSums(G)
mtt <- (mat_sub %*% G) %*% Diagonal(x = 1 / as.numeric(group_sizes))

# 行名与列名
rownames(mtt) <- rownames(mat_sub)
colnames(mtt) <- colnames(G) # 即 levels(g)

# ---- 4) pvclust ----
set.seed(1) # 可选：确保可复现
result <- pvclust(
    as.matrix(mtt),
    method.dist   = "cor",
    method.hclust = "average",
    nboot         = 1000,
    parallel      = TRUE
)

pdf("Figure3/Supplementary_Figure_4F_pvcluster_tree.pdf", width = 7, height = 9)
plot(result)
pvrect(result, alpha = 0.95)
dev.off()
##################################################################


##################################################################
## Figure 3F autograpy gene expression
# ---- 1) 取表达矩阵（RNA@data）与基因子集 ----
genes <- c(
    "SLC8A1", "ZEB2", "TLR2", "VPS13B", "ATG7", "NFKB1", "NOTCH2", "NLRP3", # nolint
    "JAK2", "ATG5", "VPS8", "MAPK8", "VPS13C", "DAPK1", "PRKAG2"
)

expr <- GetAssayData(mye_pan, assay = "RNA", slot = "data") # dgCMatrix
expr <- expr[genes, , drop = FALSE]

# ---- 2) 构造分组因子与顺序（与原顺序完全一致）----
anno_levels <- levels(mye_pan$anno)
g <- factor(mye_pan$anno, levels = anno_levels)

# 组指示稀疏矩阵：每列为一个 anno 组
G <- sparse.model.matrix(~ g - 1) # 列名类似 gLevel
colnames(G) <- levels(g)
group_sizes <- colSums(G)

# ---- 3) 计算组均值与表达比例（与 aggregate(mean) 和自定义 pos_cal 等价）----
# 平均表达： (expr %*% G) / n_g
mean_mat <- (expr %*% G) %*% Diagonal(x = 1 / as.numeric(group_sizes))
rownames(mean_mat) <- rownames(expr) # genes
colnames(mean_mat) <- levels(g) # anno levels

# 表达比例： ((expr > 0) %*% G) / n_g
pos_mat <- ((expr > 0) %*% G) %*% Diagonal(x = 1 / as.numeric(group_sizes))
rownames(pos_mat) <- rownames(expr)
colnames(pos_mat) <- levels(g)

# ---- 4) 整理为长表 ----
pos_long <- melt(as.matrix(pos_mat))
mean_long <- melt(as.matrix(mean_mat))

mer_mat <- cbind(pos_long, mean_long)
colnames(mer_mat)[c(3, 6)] <- c("pos_value", "mean_value")
mer_mat <- mer_mat[, 3:6] # 只保留 Var1 Var2 pos_value mean_value

# y 轴顺序（与原 factor(..., levels = rev(sort(...)[c(...)])) 一致）
mer_mat$Var2 <- factor(mer_mat$Var2, levels = rev(anno_levels))
# x 轴按基因给定顺序
mer_mat$Var1 <- factor(mer_mat$Var1, levels = genes)

# ---- 5) 作图 ----
p <- ggplot(mer_mat, aes(y = Var2, x = Var1)) +
    geom_point(aes(size = pos_value, fill = mean_value),
        color = "black", shape = 22
    ) +
    theme_bw() +
    scale_fill_gradientn(
        colors = c(rep("white", 10), colorRampPalette(c("white", "#ABAAAA"))(100)) # nolint
    ) +
    scale_size(range = c(0, 8)) +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave("Figure3/Figure_3F_gene_exp.pdf", p, width = 8, height = 4.2)
##################################################################
