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
bp <- readRDS("input_file/bp_sub.rds")
t_data <- readRDS("input_file/tissue_input.rds")
bp_bcr <- readRDS("input_file/bcr_processed_data.rds")
bp01_bcr <- readRDS("input_file/bp01_bcr.rds")
bp06_bcr <- readRDS("bp06_bcr.rds")
bp_tdata <- readRDS("tissue_result/tissue_bp.rds")
bp_pan <- readRDS("../blood_mao/b_cell/b_blood_mao_input.rds")
load("input_file/color.RData")
##################################################################


##################################################################
## Figure 6K DotPlot
pdf("Figure6/Figure_6K_dotplot.pdf", width = 5.2, height = 1.8)
DotPlot(bp, features = c("PPP3CA", "FKBP1A", "FKBP5"), cols = c("lightgrey", "red")) + # nolint
    coord_flip() +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1),
        axis.ticks.y = element_blank()
    )
dev.off()

pdf("Figure6/Supplementary_Figure_9B_dot.pdf", width = 5.2, height = 1.8)
DotPlot(bp, features = c("CD19", "ZEB2", "TBX21"), cols = c("lightgrey", "red")) + # nolint
    coord_flip() +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1),
        axis.ticks.y = element_blank()
    )
dev.off()
##################################################################


##################################################################
## Supplementary Figure 9A B & Plasma subcluster Ro/e
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
    cellInfo.tb = bp@meta.data, meta.cluster = bp$anno2, loc = bp$group, # nolint
    out.prefix = "Figure6/",
    pdf.width = 10, pdf.height = 8, verbose = 1
)

## OR value
or <- round(OR.B.list$OR.dist.mtx, 2)
colnames(or) <- colnames(OR.B.list$OR.dist.mtx)

p11 <- sscVis::plotMatrix.simple(or,
    out.prefix = "Figure6/Supplementary_Figure_9A_Roe_b_cluster.pdf",
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
ggsave("Figure6/Supplementary_Figure_9A_Roe_b_cluster.pdf", p11, width = 3, height = 7) # nolint
##################################################################


##################################################################
## Figure 6L GO analysis
# ==== 0) 前置：marker 计算（保持你的参数） ====
markers <- FindAllMarkers(bp, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1) # nolint

# ==== 1) 通用：基因 → GO 富集（ALL 本体），返回 data.frame ====
go_func <- function(genes_symbol) {
    df_name <- bitr(genes_symbol,
        fromType = "SYMBOL",
        toType   = c("ENTREZID"),
        OrgDb    = org.Hs.eg.db
    )
    if (nrow(df_name) == 0) {
        return(tibble())
    }
    go <- enrichGO(
        gene          = unique(df_name$ENTREZID),
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2
    )
    as_tibble(go@result)
}

# ==== 2) 通用：做一张和你风格一致的 GO 条形图 ====
plot_go <- function(go_df, terms_keep, panel_title = "") {
    # 仅保留需要显示的 term，按给定顺序，并只取 qvalue
    df <- go_df %>%
        select(Description, qvalue) %>%
        filter(Description %in% terms_keep) %>%
        mutate(Description = factor(Description, levels = rev(terms_keep))) %>%
        arrange(Description)

    # 若为空，给出一张空图避免中断
    if (nrow(df) == 0) {
        return(
            ggplot() +
                theme_void() +
                ggtitle(panel_title)
        )
    }

    # 渐变色长度与条目数一致（保持你原来的两端色）
    gradient_colors <- colorRampPalette(c("#F8E5D9", "#842827"))(nrow(df))

    ggplot(df, aes(x = Description, y = -log(qvalue, 10))) +
        geom_bar(aes(fill = Description), stat = "identity") +
        scale_fill_manual(values = gradient_colors) +
        geom_hline(yintercept = -log(0.05, 10), linetype = "dashed", color = "gray") + # nolint
        theme_bw() +
        coord_flip() +
        labs(x = "GO Biological Process", y = "Adjust p-value", title = panel_title) + # nolint
        theme(
            axis.text.x           = element_blank(),
            panel.grid.major.y    = element_blank(),
            legend.position       = "none",
            axis.text.y           = element_blank(),
            axis.ticks.y          = element_blank(),
            axis.title            = element_blank(),
            plot.title            = element_blank()
        )
}

# ==== 3) 为两个 cluster 定义要展示的术语（与你代码一致） ====
terms_bp01 <- c(
    "regulation of GTPase activity", "peptidyl-serine phosphorylation",
    "Ras protein signal transduction", "histone modification",
    "Fc receptor signaling pathway", "regulation of autophagy",
    "B cell differentiation", "Wnt signaling pathway",
    "stress-activated protein kinase signaling cascade",
    "macroautophagy", "JNK cascade", "endocytosis"
)

terms_bp06 <- c(
    "immune response-regulating signaling pathway",
    "immune response-regulating cell surface receptor signaling pathway",
    "antigen processing and presentation of peptide antigen",
    "B cell activation", "regulation of cell-cell adhesion",
    "B cell proliferation", "regulation of endocytosis",
    "phagocytosis", "tumor necrosis factor production",
    "Fc receptor signaling pathway",
    "regulation of small GTPase mediated signal transduction",
    "Ras protein signal transduction"
)

# ==== 4) 取每个 cluster 的显著基因，跑富集并出图 ====
genelist1 <- markers %>%
    filter(cluster == "bp01.Pro-BC", p_val_adj < 0.05) %>%
    pull(gene) %>%
    unique()

genelist2 <- markers %>%
    filter(cluster == "bp06.atMBC", p_val_adj < 0.05) %>%
    pull(gene) %>%
    unique()

go_res1 <- go_func(genelist1)
go_res2 <- go_func(genelist2)

p1 <- plot_go(go_res1, terms_bp01, panel_title = "05.ABC")
p2 <- plot_go(go_res2, terms_bp06, panel_title = "06.Mature BC")

# ==== 5) 拼图与保存（布局保持 1 行） ====
final <- p1 + p2 + plot_layout(nrow = 1)
ggsave("Figure6/Figure_6L_deg_go.pdf", final, width = 5, height = 4)
##################################################################


##################################################################
# 1) 子集与下采样（保持你的逻辑与随机种子）
bp_sub <- subset(bp, idents = levels(bp)[1:6])
set.seed(100)
bp_sub <- subset(bp_sub, downsample = 1000)

# 2) 取 counts 与元数据（避免直接访问 @slots）
obj_counts <- GetAssayData(bp_sub, assay = "RNA", slot = "counts")
metadata <- bp_sub@meta.data

# 3) 计算样本纯度相关指标（与原流程一致）
ent.res <- SE_fun(obj_counts)
rogue.value <- CalculateRogue(ent.res, platform = "UMI")

rogue.res <- rogue(
    obj_counts,
    labels = metadata$anno2,
    samples = metadata$orig.ident,
    platform = "UMI",
    span = 1.2
)

# 4) 整理为长表（替代 melt；稳健去 NA；因子顺序与 bp_sub$anno2 一致）
rogue_df <- rogue.res %>%
    as.data.frame() %>%
    rownames_to_column("rowname") %>%
    pivot_longer(
        cols = -rowname,
        names_to = "variable",
        values_to = "ROGUE"
    ) %>%
    drop_na(ROGUE) %>%
    mutate(
        variable = factor(variable, levels = sort(unique(bp_sub$anno2)))
    )

# 5) 画图（样式保持不变；中位数线等同）
p <- ggplot(rogue_df, aes(x = variable, y = ROGUE, color = variable)) +
    geom_boxplot(width = 0.6, outlier.colour = "white", outlier.size = 0) +
    scale_color_manual(values = mycolor$bp_type) +
    geom_hline(
        yintercept = median(rogue_df$ROGUE, na.rm = TRUE),
        linetype = "dashed", color = "gray"
    ) +
    theme_bw() +
    labs(title = "Abc") + # 轴标题下一行会被去掉，视觉与原来一致
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", angle = 30, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, color = "black")
    )

ggsave("Figure6/Figure_6I_rogue.pdf", p, width = 2.9, height = 2.7)
##################################################################


##################################################################
## Figure 6J
do.tissueDist <- function(cellInfo.tb = cellInfo.tb, # meta.data
                          meta.cluster = cellInfo.tb$meta.cluster, # 纵坐标，可以是不同的分群  # nolint
                          loc = cellInfo.tb$loc, # 不同的分组，可以是肿瘤，癌旁，正常等
                          out.prefix, # 输出文件的名字
                          pdf.width = 3,
                          pdf.height = 5,
                          verbose = 0 # 如果等于1，则返回对象，可以对对象稍作处理重新画图
) {
    ## input data
    library(data.table)
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

bp_01 <- subset(bp_bcr, anno2 == "bp01.Pro-BC")

mat <- bp_01@meta.data %>%
    as_tibble(rownames = "cell") %>%
    transmute(
        group,
        c_call,
        c_call1 = case_when(
            c_call %in% c("IGHA1", "IGHA2") ~ "IGHA",
            c_call %in% c("IGHG1", "IGHG2", "IGHG3") ~ "IGHG",
            TRUE ~ c_call
        ),
        # 固定显示顺序；即使某类当前不存在，也保留水平
        c_call1 = fct(c_call1, levels = c("IGHD", "IGHM", "IGHG", "IGHA"))
    ) %>%
    arrange(c_call1)

OR.B.list <- do.tissueDist(
    cellInfo.tb = mat, meta.cluster = mat$c_call1, loc = mat$group, # nolint
    out.prefix = "Figure6/",
    pdf.width = 10, pdf.height = 8, verbose = 1
)

## OR value
or <- round(OR.B.list$OR.dist.mtx, 2)
colnames(or) <- colnames(OR.B.list$OR.dist.mtx)
or <- or[c(2, 4, 3, 1), ]
p11 <- sscVis::plotMatrix.simple(or,
    out.prefix = "Figure6/Figure_6J_Roe_bcr.pdf",
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
    palatte = rev(brewer.pal(9, "YlGnBu")),
    pdf.width = 3, pdf.height = 5
)
p11 <- as.ggplot(p11) + theme(legend.position = "bottom")
ggsave("Figure6/Figure_6J_Roe_bcr.pdf", p11, width = 3, height = 7) # nolint
##################################################################


##################################################################
# -------- UMAP 按 c_call 上色（更稳健的因子重排 + 颜色长度自适应）--------
# 原：unique(...)[c(3,2,5,4,1,6,7)] 太脆弱；下面写法在缺失/新增水平时也不崩
lv <- levels(factor(bp06_bcr$c_call))
idx <- intersect(c(3, 8, 4:7, 1:2), seq_along(lv))
bp06_bcr$c_call <- factor(bp06_bcr$c_call, levels = lv[idx])

cols0 <- c(
    "#7F7CB9", "#F09DB6", "#C7DEEF", "#70AED6", "#3272B4",
    "#213367", "#F6C55F", "#F0993E"
)
# 颜色与现有水平数对齐（多余截断，不足则循环补齐）
n_lv <- nlevels(bp06_bcr$c_call)
cols <- if (length(cols0) >= n_lv) cols0[seq_len(n_lv)] else rep_len(cols0, n_lv) # nolint

pdf("Figure6/Supplementary_Figure_9E_bp_iso.pdf", width = 2, height = 2)
DimPlot(
    bp06_bcr,
    label = FALSE, group.by = "c_call",
    cols = cols,
    pt.size = 0.1, raster = FALSE
) + NoLegend() + NoAxes() + labs(title = "")
dev.off()

# -------- 更健壮的密度图函数（保持你的视觉参数不变）--------
Edensity <- function(data = seu_obj, features = genes) {
    p <- lapply(features, function(z) {
        p <- plot_density(data,
            reduction = "umap", z, method = "wkde",
            size = 0.2
        ) +
            coord_fixed() + theme_void() +
            theme(legend.position = "none", plot.title = element_blank()) + # nolint
            scale_color_gradientn(colours = c("lightblue", "lightyellow", "#EF3B36")) # nolint
        return(p)
    })
    pp <- patchwork::wrap_plots(p, ncol = 1)
    return(pp)
}
pdf("Figure6/Supplementary_Figure_9E_bp_density2.pdf", width = 1, height = 1)
Edensity(bp06_bcr, c("mu_freq_seq_r")) # nolint
dev.off()
##################################################################


##################################################################
## Figure 6D Fraction of Aber-MBC
# ---- 1) 频率矩阵 ----
meta_df <- as_tibble(bp@meta.data) |>
    select(orig.ident, anno2)

mat <- meta_df |>
    count(orig.ident, anno2, name = "n") |>
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
            rep(c("CAMR", "HC", "SAF"), times = c(18, 5, 6))[idx],
            levels = c("HC", "SAF", "CAMR")
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
mt1 <- mat |> filter(anno2 == "bp01.Pro-BC")
p2 <- make_box(mt1, ymax = 0.125)
ggsave("Figure6/Figure_6D_Aber-MBC_fraction.pdf", p2, width = 4, height = 6) # nolint
##################################################################


##################################################################
## Figure 6F Fraction in tissue
# 1) 取 meta，并确认每个样本(orig.ident)的分组信息
meta <- bp_tdata@meta.data %>%
    as_tibble(rownames = "cell") %>%
    select(orig.ident, label1, group)

# 2) 计算每个样本内 t_bp04_atMBC 的比例（样本内归一化）
#    若某些样本没有该亚群，补 0（更稳健，避免手动删除行）
all_samples <- meta %>% distinct(orig.ident, group)

prop_by_sample <- meta %>%
    count(orig.ident, label1, name = "n") %>%
    group_by(orig.ident) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    filter(label1 == "t_bp04_atMBC") %>%
    select(orig.ident, prop) %>%
    right_join(all_samples, by = "orig.ident") %>% # 补齐缺失样本
    mutate(
        value = if_else(is.na(prop), 0, prop),
        group = factor(group, levels = c("HC", "CAMR"))
    )

# 3) 画图
p <- ggplot(prop_by_sample, aes(x = group, y = value, color = group)) +
    geom_boxplot(outlier.size = 0, outlier.alpha = 0, width = 0.7) +
    stat_compare_means() +
    scale_color_manual(values = mycolor$group) +
    theme_classic() +
    theme(legend.position = "none")

ggsave("Figure6/Figure_6F_Fraction_in_tissue.pdf", p, width = 2.1, height = 3)
##################################################################


##################################################################
## Figure 6H Fraction in big data
# 统计并按 group 行归一化（等价于 table(..)/table(group)）
mat <- bp_pan@meta.data %>%
    as_tibble() %>%
    count(group, anno, name = "n") %>%
    group_by(group) %>%
    mutate(value = n / sum(n)) %>%
    ungroup() %>%
    mutate(
        group = factor(group, levels = c("N", "T", "I", "A")),
        # 反转 anno 的显示顺序（与原代码一致）
        anno  = factor(anno, levels = rev(levels(factor(anno))))
    )

p <- ggplot(mat, aes(x = value, y = anno, fill = group)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    theme_classic() +
    scale_fill_manual(values = pal_futurama(alpha = 0.8)(4)[c(3, 4, 1, 2)]) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(axis.line.y = element_blank())

ggsave("Figure6/Figure_6H_percent.pdf", p, width = 4.6, height = 3)
##################################################################


##################################################################
## Supplementary Figure 9D pvclust
# --- 1) 构建合并对象 ---
bp_sub <- subset(bp, idents = levels(bp)[c(1, 6)])
bp_sub$batch <- "batch2"
bp_sub$anno <- bp_sub$anno2

bp_pan_sub <- subset(bp_pan, downsample = 1000)
bp_pan_sub$batch <- "batch1"

merge_data <- merge(bp_sub, bp_pan_sub)

merge_data <- NormalizeData(
    merge_data,
    normalization.method = "LogNormalize", scale.factor = 10000
)
merge_data <- FindVariableFeatures(merge_data, selection.method = "vst", nfeatures = 2000) # nolint
all.genes <- rownames(merge_data)
merge_data <- ScaleData(merge_data, features = all.genes)
merge_data <- RunPCA(merge_data, features = VariableFeatures(object = merge_data)) # nolint
merge_data <- RunFastMNN(object.list = SplitObject(merge_data, split.by = "batch")) # nolint

# --- 2) 取 MNN 重构后的表达矩阵（data 槽），并对齐分组 ---
mat <- merge_data@assays$mnn.reconstructed@data # 通常为 dgCMatrix: genes x cells
anno <- merge_data$anno
stopifnot(ncol(mat) == length(anno))

# --- 4) 计算各 anno 的基因平均表达（一次性向量化） ---
# 思路 A（简单稳健）：对转置矩阵按组求和再除以计数
mat_t <- t(as.matrix(mat)) # cells x genes (dense 矩阵供 pvclust) # nolint
grp <- factor(anno)
sum_by_grp <- rowsum(mat_t, group = grp) # (levels(grp)) x genes 的“和” # nolint
n_by_grp <- as.vector(table(grp)) # 每组样本数
mtt <- sweep(sum_by_grp, 1, n_by_grp, `/`) # 求平均（逐行除以样本数）
mtt <- t(mtt) # genes x groups

# 行列名
rownames(mtt) <- genes_use
colnames(mtt) <- levels(grp)

# --- 5) pvclust 聚类 ---
# pvclust 需要普通 matrix
set.seed(1)
result <- pvclust(as.matrix(mtt),
    method.dist = "cor",
    method.hclust = "average",
    nboot = 1000,
    parallel = TRUE
)

pdf("Figure6/Supplementary_Figure_9D_pvclust.pdf", width = 7, height = 9)
plot(result)
pvrect(result, alpha = 0.95)
dev.off()
##################################################################


##################################################################
# 1) 基因向量 ——可只改这里
gene <- c(
    "FKBP1A", "FKBP5", "PPP3CA", "PPP3CC", "PPP3R1", "NFATC2", "NFATC3",
    "REL", "MAPK1", "MAPK14", "MAP3K2", "MAP3K5",
    "RASGRP1", "RB1CC1", "TAOK3", "CCDC88C",
    "PTK2B",
    "ATM", "UVRAG", "WAC", "PIP4K2A",
    "RPTOR", "VMP1", "VPS13A", "VPS13B", "VPS13C", "VPS13D",
    "BCL2", "ATG7", "NLRC3", "NLRC5",
    "DOCK8", "DOCK10", "RAPGEF6"
)

# 2) 取表达矩阵 & 仅保留存在的基因；准备分组因子
expr <- GetAssayData(bp_pan, assay = "RNA", slot = "data")
genes_use <- intersect(gene, rownames(expr))
if (length(genes_use) == 0) stop("None of the genes exist in the object.")

expr <- expr[genes_use, , drop = FALSE] # genes x cells
grp <- factor(bp_pan$anno) # 分组（行聚合用）
grp_lvls <- levels(grp)
n_by_grp <- as.vector(table(grp)) # 每组细胞数

# 3) 计算各组平均表达：一次性 rowsum + 除组内细胞数（向量化，替代 for+aggregate）
#    注意：pv 使用转置到 cells x genes 再 rowsum 按组求和 → 再除样本数 → 转回 genes x groups
expr_t <- t(as.matrix(expr)) # cells x genes（转 dense 以便运算）
sum_by_g <- rowsum(expr_t, group = grp) # groups x genes 的和
mean_mat <- t(sweep(sum_by_g, 1, n_by_grp, `/`))
colnames(mean_mat) <- grp_lvls
rownames(mean_mat) <- genes_use

# 4) 计算各组“阳性比例”（>0 的比例）：同样的向量化思路
pos_logical <- expr_t > 0 # cells x genes（逻辑矩阵）
pos_sum_by_g <- rowsum(pos_logical * 1, group = grp) # groups x genes 的“阳性个数”
pos_mat <- t(sweep(pos_sum_by_g, 1, n_by_grp, `/`))
colnames(pos_mat) <- grp_lvls
rownames(pos_mat) <- genes_use

# 5) 整成长表并对齐因子顺序
mer_mat <- list(
    pos  = as.data.frame(pos_mat) %>% rownames_to_column("gene") %>% pivot_longer(-gene, names_to = "Var2", values_to = "pos_value"), # nolint
    mean = as.data.frame(mean_mat) %>% rownames_to_column("gene") %>% pivot_longer(-gene, names_to = "Var2", values_to = "mean_value") # nolint
) %>%
    reduce(left_join, by = c("gene", "Var2")) %>%
    mutate(
        Var1 = gene,
        Var2 = factor(Var2, levels = rev(grp_lvls)) # 与原代码保持：反转分组显示顺序
    ) %>%
    select(Var1, Var2, pos_value, mean_value)

# 6) 作图
p <- ggplot(mer_mat, aes(y = Var2, x = Var1)) +
    geom_point(aes(size = pos_value, fill = mean_value), color = "black", shape = 22) + # nolint
    theme_bw() +
    scale_fill_gradientn(colors = c(rep("white", 10), colorRampPalette(c("white", "#ABAAAA"))(100))) + # nolint
    scale_size(range = c(0, 8)) +
    theme(
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave("Figure6/Figure_6M_exp.pdf", p, width = 9.5, height = 3.5)
##################################################################
