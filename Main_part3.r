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
    "ComplexHeatmap", "circlize", "stringr", "tibble"
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
load("input_file/color.RData")
##################################################################


##################################################################
## Figure 4A T cell subcluster Ro/e
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
    cellInfo.tb = nkt@meta.data, meta.cluster = nkt$anno2, loc = nkt$group, # nolint
    out.prefix = "Figure4/",
    pdf.width = 10, pdf.height = 8, verbose = 1
)

## OR value
or <- round(OR.B.list$OR.dist.mtx, 2)
colnames(or) <- colnames(OR.B.list$OR.dist.mtx)

p11 <- sscVis::plotMatrix.simple(or,
    out.prefix = "Figure4/Figure_4A_Roe_T_cluster.pdf",
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
ggsave("Figure4/Figure_4A_Roe_T_cluster.pdf", p11, width = 3, height = 7) # nolint
##################################################################


##################################################################
## Figure 4J DotPlot for expression of FKBP1A
nkt_sub <- subset(nkt, idents = levels(nkt)[1:13])
pdf("Figure4/Figure_4J_FKBP1A.pdf", width = 5.2, height = 1)
DotPlot(nkt_sub, features = "FKBP1A", cols = c("lightgrey", "red")) +
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
# ---- 并行与全局选项：仅设置一次 ----
plan("multicore", workers = 60)
options(future.globals.maxSize = 100000 * 1024^2)

# ---- 公共：制备 MSigDB 的 TERM2GENE（GO:BP + HALLMARK），只做一次复用 ----
get_term2gene <- function() {
    m_t2g <- msigdbr(species = "Homo sapiens") |>
        as.data.table()

    # 原作者有检查 GOBP 出现的类别，这里保留（无副作用）
    m_t2g[grepl("GOBP", gs_name)] |>
        pull(gs_cat) |>
        table() |>
        invisible()

    # 仅保留 GO:BP 与 HALLMARK，并选择需要的两列
    m_t2g |>
        dplyr::filter(gs_subcat == "GO:BP" | gs_cat == "H") |>
        dplyr::select(gs_name, entrez_gene)
}
TERM2GENE <- get_term2gene()

# ---- 公共：差异 + 富集 + 结果落盘的主函数（不改作图参数）----
run_de_enrich_and_plot <- function(
    obj, keep_annos, target_anno, # 数据与子集：保留的 anno3、目标类
    marker_out_path, profile_out_path, # 中间结果与富集表输出路径
    select_vec, filter_field = c("ID", "Description"), # 选择的通路列表及筛选字段
    plot_file, # 图输出路径
    palette_n, # 色板长度
    width, height, # ggsave 宽高（保持原值）
    blank_axis_text = c("all", "y") # 与原来两图一致：第一幅全空，第二幅只去掉 y 轴文字
    ) {
    filter_field <- match.arg(filter_field)
    blank_axis_text <- match.arg(blank_axis_text)

    # 1) 子集 + 打标签
    nkt_sub <- subset(obj, anno3 %in% keep_annos)
    nkt_sub$annox <- ifelse(nkt_sub$anno3 == target_anno, target_anno, "others")

    # 2) FindMarkers（参数保持不变）
    markers <- FindMarkers(
        nkt_sub,
        ident.1 = target_anno,
        ident.2 = "others",
        test.use = "wilcox",
        min.pct = 0.25,
        group.by = "annox",
        logfc.threshold = 0.25
    )
    markers$gene <- rownames(markers)

    # 与原代码保持一致的落盘方式（一个用 rds，一个用 txt）
    if (endsWith(marker_out_path, ".rds")) {
        write_rds(markers, marker_out_path)
        # 原脚本中紧接着 readRDS 再读回来，这里保持一致（虽然功能上可省略）
        markers <- readRDS(marker_out_path)
    } else {
        # 写成与原脚本完全一致的制表符文本
        utils::write.table(markers, marker_out_path, quote = FALSE, sep = "\t")
    }

    # 3) 取上调且显著
    p_row <- dplyr::filter(markers, avg_log2FC > 0, p_val_adj < 0.05)

    # 4) 基因 ID 转换
    eg <- bitr(p_row$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") # nolint
    p_row1 <- dplyr::left_join(p_row, eg, by = c("gene" = "SYMBOL")) |>
        dplyr::filter(!is.na(ENTREZID))

    # 5) 富集分析（TERM2GENE 复用）
    em <- enricher(p_row1$ENTREZID, TERM2GENE = TERM2GENE)
    x <- as.data.table(em)
    x$Description <- x$Description |> tolower()

    # 富集表落盘（保持 csv）
    utils::write.csv(x, profile_out_path, row.names = FALSE)

    # 6) 根据不同字段筛选通路（保持原行为）
    if (filter_field == "ID") {
        x_sub <- subset(x, ID %in% select_vec)
        # 与原代码一致：按 Description 的倒序设因子
        x_sub$Description <- factor(x_sub$Description, levels = rev(x_sub$Description)) # nolint
    } else {
        # Description 已转小写，传入的 select_vec 也需小写（外层已按你的两段代码分别给出）
        x_sub <- subset(x, Description %in% select_vec)
        x_sub$Description <- factor(x_sub$Description, levels = rev(x_sub$Description)) # nolint
    }

    # 7) 作图（严格保持两段代码的参数设置）
    p <- ggplot(x_sub, aes(x = -log(qvalue, 10), y = Description)) +
        geom_bar(aes(fill = Description), color = "black", stat = "identity") +
        scale_fill_manual(values = colorRampPalette(c("#F5E3DE", "#831A1F"))(palette_n)) + # nolint
        theme_classic() +
        {
            if (blank_axis_text == "all") {
                theme(
                    axis.title = element_blank(),
                    axis.text = element_blank(),
                    legend.position = "none"
                )
            } else {
                theme(
                    axis.title = element_blank(),
                    axis.text.y = element_blank(),
                    legend.position = "none"
                )
            }
        }

    ggsave(plot_file, p, width = width, height = height)
    invisible(p)
}

# --- 第一段（tn03）：保持原始选择与作图参数 ---
tn03_select_ID <- c(
    "GOBP_HISTONE_MODIFICATION", "GOBP_REGULATION_OF_GTPASE_ACTIVITY",
    "HALLMARK_UV_RESPONSE_DN", "GOBP_PEPTIDYL_SERINE_MODIFICATION",
    "HALLMARK_MITOTIC_SPINDLE", "GOBP_REGULATION_OF_SMALL_GTPASE_MEDIATED_SIGNAL_TRANSDUCTION", # nolint
    "GOBP_LEUKOCYTE_CELL_CELL_ADHESION", "GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY", # nolint
    "GOBP_PEPTIDYL_LYSINE_MODIFICATION", "GOBP_T_CELL_DIFFERENTIATION",
    "GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY", "GOBP_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY", # nolint
    "GOBP_LEUKOCYTE_PROLIFERATION", "GOBP_REGULATION_OF_RAS_PROTEIN_SIGNAL_TRANSDUCTION", # nolint
    "GOBP_ALPHA_BETA_T_CELL_ACTIVATION", "GOBP_REGULATION_OF_WNT_SIGNALING_PATHWAY", # nolint
    "GOBP_HISTONE_METHYLATION", "GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_2_PRODUCTION", # nolint
    "HALLMARK_IL2_STAT5_SIGNALING", "GOBP_FC_RECEPTOR_SIGNALING_PATHWAY",
    "GOBP_CALCINEURIN_MEDIATED_SIGNALING"
)

run_de_enrich_and_plot(
    obj = nkt,
    keep_annos = c("tn01", "tn02", "tn03", "tn04", "tn05"),
    target_anno = "tn03",
    marker_out_path = "temp/tn03_marker.rds", # 与原代码一致
    profile_out_path = "temp/tn03_bp_hallmarker_profile.csv",
    select_vec = tn03_select_ID,
    filter_field = "ID", # 第一段按 ID 过滤（保持不变）
    plot_file = "Figure4/Figure_4I_tn03_go_bp.pdf",
    palette_n = 21, # 与原色板长度一致
    width = 4, height = 6,
    blank_axis_text = "all" # 第一段 axis.text 全空
)

# --- 第二段（tn10）：保持原始选择与作图参数 ---
tn10_select_desc_lower <- c(
    "gobp_histone_modification", "gobp_peptidyl_serine_modification",
    "gobp_regulation_of_small_gtpase_mediated_signal_transduction",
    "gobp_immune_response_regulating_signaling_pathway",
    "gobp_regulation_of_gtpase_activity", "hallmark_uv_response_dn",
    "gobp_protein_dephosphorylation", "gobp_antigen_receptor_mediated_signaling_pathway", # nolint
    "gobp_regulation_of_cell_matrix_adhesion", "gobp_t_cell_receptor_signaling_pathway", # nolint
    "gobp_regulation_of_ras_protein_signal_transduction", "gobp_peptidyl_lysine_modification", # nolint
    "hallmark_mitotic_spindle", "gobp_leukocyte_proliferation",
    "gobp_peptidyl_threonine_modification", "gobp_calcineurin_mediated_signaling", # nolint
    "gobp_stress_activated_protein_kinase_signaling_cascade",
    "gobp_regulation_of_wnt_signaling_pathway", "gobp_alpha_beta_t_cell_activation", # nolint
    "hallmark_pi3k_akt_mtor_signaling", "gobp_regulation_of_jnk_cascade",
    "gobp_regulation_of_autophagy", "gobp_positive_regulation_of_cytokine_production", # nolint
    "gobp_regulation_of_extrinsic_apoptotic_signaling_pathway"
)

run_de_enrich_and_plot(
    obj = nkt,
    keep_annos = c("tn06", "tn07", "tn09", "tn10"),
    target_anno = "tn10",
    marker_out_path = "temp/tn10_marker.txt", # 与原代码一致（写 txt）
    profile_out_path = "temp/tn10_bp_hallmarker_profile.csv",
    select_vec = tn10_select_desc_lower,
    filter_field = "Description", # 第二段按小写 Description
    plot_file = "Figure4/Figure_4I_tn10_go_bp.pdf",
    palette_n = 24, # 与原色板长度一致
    width = 1.5, height = 7.5,
    blank_axis_text = "y" # 第二段仅去掉 y 轴文字
)
##################################################################


##################################################################
## Figure 4B Fraction of ANK3+CD4+ Tm and ZEB2hi T
# ---- 1) 频率矩阵 ----
meta_df <- as_tibble(nkt@meta.data) |>
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
mt1 <- mat |> filter(anno2 == "tn03_ANK3+CD4+ Tm")
p2 <- make_box(mt1, ymax = 0.04)

mt2 <- mat |> filter(anno2 == "tn10_ZEB2+ T")
p3 <- make_box(mt2, ymax = 0.03)

final <- p2 + plot_spacer() + p3 + plot_layout(nrow = 1, widths = c(1, 0.1, 1, 0.1, 1)) # nolint
ggsave("Figure4/Figure_4B_ANK3_ZEB2_T_fraction.pdf", final, width = 8, height = 6) # nolint
##################################################################


##################################################################
## Figure 4C correlation of ANK3+CD4+ Tm and ZEB2hi T
# 1) 计算每个样本内各 anno2 的比例
mat <- nkt@meta.data |>
    as.data.frame() |>
    count(orig.ident, anno2, name = "n") |>
    group_by(orig.ident) |>
    mutate(prop = n / sum(n)) |>
    ungroup() |>
    filter(anno2 %in% c("tn03_ANK3+CD4+ Tm", "tn10_ZEB2+ T")) |>
    select(orig.ident, anno2, prop) |>
    pivot_wider(names_from = anno2, values_from = prop, values_fill = 0) |>
    # 添加分组：HC 5 个，SAF 6 个，CAMR 18 个
    mutate(group = factor(
        c(rep("CAMR", 18), rep("HC", 5), rep("SAF", 6)),
        levels = c("HC", "SAF", "CAMR")
    ))

# 2) 绘图（保持原参数）
p1 <- ggplot(mat, aes(x = `tn03_ANK3+CD4+ Tm`, y = `tn10_ZEB2+ T`)) +
    geom_point(aes(color = group), size = 1) +
    scale_color_manual(values = mycolor$group) +
    stat_cor() +
    geom_smooth(method = "lm", linewidth = 0.2, color = "black") +
    theme_classic() +
    theme(
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none"
    )
##################################################################


##################################################################
## Figure 4D Fraction of ANK3+CD4+ Tm in tissue data
# ---- 1) 频率矩阵 ----
meta_df <- as_tibble(nkt_tdata@meta.data) |>
    select(orig.ident, label1)

mat <- meta_df |>
    count(orig.ident, label1, name = "n") |>
    group_by(orig.ident) |>
    mutate(value = n / sum(n)) |>
    ungroup()

all_labels <- sort(unique(mat$label1))
mat <- mat %>%
    group_by(orig.ident) %>%
    tidyr::complete(label1 = all_labels, fill = list(n = 0, value = 0)) %>%
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

mt1 <- mat |> filter(label1 == "tn09_ANK3+ Tm")
p2 <- make_box(mt1, ymax = 0.08)
ggsave("Figure4/Figure_4C_ANK3_T_tissue_fraction.pdf", p2, width = 3.5, height = 6) # nolint
##################################################################


##################################################################
## Figure 4H Heatmap
# 基因列表（原样）
genes <- c(
    "FKBP5", "LPCAT1", "PPP3CA", "PPP3CC", "PPP3R1", "NFATC2", "NFATC3", # calcineurin # nolint
    "REL", "MAPK1", "MAPK14", "MAP3K2", "MAP3K5", # stress
    "RASGRP1", "RB1CC1", "TAOK3", "TNIK", "CCDC88C", # JNK
    "PLCB1", "PTK2B", # autophagy   # nolint
    "ATM", "UVRAG", "WAC", "IFI16", "PIP4K2A", # autophagy   # nolint
    "RPTOR", "VMP1", "VPS13A", "VPS13B", "VPS13C", "VPS13D", # autophagy   # nolint
    "BCL2", "ATG7", "NLRC3", "NLRC5", # autophagy   # nolint
    "DOCK8", "DOCK10", "RAPGEF6", "PTPRC", "CD247"
)

# 1) 子集（保持：取按字典序排序后前 13 个 anno3）
nkt_sub <- subset(nkt, anno3 %in% head(sort(unique(nkt$anno3)), 13))

# 2) 计算各 anno2 的平均表达（等价于原来 aggregate(..., mean)）
#    AverageExpression 默认使用 slot="data"，与原代码一致
ae <- suppressMessages(
    AverageExpression(
        nkt_sub,
        assays = "RNA",
        features = genes,
        group.by = "anno2",
        slot = "data",
        return.seurat = FALSE
    )
)$RNA # 基因 x 组 的矩阵

# 按出现顺序固定列顺序，避免因因子/字母序造成的列重排
group_levels <- unique(as.character(nkt_sub$anno2))
ae <- ae[, intersect(group_levels, colnames(ae)), drop = FALSE]

# 行顺序按输入 genes 固定（若有缺失基因则自动剔除）
ae <- ae[intersect(genes, rownames(ae)), , drop = FALSE]

# 3) 按基因对各组做 Z-score（与原来 scale(t(matr)) 再 t() 等价）
matr <- t(scale(t(ae)))
matr[is.na(matr)] <- 0

# 4) 颜色映射保持不变
col1 <- colorRamp2(c(-1, 0, 1), colors = viridis(3, option = "A"))

# 5) 绘图（所有 Heatmap 参数与原代码完全一致）
pdf("Figure4/Figure_4H_Heatmap.pdf", width = 4, height = 8)
Heatmap(
    matr, # 这里已是 基因 x 组（等价于原来的 t(matr)）
    name = "Enrichment\nscore",
    cluster_columns = TRUE,
    cluster_rows = FALSE,
    row_title_rot = 0,
    row_names_gp = gpar(fontsize = 10),
    row_title_gp = gpar(fontsize = 11),
    column_names_gp = gpar(fontsize = 10),
    column_title_gp = gpar(fontsize = 11),
    row_names_side = "left",
    row_title_side = "right",
    column_km = 2,
    row_km = 1,
    col = col1
)
dev.off()
##################################################################


##################################################################
## Figure 4G Fraction of T cells in big blood data
# 1) 计算各 group 内的比例
anno_levels <- levels(t_pan$anno)
meta_df <- as_tibble(t_pan@meta.data) %>%
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

ggsave("Figure4/Figure_4G_pan_t_percent.pdf", p, width = 4.6, height = 3)
##################################################################


##################################################################
## Supplementary Figure 5D pvcluster tree
# ---- 1) 子集 ----
nkt_sub <- subset(nkt, idents = levels(nkt)[c(3, 9)])
nkt_sub$batch <- "batch2"
nkt_sub$anno <- nkt_sub$anno2

# ---- 2) 下采样并合并（保持参数不变）----
t_pan@active.ident <- factor(t_pan$anno)
t_pan_sub <- subset(t_pan, downsample = 1000)
t_pan_sub$batch <- "batch1"

merge_data <- merge(nkt_sub, t_pan_sub)

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

pdf("Figure4/Supplementary_Figure_5D_pvcluster_tree.pdf", width = 7, height = 9)
plot(result)
pvrect(result, alpha = 0.95)
dev.off()
##################################################################


##################################################################
## Supplementary Figure 5E autograpy gene expression
# ---- 1) 取表达矩阵（RNA@data）与基因子集 ----
genes <- c(
    "FKBP5", "LPCAT1", "PPP3CA", "PPP3CC", "PPP3R1", "NFATC2", "NFATC3", # calcineurin # nolint
    "REL", "MAPK1", "MAPK14", "MAP3K2", "MAP3K5", # stress
    "RASGRP1", "RB1CC1", "TAOK3", "TNIK", "CCDC88C", # JNK
    "PLCB1", "PTK2B", # autophagy
    "UVRAG", "WAC", "IFI16", "PIP4K2A", # autophagy
    "RPTOR", "VMP1", "VPS13A", "VPS13B", "VPS13C", "VPS13D", # autophagy
    "BCL2", "ATG7", "NLRC3", "NLRC5", # autophagy
    "DOCK8", "DOCK10", "RAPGEF6"
)

expr <- GetAssayData(t_pan, assay = "RNA", slot = "data") # dgCMatrix
expr <- expr[genes, , drop = FALSE]

# ---- 2) 构造分组因子与顺序（与原顺序完全一致）----
anno_levels <- levels(t_pan$anno)
g <- factor(t_pan$anno, levels = anno_levels)

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

ggsave("Figure4/Supplementary_Figure_5E_gene_exp.pdf", p, width = 10, height = 3.5) # nolint
##################################################################


##################################################################
## Figure 4M DEG
# —— 1) 取子集（Seurat 对象用 subset 更稳）——
nkt_sub <- subset(nkt, anno3 %in% "tn10")

# —— 2) 抽象一个通用的 FindMarkers 封装 ——
# contrast_name 用来在列名上打标签（如 "HC" / "SAF"）
run_markers <- function(obj, ident1, ident2, group_by = "group",
                        test.use = "wilcox",
                        min.pct = 0, logfc.threshold = 0,
                        contrast_name = ident2) {
    res <- FindMarkers(
        obj,
        ident.1 = ident1,
        ident.2 = ident2,
        group.by = group_by,
        test.use = test.use,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold
    ) %>%
        rownames_to_column("gene") %>%
        # 给关键列加上后缀以便合并不冲突
        rename(
            !!str_glue("avg_log2FC_{contrast_name}") := avg_log2FC,
            !!str_glue("pct.1_{contrast_name}") := pct.1,
            !!str_glue("pct.2_{contrast_name}") := pct.2,
            !!str_glue("p_val_{contrast_name}") := p_val,
            !!str_glue("p_val_adj_{contrast_name}") := p_val_adj
        )
    res
}

# —— 3) 首轮（筛选基因集合）：min.pct=0.25, logfc=0.25 ——
contrasts <- c("HC", "SAF")

genes_keep <-
    map(contrasts, ~ run_markers(
        obj = nkt_sub, ident1 = "CAMR", ident2 = .x,
        group_by = "group", test.use = "wilcox",
        min.pct = 0.25, logfc.threshold = 0.25,
        contrast_name = .x
    )) %>%
    map(~ pull(.x, gene)) %>%
    reduce(union)

# —— 4) 第二轮（无阈值统计）：min.pct=0, logfc=0，并合并到一张表 ——
mark_wide <-
    map(contrasts, ~ run_markers(
        obj = nkt_sub, ident1 = "CAMR", ident2 = .x,
        group_by = "group", test.use = "wilcox",
        min.pct = 0, logfc.threshold = 0,
        contrast_name = .x
    )) %>%
    # 先把每个对比都过滤到 genes_keep，再按 gene 内连接到宽表
    map(~ filter(.x, gene %in% genes_keep)) %>%
    reduce(~ inner_join(.x, .y, by = "gene"))

# —— 5) 复现你原来的派生列（si / la）与过滤 ——
genes_label <- c(
    "FKBP5", "JUND", "FOS", "THEMIS", "JUN", "DUSP1",
    "GZMB", "CD69", "NFKB1", "IL7R", "STAT1"
)

mark <- mark_wide %>%
    mutate(
        # 与原逻辑一致：若两边 log2FC 都 < 0，用 SAF 的 pct.2，否则用 HC 的 pct.1
        si = if_else(
            .data$avg_log2FC_HC < 0 & .data$avg_log2FC_SAF < 0,
            .data$`pct.2_SAF`,
            .data$`pct.1_HC`
        ),
        la = if_else(gene %in% genes_label, gene, NA_character_)
    ) %>%
    filter(gene != "FP671120.4")

# —— 6) 作图 ——
p <- ggplot(mark, aes(x = avg_log2FC_SAF, y = avg_log2FC_HC)) +
    geom_vline(xintercept = 0, linewidth = 0.8, color = "black", linetype = "dashed") + # nolint
    geom_hline(yintercept = 0, linewidth = 0.8, color = "black", linetype = "dashed") + # nolint
    geom_point(aes(fill = avg_log2FC_SAF, size = si), shape = 21) +
    scale_fill_gradientn(colors = colorRampPalette(c("#594C99", "white", "#752C2A"))(100)) + # nolint
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        panel.border = element_rect(size = 1.5, color = "black")
    )

ggsave("Figure4/Figure_4M_deg_vol.pdf", p, width = 4.5, height = 4.5)
##################################################################
