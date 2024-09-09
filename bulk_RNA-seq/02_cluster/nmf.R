# environment
workdir = 'D://_SFT/'
setwd(workdir)

# packages
library(tidyverse)
library(NMF)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(GSEABase)
library(GSVA)
library(patchwork)
library(ggrepel)

# function
source("../functions.R")

# step 1
sft_cluster_mat = read_tsv("normalized_scale.data.tsv") %>% 
    column_to_rownames("gene_id") %>%
    as.matrix()
dim(sft_cluster_mat)

if(nrow(sft_cluster_mat) == nrow(filter(sample_rna_clean, is.na(Status)))){
    message("Sample load successful")
}


## 用NMF做聚类，选择最佳k值
runNMF1 = function(exp_matrix, gene_num = seq(1000, 2000, 100), rank = 3:6) {
    exp_by_mad = exp_matrix[rev(order(apply(exp_matrix, 1, mad))),]
    
    res = lapply(gene_num, function(gene) {
        nmf_res = nmf(exp_by_mad[1:gene,], rank = rank, seed = 12345, nrun = 20, .opt = 'vp100')
        coph = nmf_res$measures$cophenetic
        names(coph) = paste0("rank", rank)
        
        png(paste0("gene", gene, "rank.png"))
        par(cex.axis = 1.5)
        plot(rank, coph, xlab = "", ylab = "", type = "b", col = "red", lwd = 4, xaxt = "n")
        axis(side = 1, at = 1:length(rank))
        title(xlab = "number of clusters", ylab = "Cophenetic coefficient", cex.lab = 1.5)
        dev.off()
        
        return(coph)
    })
    
    res_matrix = do.call(rbind, res)
    colnames(res_matrix) = paste0("rank", rank)
    rownames(res_matrix) = paste0("gene", gene_num)
    write.csv(res_matrix, file = "coph_rank_gene.csv")
}

runNMF1(sft_cluster_mat)

coph_matrix = read.csv("coph_rank_gene.csv", row.names = 1)
coph_data = as.data.frame(as.table(as.matrix(coph_matrix)))

# visualization 
colnames(coph_data) = c("Gene_Number", "Rank", "Cophenetic")
ggplot(coph_data, aes(x = Rank, y = Cophenetic, group = Gene_Number, color = Gene_Number)) +
    geom_line() +
    geom_point() +
    labs(title = "Cophenetic Coefficients for Different k Values",
         x = "Number of Clusters (k)",
         y = "Cophenetic Coefficient") +
    theme_minimal()

# step 2
## 用NMF做聚类，自动选择最佳聚类结果
runNMF2 = function(exp_matrix, k = 4, gene_num = seq(1000, 2000, 100)) {
    exp_by_mad = exp_matrix[rev(order(apply(exp_matrix, 1, mad))),][1:max(gene_num),]
    
    res = lapply(gene_num, function(num_genes) {
        selected_genes = rownames(exp_by_mad)[1:num_genes]
        nmf_res = nmf(exp_by_mad[1:num_genes,], rank = k, seed = 12345, nrun = 20, .opt = 'vp100')
        list(coph = cophcor(nmf_res), cls = as.numeric(predict(nmf_res)), genes = selected_genes)
    })
    
    names(res) = paste0("gene", gene_num)
    return(res)
}

nmf_results = runNMF2(sft_cluster_mat)
cluster_results = sapply(nmf_results, `[[`, 'cls')
rownames(cluster_results) = colnames(sft_cluster_mat)
write.table(cluster_results, file = "cluster_results.txt", sep = "\t")


# step 3
## 留一法做重复
leave_one_out = function(index, exp_matrix) {
    tmp_matrix = exp_matrix[, -index]
    nmf_res = runNMF2(tmp_matrix)
    cls_matrix = sapply(nmf_res, `[[`, 'cls')
    gene_list = unlist(lapply(nmf_res, `[[`, 'genes'))
    rownames(cls_matrix) = colnames(tmp_matrix)
    coph_values = sapply(nmf_res, `[[`, 'coph')
    
    best_cls = cls_matrix[, which.max(coph_values)]
    additional_cls = setNames("N", setdiff(colnames(exp_matrix), names(best_cls)))
    final_cls = c(best_cls, additional_cls)[colnames(exp_matrix)]
    
    return(list(cls = final_cls, genes = gene_list))
}

final_nmf_results = sapply(1:ncol(sft_cluster_mat), 
                           leave_one_out, exp_matrix = sft_cluster_mat, simplify = FALSE)

# 提取样本的聚类结果
final_cls_matrix = do.call(cbind, lapply(final_nmf_results, `[[`, 'cls'))
dim(final_cls_matrix)

# 保存最终的聚类结果矩阵
write.table(final_cls_matrix, file = "final_cls_matrix.txt", sep = "\t", quote = FALSE)

## 生成一个方阵，用于统计聚类结果
count_clusters = function(cluster_vector) {
    stat_matrix = matrix(0, length(cluster_vector), length(cluster_vector), 
                         dimnames = list(names(cluster_vector), names(cluster_vector)))
    unique_clusters = unique(cluster_vector)
    
    for (cluster in unique_clusters) {
        cluster_names = names(cluster_vector)[which(cluster_vector == cluster)]
        stat_matrix[cluster_names, cluster_names] = 1
    }
    
    return(stat_matrix)
}

final_stat_matrix = matrix(0, nrow(final_cls_matrix), nrow(final_cls_matrix))
for (i in 1:ncol(final_cls_matrix)) {
    tmp_matrix = count_clusters(final_cls_matrix[, i])
    final_stat_matrix = final_stat_matrix + tmp_matrix
}
write.table(final_stat_matrix, file = "final_stat_matrix.txt", sep = "\t", quote = FALSE)

# ## 统计聚类结果
final_stat_matrix = read.table("final_stat_matrix.txt", sep = "\t", header = TRUE)

sft_cluster_info = sft_cluster_info %>%
    mutate(Alias = recode(cluster, `1` = "Inflamed",
                          `2` = "Classical",
                          `3` = "Neural-like",
                          `4` = "Migratory")) %>%
    write_tsv("sample_anno_cluster.tsv")

# 一致性矩阵可视化
nmf_plot_idx = sft_cluster_info %>% filter(Location == "CNS", is.na(Status)) %>%
    mutate(idx = paste(Batch, PatientID, sep = "_")) %>% pull(idx)

nmf_heat_anno = sft_cluster_info %>% 
    filter(Location == "CNS", is.na(Status)) %>% 
    mutate(idx = paste(Batch, PatientID, sep = "_")) %>%
    mutate(Cluster = recode(cluster, `1` = "Inflamed",
                            `2` = "Classical",
                            `3` = "Neural-like",
                            `4` = "Migratory"),
           Cluster = factor(Cluster, levels = c("Classical","Neural-like","Inflamed","Migratory")),
           Grade = paste("WHO", Grade)) %>%
    column_to_rownames("idx") %>%
    dplyr::select(Grade, Cluster)


pdf("NMF_clustering_legend.pdf", width = 5.5, height = 5)
cluster_groups = pheatmap::pheatmap(
    final_stat_matrix[nmf_plot_idx, nmf_plot_idx], 
    scale = "none", cluster_rows = TRUE, show_rownames = F, 
    show_colnames = F, cluster_cols = TRUE, 
    annotation_col = nmf_heat_anno,
    annotation_row = nmf_heat_anno,
    annotation_colors = list(
        Grade = custom_colors$Grade,
        Cluster = custom_colors$Cluster
    ),
    treeheight_row = 12, treeheight_col = 12,
    fontsize_row = 3, border = NA, 
    clustering_distance_rows = "correlation", 
    breaks = seq(0, 65, length = 100), 
    color = colorRampPalette(c("grey90", rb_two[2]))(100))
dev.off()

# step4
# genes and features
sft_cluster_mad2k = sft_cluster_mat[rev(order(apply(sft_cluster_mat, 1, mad))),][1:2000,] %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>%
    separate(col = ID, into = c("ENSEMBL", "SYMBOL")) %>%
    dplyr::select(-SYMBOL) %>%
    column_to_rownames("ENSEMBL") %>%
    as.matrix() 

# GSVA NES score matrix
sft_seurat_nmf_nes = bind_rows(
    as.data.frame(OneStepGSVA(cat = "H", matrix = sft_cluster_mad2k)),
    as.data.frame(OneStepGSVA(cat = "C2", matrix = sft_cluster_mad2k)),
    as.data.frame(OneStepGSVA(cat = "C5", subcat = "MF", matrix = sft_cluster_mad2k)),
    as.data.frame(OneStepGSVA(cat = "C5", subcat = "BP", matrix = sft_cluster_mad2k)),
    as.data.frame(OneStepGSVA(cat = "C5", subcat = "CC", matrix = sft_cluster_mad2k)),
) 

# NES score to longer
sft_seurat_nmf_nes_longer = sft_seurat_nmf_nes %>%
    rownames_to_column("terms")  %>%
    pivot_longer(cols = contains("batch"), names_to = "sample", values_to = "NES") %>%
    separate(sample, sep = "_", into = c("Batch", "PatientID")) %>%
    left_join(rownames_to_column(sft_cluster_mad_anno, "PatientID"))

# significant NES and terms
sft_seurat_nmf_nes_tests = map_df(unique(sft_seurat_nmf_nes_longer$terms), ~ NFStest(.x, sft_seurat_nmf_nes_longer)) %>%
    pivot_longer(cols = contains("p"), names_to = "p_type", values_to = "p_value") %>%
    group_by(p_type) %>%
    mutate(p_adj = p.adjust(p_value, method = "bonferroni")) 

sft_seurat_nmf_nes_tests %>% filter(p_adj < 0.05) %>%
    mutate(p_type = recode(p_type,
                           "c1p" = "Cluster1vsOther",
                           "c2p" = "Cluster2vsOther",
                           "c3p" = "Cluster3vsOther",
                           "c4p" = "Cluster4vsOther",
                           "g1p" = "Grade1vsOther",
                           "g2p" = "Grade1vsOther",
                           "g3p" = "Grade1vsOther")) %>%
    write_csv("cluster_gsva_pvalue005.csv")

# top n significant NES and terms
SigNEStermsHeatmap(20)

# step5
# DE genes among clusters
library(edgeR)
library(limma)

cluster_countdata = sft_cluster_mat
dim(cluster_countdata)

cluster_coldata = cutree(cluster_groups$tree_row, k = 4) %>%
    as.data.frame() %>%
    set_names("cluster") 

cluster_countdata = cluster_countdata[rowSums(cluster_countdata) > 0,]
cluster_group = factor(cluster_coldata$cluster)

dge = DGEList(counts = cluster_countdata, group = cluster_group) %>%
    calcNormFactors()

comparisons = list(
    "cluster1_vs_rest" = 1,
    "cluster2_vs_rest" = 2,
    "cluster3_vs_rest" = 3,
    "cluster4_vs_rest" = 4
)

compare_clusters = function(cluster_idx, rest_idx, dge, comparison_name) {
    # group = factor(c(rep("cluster", length(cluster_idx)), rep("rest", length(rest_idx))))
    group = factor(c(rep("rest", length(rest_idx)), rep("cluster", length(cluster_idx))), levels = c("rest", "cluster"))
    design = model.matrix(~ group)
    
    # v = voom(dge[, c(cluster_idx, rest_idx)], design)
    v = voom(dge[, c(rest_idx, cluster_idx)], design)
    
    fit = lmFit(v, design)
    fit = eBayes(fit)
    all_genes_results = topTable(fit, coef = 2, number = Inf, sort.by = "none")
    all_genes_results$comparison = comparison_name
    return(all_genes_results)
}

cluster_res = lapply(names(comparisons), function(comparison_name) {
    cluster_num = comparisons[[comparison_name]]
    cluster_idx = which(cluster_group == cluster_num)
    rest_idx = setdiff(seq_along(cluster_group), cluster_idx)
    compare_clusters(cluster_idx, rest_idx, dge, comparison_name)
}) %>% do.call(rbind, .)


# volcano
cluster_volcano(cluster_res, comp = "cluster1_vs_rest") + 
    cluster_volcano(cluster_res, comp = "cluster2_vs_rest") +
    cluster_volcano(cluster_res, comp = "cluster3_vs_rest") + 
    cluster_volcano(cluster_res, comp = "cluster4_vs_rest")

# heatmap of clusters
nfm_deg_idx = cluster_res %>% 
    rownames_to_column("geneID") %>%
    separate(geneID, sep = "-", into = c("gene_id", "gene_name")) %>%
    mutate(gene_name = case_when(
        comparison == "cluster2_vs_rest" ~ str_replace(gene_name, "1$", ""),
        comparison == "cluster3_vs_rest" ~ str_replace(gene_name, "2$", ""),
        comparison == "cluster4_vs_rest" ~ str_replace(gene_name, "3$", ""),
        TRUE ~ gene_name)) %>% 
    filter(adj.P.Val < 0.05, logFC>log2(1.5)) %>%
    group_by(comparison) %>%
    arrange(desc(logFC), adj.P.Val) %>%
    dplyr::slice(1:20) %>%
    pull(gene_name) %>% unique()


nfm_deg_anno = sft_cluster_info %>% 
    filter(Location == "CNS", is.na(Status)) %>% 
    mutate(Cluster = recode(cluster, `1` = "Inflamed",
                            `2` = "Classical",
                            `3` = "Neural-like",
                            `4` = "Migratory"),
           Cluster = factor(Cluster, levels = c("Neural-like","Classical","Inflamed","Migratory")),
           Grade = paste("WHO", Grade)) %>%
    column_to_rownames("PatientID") %>%
    dplyr::select(Grade, Cluster)

pdf("cluster_DEG_heatmap.pdf", height = 9, width = 7)
sft_tpm %>%
    filter(geneName %in% nfm_deg_idx) %>%
    .[!duplicated(.$geneName),] %>%
    column_to_rownames("geneName") %>% 
    dplyr::select(any_of(filter(sft_cluster_info, is.na(Status), Location == "CNS")$PatientID)) %>%
    .[,rownames(nfm_deg_anno[order(nfm_deg_anno$Cluster),])] %>%
    pheatmap::pheatmap(scale = "row", cluster_cols = F,
                       show_colnames = F,
                       annotation_col = nfm_deg_anno,
                       annotation_colors = list(
                           Grade = custom_colors$Grade,
                           Cluster = custom_colors$Cluster
                       ),
                       fontsize_row = 7, border_color = NA,
                       treeheight_row = 12, treeheight_col = 12,
                       clustering_distance_rows = "correlation",
                       breaks = seq(-1.5,1.5, length = 100), 
                       color = colorRampPalette(c("blue","white","red"))(100))
dev.off()

# out
cluster_res %>% 
    rownames_to_column("geneID") %>%
    separate(geneID, sep = "-", into = c("gene_id", "gene_name")) %>%
    filter(abs(logFC) > log2(1.5), adj.P.Val < 0.05)  %>%
    mutate(gene_name = case_when(
        comparison == "cluster2_vs_rest" ~ str_replace(gene_name, "1$", ""),
        comparison == "cluster3_vs_rest" ~ str_replace(gene_name, "2$", ""),
        comparison == "cluster4_vs_rest" ~ str_replace(gene_name, "3$", ""),
        TRUE ~ gene_name)) %>% 
    write_csv("Sample_cluster_DEG.csv")
