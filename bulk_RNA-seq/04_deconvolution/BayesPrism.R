# environment
workdir = "/media/data5/hxy/SFT/RNAseq/3_BayesPrism/"
setwd(workdir)

# packages
library(tidyverse)
library(Seurat)
library(BayesPrism)

# input data
message("*** Loading RDS ***")
sc_data = readRDS("SFT_res06_final2.rds")
cell_expr_matrix = as.data.frame(Matrix::t(sc_data@assays$RNA@counts))
cell_type_label = sc_data@meta.data[["cell_type2"]]
cell_state_label = sc_data@meta.data[["cell_state"]]

message("*** Loading Expr ***")
bulk_expr_matrix = read_tsv("../1_BatchRecor/normalized_scale.data.tsv") %>% 
    column_to_rownames("gene_id") %>%
    as.matrix()

colnames(bulk_expr_matrix) = str_split(colnames(bulk_expr_matrix), "_", simplify = T)[,2]
rownames(bulk_expr_matrix) = str_split(rownames(bulk_expr_matrix), "-", simplify = T)[,2]

bulk_expr_matrix = bulk_expr_matrix %>%
    as.data.frame() %>%
    mutate(gene_id = rownames(bulk_expr_matrix))

message("*** Loading Sample info ***")
#sample_rna_clean = read_tsv("../RNA-seq_all_clean.tsv") %>%
#   mutate(Grade = as.character(Grade))
sample_rna_clean = read_tsv("../1_BatchRecor/cluster/sample_cluster.tsv") %>%
   mutate(Grade = as.character(Grade),
        cluster = as.character(cluster))

# matrix overlap checking
table(bulk_expr_matrix$gene_id %in% colnames(cell_expr_matrix))
table(colnames(cell_expr_matrix) %in% bulk_expr_matrix$gene_id)

# prepare matrix
keep_gene = intersect(bulk_expr_matrix$gene_id, colnames(cell_expr_matrix))

bulk_expr_matrix_keep = bulk_expr_matrix %>%
    dplyr::filter(gene_id %in% keep_gene) %>%
    .[!duplicated(.$gene_id),] %>%
    dplyr::select(-gene_id) %>%
    as.matrix %>% t()

cell_expr_matrix_keep = cell_expr_matrix %>% dplyr::select(all_of(keep_gene))

if (ncol(bulk_expr_matrix_keep) != ncol(cell_expr_matrix_keep)){
    stop("Error: mismatch genes in cell matrix and bulk matrix")
} else if(nrow(cell_expr_matrix_keep) != length(cell_type_label)){
    stop("Error: mismatch cell type in cell matrix and cell type lable")
} else {
    message("*** Precessed ***")
}

# quality control
# correlation of cell types
message("*** Plot correlation of cell type ***")
plot.cor.phi(
    input = cell_expr_matrix_keep,
    input.labels = cell_type_label,
    title="cell type correlation",
    pdf.prefix="sft.cor.cs_type", 
    cexRow=0.2, cexCol=0.2,
    cex.axis = 1.5,
    cex.labels = 1.5,
    margins=c(2,2)
)

message("*** Plot correlation of cell states ***")
plot.cor.phi(
    input = cell_expr_matrix_keep,
    input.labels = cell_state_label,
    title="cell type correlation",
    pdf.prefix="sft.cor.cs_state", 
    cexRow=0.2, cexCol=0.2,
    cex.axis = 1.5,
    cex.labels = 1.5,
    margins=c(2,2)
)

# outlier of cell epxr matrix
message("*** Plot cell type outliers of scRNA-seq ***")
sc_stat = plot.scRNA.outlier(   
    input = cell_expr_matrix_keep, 
    cell.type.labels = cell_type_label,
    species="hs",
    return.raw=TRUE, 
    pdf.prefix="sft.sc.stat_type"
)

message("*** Plot cell state outliers of scRNA-seq ***")
sc_stat = plot.scRNA.outlier(   
    input = cell_expr_matrix_keep, 
    cell.type.labels = cell_state_label,
    species="hs",
    return.raw=TRUE, 
    pdf.prefix="sft.sc.stat_state"
)

# outlier of bulk 
message("*** plot outliers of bulk rna-seq ***")
bk_stat = plot.bulk.outlier(
    bulk.input = bulk_expr_matrix_keep, 
    sc.input = cell_expr_matrix_keep, 
    cell.type.labels = cell_type_label,
    species="hs",
    return.raw=TRUE,
    pdf.prefix="sft.bk.stat_type"
)
bk_stat = plot.bulk.outlier(
    bulk.input = bulk_expr_matrix_keep, 
    sc.input = cell_expr_matrix_keep, 
    cell.type.labels = cell_state_label,
    species="hs",
    return.raw=TRUE,
    pdf.prefix="sft.bk.stat_state"
)


# filter outlier
message("*** filter scRNA-seq ***")
cell_expr_filtered = cleanup.genes(
    input = cell_expr_matrix_keep,
    input.type = "count.matrix",
    species = "hs", 
    gene.group = c( "Rb","Mrp","other_Rb","chrM", "MALAT1","chrX","chrY","hb","act"),
    exp.cells = 5
)

# correlation gene types of ciell and bulk
message("*** Plot gene correlation ***")
plot.bulk.vs.sc(
    sc.input = cell_expr_filtered,
    bulk.input = bulk_expr_matrix_keep,
    pdf.prefix="sft.bk.vs.sc"
)
 
# selecting signature genes for BayesPrism
# select protein coding genes
message("*** filter protein-coding genes ***")
cell_expr_filtered_pc = select.gene.type(
    cell_expr_filtered,
    gene.type = "protein_coding"
)

# BayesPrism
message("*** Start BayesPrism ***")
myPrism = new.prism(
  reference = cell_expr_filtered_pc, 
  mixture = bulk_expr_matrix_keep,
  input.type="count.matrix", 
  cell.type.labels = cell_type_label, 
  cell.state.labels = cell_state_label,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

bp_res = run.prism(prism = myPrism, n.cores = 80)
bp_res
slotNames(bp_res)

message("*** Saving RData ***")
save(bp_res, file="sft.bp.res.rdata")

theta = get.fraction(bp = bp_res, 
    which.theta="final",
    state.or.type="type") 
write.csv(theta, file="theta_final_type.csv")

get.fraction(bp = bp_res, 
    which.theta="first",
    state.or.type="type") %>%
    write.csv(file="theta_first_type.csv")

theta_state = get.fraction(bp = bp_res, 
    which.theta="first",
    state.or.type="state") %>%
    write.csv(file="theta_first_state.csv")

theta_cv = bp_res@posterior.theta_f@theta.cv
head(theta_cv)
 
# extract posterior mean of cell type-specific gene expression count matrix Z
message("*** Saving z-matrix of cell types ***")
cell_type = colnames(theta_cv)
for (i in cell_type){
    print(i)
    if (i == "T/NK"){
    z_matrix = get.exp(bp=bp_res, state.or.type="type", cell.name = i) %>%
        as.data.frame() %>%
        write_tsv("z_matrix_T-NK.tsv")
    } else {
    z_matrix = get.exp(bp=bp_res, state.or.type="type", cell.name = i) %>%
        as.data.frame() %>%
        write_tsv(paste0("z_matrix_",i,".tsv"))
    }
}

## visualization
theta_state = read_csv("theta_first_type.csv") %>%
  column_to_rownames("...1") %>%
  as.data.frame() %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("cell_type") %>%
  pivot_longer(cols = -cell_type,
               names_to = "sample", values_to = "ratio") %>%
  mutate(sample = str_split(sample, "-", simplify = T)[,1]) %>%
  left_join(sft_cluster_info, by = c("sample" = "PatientID")) %>%
  mutate(sample = as.character(sample),
         Grade = paste0("WHO ", Grade),
         cluster = factor(as.character(cluster), levels = c("4","3","2","1")),
         cell_type = dplyr::recode(cell_type, 
                                   "SFT0" = "SFT_classical",
                                   "SFT1" = "SFT_inflamed",
                                   "SFT2" = "SFT_neural", 
                                   "SFT3" = "SFT_migratory",
                                   "SFT4" = "SFT_angiogenic",
                                   "Oligo"= "Fibroblast", 
                                   "Neuron" = "Unknown")) %>%
  mutate(Cluster = recode(cluster, `1` = "Inflamed",
                          `2` = "Classical",
                          `3` = "Neural-like",
                          `4` = "Migratory"),
         Cluster = factor(Cluster, levels = c("Classical","Neural-like","Inflamed","Migratory"))) %>%
  mutate(cell_label = ifelse(str_detect(cell_type, "SFT"), "Tumor", "Non-tumor")) %>%
  filter(is.na(Status), Location == "CNS")


write_tsv(theta_state, "cell_state_sample_ratio.tsv")

# cell ratio test
cell_grade_sig = map(unique(theta_state$cell_type), MulitCompare) %>% bind_rows() %>%
  dplyr::select(-mean_diff) %>% 
  # filter(p_adj < 0.05, str_detect(cell_type, "SFT")) %>%
  mutate(sig = case_when(
    p_adj<0.05 & p_adj>0.01 ~ "*",
    p_adj<0.01 & p_adj>0.001 ~ "**",
    p_adj<0.001 ~ "***",
    TRUE ~ ""
  ))

cell_grade_sig %>%
  filter(p_adj < 0.05, !str_detect(cell_type, "SFT"))

# ratio visualization
# boxplot
theta_state %>%
  filter(cell_label == "Tumor") %>%
  mutate(cell_label = ifelse(str_detect(cell_type, "SFT"), "Tumor", "Non-tumor")) %>%
  ggplot() +
  geom_hline(yintercept = seq(0.1,0.8,0.2), linetype = "dashed", linewidth = 0.5, color = "grey") +
  geom_boxplot(aes(x = cell_type, y = ratio, color = Cluster), size = 0.4, 
               outlier.colour = NA, 
               outlier.size = 0) +
  geom_jitter(aes(x = cell_type, y = ratio, color = Cluster),
              size = 0.35, alpha = 0.6,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.5)) +
  #facet_grid(cols = vars(cell_label), space = "free", scales = "free") +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "cell ratio") +
  scale_color_manual(values = custom_colors$Cluster) +
  theme_bw() +
  mytheme +
  theme(
    legend.position = c(0.8,0.8),
    legend.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  )

ggsave("cell_cluster_ratio.pdf", 
       units = "mm", height = 90, width = 100)


theta_state %>%
  filter(cell_label != "Tumor", cell_type != "Unknown") %>%
  mutate(cell_label = ifelse(str_detect(cell_type, "SFT"), "Tumor", "Non-tumor")) %>%
  ggplot() +
  geom_hline(yintercept = seq(0.05,0.25,0.05), linetype = "dashed", linewidth = 0.5, color = "grey") +
  geom_boxplot(aes(x = cell_type, y = ratio, color = Cluster), size = 0.4, 
               outlier.colour = NA, 
               outlier.size = 0) +
  geom_jitter(aes(x = cell_type, y = ratio, color = Cluster),
              size = 0.35, alpha = 0.6,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.5)) +
  #facet_grid(cols = vars(cell_label), space = "free", scales = "free") +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "cell ratio", limits = c(0, 0.27)) +
  scale_color_manual(values = custom_colors$Cluster) +
  theme_bw() +
  mytheme +
  theme(
    legend.position = c(0.85,0.8),
    legend.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  )

ggsave("microenv_cluster_ratio.pdf", 
       units = "mm", height = 90, width = 100)

# barplot
ggplot(theta_state) +
  geom_bar(aes(x = as.character(sample), y = ratio, fill = cell_type),
           stat = "identity", width = 0.7,size = 0.5)+
  facet_grid(vars(Cluster), space = "free", scale = "free") +
  theme_classic() +
  scale_fill_manual(values = custom_colors$Cell) +
  scale_x_discrete(labels = NULL) +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  mytheme +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
# ggsave("3_results/3_RNA-seq/4_singleCell/final90/cell_state_sample_ratio.pdf",
#        width = 110, height = 210, units = "mm")

# heatmap
theta_state_mat = theta_state %>%
  #    filter(cell_type != "T cell",cell_type != "Unknown" ) %>%
  dplyr::select(cell_type, sample, ratio) %>%
  pivot_wider(names_from = sample, values_from = ratio) %>%
  column_to_rownames("cell_type") %>%
  as.matrix() 


pdf("cell_ratio_heatmap.pdf", height = 5, width = 7)
Heatmap(theta_state_mat, 
        column_split = theta_state[match(colnames(theta_state_mat), theta_state$sample), ]$Cluster,
        row_split = theta_state[match(rownames(theta_state_mat), theta_state$cell_type), ]$cell_label,
        top_annotation = HeatmapAnnotation(Cluster = theta_state[match(colnames(theta_state_mat), theta_state$sample), ]$Cluster,
                                           #ratio = anno_barplot(t(apply(theta_state_mat, 2, function(x) x/sum(x)))),
                                           #Location = theta_state[match(colnames(theta_state_mat), theta_state$sample), ]$Location,
                                           #Batch = theta_state[match(colnames(theta_state_mat), theta_state$sample), ]$Batch),
                                           Grade = theta_state[match(colnames(theta_state_mat), theta_state$sample), ]$Grade,
                                           col = custom_colors,
                                           annotation_name_side = "left",
                                           annotation_name_gp = gpar(col = "grey20", fontsize = 8, fontface = "bold"), # row legend title
                                           simple_anno_size = unit(0.3, "cm"), # row legend height
                                           gap = unit(c(rep(2,2),3), "points")),
        left_annotation = rowAnnotation(Cell = theta_state[match(rownames(theta_state_mat), theta_state$cell_type), ]$cell_type,
                                        Ratio = row_anno_boxplot(theta_state_mat, 
                                                                 gp = gpar(fill = alpha(rb_two[1], 0.5)),
                                                                 axis_param = list(direction = "reverse")),
                                        col = custom_colors,
                                        annotation_name_gp = gpar(col = "grey20", fontsize = 8, fontface = "bold"), # row legend title
                                        simple_anno_size = unit(0.3, "cm"), # row legend height
                                        gap = unit(c(2,3), "points")),
        row_gap = unit(1.5, "mm"), column_gap = unit(1.5, "mm"), border = TRUE,
        col = colorRampPalette(c("grey90", rb_two[2]))(50),
        column_labels = rep("", ncol(theta_state_mat)),
        row_names_gp =  gpar(fontsize = 8, fontface = "bold"),
        name = "Cell ratio") %>%
  draw(legend_grouping = "original", merge_legends = TRUE)
dev.off()

