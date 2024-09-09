# environment
workdir = "/media/data5/hxy/SFT/RNAseq/2_Brain_SFT/MeningiomaBrainSFT"
setwd(workdir)

# packages
library(tidyverse)
library(Seurat)
library(reshape)
library(harmony)
library(patchwork)

# sample info
sft_data = read_tsv("../../0_RawExpr/sample_rna_cns.tsv") %>%
    column_to_rownames("PatientID") %>%
    mutate(Tissue = "SFT") %>%
    dplyr::select(Tissue, Grade, Cluster) 

gtex_data  = read_tsv("../../0_RawExpr/GTEx_sample_info.tsv", col_names = c("Sample", "Cluster")) %>%
    column_to_rownames("Sample") %>%
    mutate(Grade = Cluster, Tissue = "Normal") %>%
    filter(Grade != "cerebellum", Grade != "ch")

pnas_info = read_csv("../Meningioma/Bulk_method/2019PNAS/pnas_sample_info.csv") %>%
    filter(!is.na(Identifyers)) %>%
    dplyr::select(-Cohort,-Sample) %>%
    column_to_rownames("Identifyers") %>%
    dplyr::rename(Grade = Pathology) %>%
    mutate(Grade = recode(Grade, "WHO I" = "WHO1", "WHO II" = "WHO2", "WHO III" = "WHO3"),
        Tissue = "Meningioma1",
        Cluster = Grade) %>%
    dplyr::select(Tissue, Grade, Cluster)

sa_info = read_csv("../Meningioma/Bulk_method/2022SA/GSE189672_Clinical_data.csv") %>%
    column_to_rownames("Sample") %>%
    dplyr::select(-Gender) %>%
    dplyr::rename(Grade = `WHO grade`) %>%
    mutate(Grade = recode(Grade, "WHO I" = "WHO1", "WHO II" = "WHO2", "WHO III" = "WHO3"),
        Tissue = "Meningioma2",
        Cluster = Grade) %>%
    dplyr::select(Tissue, Grade, Cluster)


colData = bind_rows(sft_data, gtex_data, pnas_info, sa_info) %>%
    rownames_to_column("Sample") %>%
    mutate(batch_id = paste(Tissue, Sample, sep = "_")) %>%
    column_to_rownames("Sample") 

write_tsv(rownames_to_column(colData, "Sample"), "MeningiomaBrainSFT_colData.tsv")

# epxression matrix
protein_map = read_tsv("../../0_RawExpr/protein_id.map", col_names = c("geneID", "geneName")) %>%
    mutate(gene_id = paste(geneID, geneName, sep = "|")) %>%
    mutate(gene = str_split(geneID,"\\.", simplify = T)[,1])

sft_expr = read_csv("../../0_RawExpr/gene_count_matrix.csv") %>%
    dplyr::select(gene_id, any_of(rownames(sft_data))) %>%
    filter(gene_id %in% protein_map$gene_id) %>%
    mutate(mean = rowMeans(across(2:(nrow(sft_data))))) %>%
    filter(mean > 3) %>%
    dplyr::select(-mean) %>%
    separate(gene_id, sep = "\\|", into = c("geneID","geneName")) %>%
    dplyr::select(-geneName) %>%
    separate(geneID, sep = "\\.", into = c("geneID", "version")) %>%
    dplyr::select(-version)

gtex_expr = read_tsv("../BrainNormal/GTEx_total_count_clean.tsv") %>%
    dplyr::select(geneID, any_of(rownames(gtex_data))) %>%
    separate(geneID, sep = "\\.", into = c("geneID", "version")) %>%
    filter(geneID %in% protein_map$gene) %>%
    dplyr::select(geneID, any_of(rownames(gtex_data))) 


pnas_expr = read_tsv("../Meningioma/Bulk_method/2019PNAS/pnas_htseq_merge.tsv") 
colnames(pnas_expr) = c("geneID",(str_split(colnames(pnas_expr), "_", simplify = T)[,2])[-1])

pnas_expr = pnas_expr %>%
    dplyr::select(geneID, any_of(rownames(pnas_info))) %>%
    mutate(mean = rowMeans(across(2:(nrow(pnas_info))))) %>%
    filter(mean > 3) %>%
    dplyr::select(-mean)

sa_expr = read_csv("../Meningioma/Bulk_method/2022SA/GSE189672_RNAseq_raw_counts.csv") %>%
    dplyr::rename(geneID = "...1") %>%
    mutate(geneID = str_split(geneID, "\\.", simplify = T)[,1])
colnames(sa_expr) = str_replace(colnames(sa_expr), "01-0", "") 

sa_expr = sa_expr %>%
    dplyr::select(geneID, any_of(rownames(sa_info))) %>%
    mutate(mean = rowMeans(across(2:(nrow(sa_info))))) %>%
    filter(mean > 3) %>%
    dplyr::select(-mean)

merged_expr = Reduce(full_join, list(sft_expr, gtex_expr, pnas_expr, sa_expr))
merged_expr[is.na(merged_expr)] = 0


# build object
sceList = lapply(c("SFT", "Normal", "Meningioma1", "Meningioma2"), function(batch) {
    samples = colData %>%
        filter(Tissue == batch) %>%
        rownames()

  sce = merged_expr %>%
    select(geneID, any_of(samples)) %>%
    column_to_rownames("geneID") %>%
    as.matrix() %>%
    CreateSeuratObject(counts = ., project = batch) %>%
    NormalizeData() %>%
    FindVariableFeatures(., selection.method = "vst", nfeatures = 2000)
    
    return(sce)
})

# merge batches
message("*** merging data ***")
sce.all = merge(sceList[[1]], y = sceList[-1], add.cell.ids = c("SFT", "Normal", "Meningioma1", "Meningioma2"))
sce.all$grade = colData[match(names(sce.all$orig.ident), colData$batch_id), "Grade"]
sce.all$cluster = colData[match(names(sce.all$orig.ident), colData$batch_id), "Cluster"]

# normailiza again
message("*** normalizaing data ***")
sce.all = NormalizeData(sce.all, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce.all = FindVariableFeatures(sce.all)
sce.all = ScaleData(sce.all)
sce.all = RunPCA(sce.all, features = VariableFeatures(object = sce.all))

# batch re-correction
message("*** batch removal ***")
sce.all = RunHarmony(sce.all, group.by.vars = "orig.ident") 
print(sce.all[["harmony"]], dims = 1:5, nfeatures = 5)

# cluster
message("*** cluster ***")
sce.all = RunUMAP(sce.all, reduction = "harmony", dims = 1:30)
sce.all = FindNeighbors(sce.all, reduction = "harmony", dims = 1:30)
sce.all = FindClusters(sce.all, resolution = 0.5)

# PCA and UMP
message("*** visualing ***")
p1 = DimPlot(sce.all, reduction = "pca", group.by = "orig.ident")
p2 = DimPlot(sce.all, reduction = "umap", group.by = "orig.ident")
p3 = DimPlot(sce.all, reduction = "pca", group.by = "grade")
p4 = DimPlot(sce.all, reduction = "umap", group.by = "grade")
p5 = DimPlot(sce.all, reduction = "pca", group.by = "cluster")
p6 = DimPlot(sce.all, reduction = "umap", group.by = "cluster")
ggsave("MeningiomaBrainSFT_PCA_plot.pdf", plot = (p1+p2)/(p3+p4)/(p5+p6), unit = "mm", height = 500, width = 500)


# final expression count
message("*** saving ***")
sce.all = JoinLayers(object = sce.all)
corrected_matrix = as.data.frame(sce.all@assays$RNA$scale.data) %>%
    rownames_to_column("gene_id")
write_tsv(corrected_matrix, "MeningiomaBrainSFT_corrected_scale.data.tsv")

normalized_matrix = as.data.frame(sce.all@assays$RNA$data) %>%
    rownames_to_column("gene_id")
write_tsv(normalized_matrix, "GTEx_sft_normalized_scale.data.tsv")


# special gene
sm_specific_sets = protein_map %>% filter(geneName %in% c("FER","FRK","FES","STAT6")) %>%
   separate(geneID, sep = "\\.", into = c("geneID", "version"))

sm_merged_epxr_longer = normalized_matrix %>%
    filter(gene_id %in% sm_specific_sets$geneID) %>%
    pivot_longer(cols = -gene_id, names_to = "sample", values_to = "scaled_expr") %>%
    mutate(Tissue = str_split_fixed(sample, "_",2)[,1],
        Sample = str_split_fixed(sample, "_",2)[,2]) %>%
    left_join(colData, by = c("sample" = "batch_id", "Tissue" = "Tissue")) %>%
    left_join(sm_specific_sets, by = c("gene_id" = "geneID"))


write_tsv(sm_merged_epxr_longer, "sm_merged_epxr_longer.tsv")

save(list = ls(), file = "MeningiomaBrainSFT_Seurat_batch.RData")
