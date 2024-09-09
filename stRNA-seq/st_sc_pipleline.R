# environment
options(future.globals.maxSize = 8 * 1024^3)

workdir = '/media/data5/hxy/SFT/RNAseq/4_Spatial'
setwd(workdir)
args = commandArgs(trailingOnly = TRUE)
sample_prefix = args[1]

# data
scRNAseq_ref = '../3_BayesPrism/SFT_raw.rds'
scRNAseq_ref_label = "cell_type2"

used_thread = 60
custom_colors$Cell = c("SFT_classical" = "#F8766D", "SFT_inflamed" = "#7CAE00", "SFT_neural" = "#00BFC4",
                       "SFT_migratory" = "#C77Cff", "SFT_angiogenic" = "#E68613", "Peri" = "#8494FF",
                       "Endo" = "#00B8E7", "Macro" = "#CD9600", "T cell" ="#FF68A1",
                       "Fibroblast" = "#FF61CC", "Unknown" = "#00BE67")

# packages
library(hdf5r)
library(arrow)
library(Seurat)
library(tidyverse)
library(patchwork)
# cluster sort
library(ape)
# faster in DEG of clusteres
library(presto)
# deconvolution
library(spacexr)
library(Matrix)

# data
message("*** Sample is ", sample_prefix, " ***")

message("*** Reading Data ***")
tmp1 = Load10X_Spatial(paste0(sample_prefix, "/binned_outputs/square_008um"))

# total sequencing depth
p1 = SpatialFeaturePlot(tmp1, features = "nFeature_Spatial") + theme(legend.position = "right")
p2 = SpatialFeaturePlot(tmp1, features = "nCount_Spatial") + theme(legend.position = "right")
ggsave(paste0(sample_prefix, "_sequencing_depth.pdf"), plot = p1+p2, unit = "mm", width = 80, height = 80)

# data clean
# mitochondria
message("*** Cleaning Data ***")
tmp1[["percent.mt"]] = PercentageFeatureSet(object = tmp1, pattern = "^MT-")

# gene expressed less than 200 cells
# cells expressed less than 200 or more than 4000 genes
# mitochondria precentage less than 10%
tmp1 = subset(tmp1, subset = rowSums(tmp1@assays$Spatial@layers$counts > 0) >= 200 &
     nFeature_Spatial >= 200 & nFeature_Spatial <= 4000 & percent.mt < 10)

# normalization and identify feature
message("*** find feature ***")
DefaultAssay(tmp1) = "Spatial"
tmp1 = NormalizeData(tmp1)
tmp1 = FindVariableFeatures(tmp1)
tmp1 = ScaleData(tmp1)

# select 20% | 50000 cells as sketch assay
message("*** Sketch data ***")
DefaultAssay(tmp1) = "Spatial"
tmp1 = SketchData(
  object = tmp1,
  ncells = 0.2*(ncol(tmp1)),
  method = "LeverageScore",
  sketched.assay = "sketch")

# find features and clustering in sketched data (too large in full dataset)
message("*** Sketch uMAP ***")
DefaultAssay(tmp1) = "sketch"

tmp1 = FindVariableFeatures(tmp1)
tmp1 = ScaleData(tmp1)
tmp1 = RunPCA(tmp1, assay = "sketch", reduction.name = "pca.sketch")
tmp1 = FindNeighbors(tmp1, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
tmp1 = FindClusters(tmp1, cluster.name = "seurat_cluster.sketched", resolution = 3)
tmp1 = RunUMAP(tmp1, reduction = "pca.sketch", reduction.name = "umap.sketch",
                  return.model = T, dims = 1:50)

# project to full sets
message("*** Project data ***")
tmp1 = ProjectData(
  object = tmp1,
  assay = "Spatial",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched") # ref data
)

# umap compare of sketched and full data
message("*** sketch and project compare ***")
DefaultAssay(tmp1) = "sketch"
Idents(tmp1) = "seurat_cluster.sketched"
p1 = DimPlot(tmp1, reduction = "umap.sketch", label = F) + 
    ggtitle("Sketched clustering") +
    theme(legend.position = "bottom")

DefaultAssay(tmp1) = "Spatial"
Idents(tmp1) = "seurat_cluster.projected"
p2 = DimPlot(tmp1, reduction = "full.umap.sketch", label = F) + 
    ggtitle("Projected clustering (full dataset)") + 
    theme(legend.position = "bottom")
ggsave(paste0(sample_prefix, "_umap_compare.pdf"), plot = p1 + p2, unit = "mm", width = 300, height = 300)

# cluster distribution visualization
p = SpatialDimPlot(tmp1, label = T, repel = T, label.size = 4)
ggsave(paste0(sample_prefix, "_cluster_distribution.pdf"), plot = p, unit = "mm", width = 300, height = 120)

# save data
message("*** Saving Data ***")
saveRDS(tmp1, file = paste0(sample_prefix, ".RDS"))

# identification cluster top gene
message("*** identify cluster trees ***")
DefaultAssay(tmp1) = "Spatial"
Idents(tmp1) = "seurat_cluster.projected"
object_subset = subset(tmp1, cells = Cells(tmp1[["Spatial"]]), downsample = 0.1*(ncol(tmp1)))

# sort by similarity of cluster (ape packages)
DefaultAssay(object_subset) = "Spatial"
Idents(object_subset) = "seurat_cluster.projected"
object_subset = BuildClusterTree(object_subset, assay = "Spatial", reduction = "full.pca.sketch", reorder = T)

# identify differenctial expressed genes of clusters (faster when using prest)
message("*** identify cluster markers ***")
markers = FindAllMarkers(object_subset, assay = "Spatial", only.pos = TRUE)

message("*** output cluster markers ***")
markers %>%
    as.data.frame() %>%
    filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
    write_tsv(paste0(sample_prefix, "_sketch_subset_markers.tsv"))

# top gene expression heatmap
top5 = markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() 

object_subset = ScaleData(object_subset, assay = "Spatial", features = top5$gene)
p = DoHeatmap(object_subset, assay = "Spatial", features = top5$gene, size = 2.5) + 
    theme(axis.text = element_text(size = 5.5)) + NoLegend()
ggsave(paste0(sample_prefix, "_cluster_top5gene.pdf"), plot = p, unit = "mm", width = 800, height = 400)

# save seurat object with cluster information
message("*** Saving Data ***")
saveRDS(tmp1, file = paste0(sample_prefix, ".RDS"))

# merge with scRNA-seq
message("*** Load and precess scRNA-seq ***")
ref = readRDS(scRNAseq_ref)

Idents(ref) = scRNAseq_ref_label
counts = ref[["originalexp"]]$counts
cluster = as.factor(ref$cell_type2)
nUMI = ref$nCount_originalexp
levels(cluster) = gsub("/", "-", levels(cluster)) # NK/T
cluster = droplevels(cluster)

# RCTD preparation
message("***Start RCTD (spacexr package) ***")
reference = Reference(counts, cluster, nUMI)

# full object thus not need to project
counts_hd = tmp1[["Spatial"]]$counts
tmp1_cells_hd = colnames(tmp1[["Spatial"]])
coords = GetTissueCoordinates(tmp1)[tmp1_cells_hd, 1:2]

# construct RCTD query object
message("*** RCTD query construction ***")
query = SpatialRNA(coords, counts_hd, colSums(counts_hd))

# start RCTD analysis
message("*** RCTD analysis ***")
RCTD = create.RCTD(query, reference, max_cores = used_thread)
RCTD = run.RCTD(RCTD, doublet_mode = "doublet") 

# add RCTD to ST object
message("*** add RCTD results ***")
tmp1 = AddMetaData(tmp1, metadata = RCTD@results$results_df)

# adding deconvolution information of scRNA-seq
message("*** cell types calculated by RCTD ***")
tmp1$first_type = as.character(tmp1$first_type)
print(table(tmp1$first_type))

tmp1$first_type[is.na(tmp1$first_type)] = "Unknown"

# save seurat object with deconvolution information
message("*** Saving Data ***")
saveRDS(tmp1, file = paste0(sample_prefix, ".RDS"))

# visual pca and umap
message("*** umap visualization ***")

DefaultAssay(tmp1) = "Spatial"
Idents(tmp1) = "first_type"

p = DimPlot(tmp1, reduction = "full.umap.sketch", label = F, cols = used_col) + 
    ggtitle("Projected clustering (full dataset)") + 
    theme(legend.position = "bottom")
ggsave(paste0(sample_prefix, "_scumap_compare.pdf"), plot = p, unit = "mm", width = 150, height = 150)

# cell position visualization
message("*** Visual cell position ***")
DefaultAssay(tmp1) = "Spatial"
Idents(tmp1) = "first_type"

cells = CellsByIdentities(tmp1)

#p = SpatialDimPlot(tmp1, cells.highlight = cells[selected_names], 
#    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T, ncol = 4)
p = SpatialDimPlot(tmp1, cells.highlight = cells, 
    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T, ncol = 4)
ggsave(paste0(sample_prefix, "_tumor_cell_position.pdf"), plot = p, unit = "mm", width = 600, height = 600)

message("*** cell marker visualization ***")
DefaultAssay(tmp1) ="Spatial"
Idents(tmp1) = "first_type"

CellTypeMarkerPlot = function(obj = tmp1, cell = "SFT0"){
    markers_CellType = FindMarkers(obj, ident.1 = cell, only.pos = TRUE)

    markers_CellType %>%
        as.data.frame() %>%
        filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
        rownames_to_column("gene_name") %>%
        write_tsv(paste0(sample_prefix, "_",cell, "_markergene.tsv"))

    top5gene = markers_CellType %>%
        filter(avg_log2FC > 1) %>%
        head(5)
     if (nrow(top5gene) > 0) {

        p1 = FeaturePlot(obj, features = rownames(top5gene))
        p2 = SpatialFeaturePlot(obj, features = rownames(top5gene), combine = TRUE) + 
          ggtitle(paste0("Top 5 Markers in ", cell, " Cells")) +
          theme_minimal()
        p = p1 / p2
        ggsave(paste0(sample_prefix, "_",cell, "_markergene.pdf"), plot = p, unit = "mm", width = 300, height = 600)
        message("*** plot ", cell, " ***")

    } else {
        message("*** No top 5 genes for ", cell, " ***")
    }
   
}

walk(names(cells), ~CellTypeMarkerPlot(obj = tmp1, cell = .x))

p = SpatialFeaturePlot(tmp1, features = c("CD34","CD44","ALDH1A1","ENG"), combine = TRUE)
ggsave(paste0(sample_prefix, "_given_markergene.pdf"), plot = p, unit = "mm", width = 200, height = 200)

# plot
# functions
source("../functions.R")
# data
message("*** Reading Data ***")
message("*** Sample is ", sample_prefix, " ***")
st_obj = readRDS(paste0("RDS/", sample_prefix, ".RDS"))

message("*** Renaming ***")
st_obj = RenameCellType(st_obj)

#message("*** Ploting ***")
PlotCellPosition(st_obj, sample_prefix)
PlotSpecialCellPosition(st_obj, "^SFT", cell_colors, sample_prefix)
PlotCellPositionSeperate(st_obj, sample_prefix)
PlotSpecialCellPosition(st_obj, "nflamed|igratory|ngiogenic|Endo", 
                        c("red","blue","gold","green","lightgrey"), paste(sample_prefix,"_sft"))

#message("*** Markers ***")
walk(names(cell_markers), ~ PlotCellTypeMarker(obj = st_obj, cell = .x, sample_prefix = sample_prefix))
PlotCellTypeMarker(obj = st_obj, cell = "Macro", sample_prefix = sample_prefix)

message("*** Saving ***")
saveRDS(st_obj, paste0(sample_prefix,"_new.RDS"))
