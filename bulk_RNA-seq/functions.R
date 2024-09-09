################## Genome ##################
# filter snv by gnomeAD and location
SNVfilter = function(esnv){
    esnv %>% 
        mutate(gnomad40_exome_AF = ifelse(gnomad40_exome_AF == ".", 0, as.numeric(gnomad40_exome_AF)),
               gnomad40_genome_AF = ifelse(gnomad40_genome_AF == ".", 0, as.numeric(gnomad40_genome_AF)),
               gnomad40_exome_AF_eas = ifelse(gnomad40_exome_AF_eas == ".", 0, as.numeric(gnomad40_exome_AF_eas)),
               gnomad40_genome_AF_eas = ifelse(gnomad40_genome_AF_eas == ".", 0, as.numeric(gnomad40_genome_AF_eas))) %>% 
        filter(!(gnomad40_exome_AF > 0.001 | gnomad40_genome_AF > 0.001 | gnomad40_exome_AF_eas > 0.001 | gnomad40_genome_AF_eas > 0.001)) %>%
        .[,c(7,9,85:175)] %>%
        pivot_longer(cols = contains("Otherinfo"), names_to = "sample", values_to = "GT") %>%
        mutate(GT = ifelse(str_starts(GT, "\\./\\."), 0, 1),
               sample = paste0("SFT",str_pad(as.numeric(str_sub(sample, 10,12))-12, width = 3, side = "left", pad = 0))) %>%
        filter(GT == 1, ExonicFunc.refGene %in% c("frameshift deletion",
                                                  "frameshift insertion", "nonsynonymous SNV",
                                                  "startloss", "stopgain", "stoploss")) %>%
        distinct() %>%
        inner_join(filter(sample_info_merge, Seq_type == "WGS"), by = c("sample" = "SampleID"))
}



# self dotplot
SelfDotplot = function(enrich){
    ggplot(enrich, aes(y = str_wrap(Description, 30), x = enrichfold,
                       color = p.adjust, size = enrichfold, shape = Dtype)) +
        geom_point() +
        facet_grid(rows = vars(Etype), scales = "free", space = "free") +
        scale_color_gradient2(high = "#76ABD5", mid = "#8A4F80") +
        scale_size_continuous(range = c(1,4)) +
        labs(size = "Enrichment Fold", shape = "mutations", 
             x = "Enrichment Fold", y = "Description") +
        theme_bw() +
        mytheme
}

################## RNA-seq ##################
NFStest = function(term, data) {
    tmp = filter(data, terms == term)
    
    data.frame(
        "Term" = term,
        g1p = wilcox.test(pull(filter(tmp, Grade == 1), NES), pull(filter(tmp, Grade != 1), NES), alternative = "greater")$p.value,
        g2p = wilcox.test(pull(filter(tmp, Grade == 2), NES), pull(filter(tmp, Grade != 2), NES), alternative = "greater")$p.value,
        g3p = wilcox.test(pull(filter(tmp, Grade == 3), NES), pull(filter(tmp, Grade != 3), NES), alternative = "greater")$p.value,
        c1p = wilcox.test(pull(filter(tmp, cluster == 1), NES), pull(filter(tmp, cluster != 1), NES), alternative = "greater")$p.value,
        c2p = wilcox.test(pull(filter(tmp, cluster == 2), NES), pull(filter(tmp, cluster != 2), NES), alternative = "greater")$p.value,
        c3p = wilcox.test(pull(filter(tmp, cluster == 3), NES), pull(filter(tmp, cluster != 3), NES), alternative = "greater")$p.value,
        c4p = wilcox.test(pull(filter(tmp, cluster == 4), NES), pull(filter(tmp, cluster != 4), NES), alternative = "greater")$p.value
    )
}


SigNESterms = function(test_data, ptype, n_term) {
    filtered_data = test_data %>%
        filter(p_type == ptype) %>%
        filter(p_adj < 0.05) %>%
        arrange(p_adj) %>%
        slice_head(n = n_term) %>%
        dplyr::select(Term, p_type)
    
    left_join(filtered_data, test_data, by = c("Term", "p_type"))
}

SigNEStermsHeatmap = function(n_terms){
    sft_seurat_nmf_nes_top = map_dfr(unique(sft_seurat_nmf_nes_tests$p_type), 
                                     ~ SigNESterms(sft_seurat_nmf_nes_tests, .x, n_terms))
    
    sft_nes_top = sft_seurat_nmf_nes[filter(sft_seurat_nmf_nes_top, p_type %in% c("c1p", "c2p", "c3p", "c4p"))$Term, ] %>%
        distinct()
    #tmp = sft_seurat_nmf_nes[sft_seurat_nmf_nes_top$Term, ]
    colnames(sft_nes_top) = str_split(colnames(sft_nes_top), pattern = "_", simplify = T)[,2]
    rownames(sft_nes_top)
    
    pdf("h_go_2k.pdf", width = 15, height = 6)
    sft_nes_top[,rownames(arrange(sft_cluster_mad_anno, desc(cluster)))] %>%
        pheatmap::pheatmap(scale = "row", annotation_col = sft_cluster_mad_anno, cluster_cols = F,
                           cutree_rows = 4, cutree_cols = 4, show_colnames = F, show_rownames = T, border = NA)
    dev.off()
    
    # pdf("c7_2k_cluster.pdf", width = 10, height = 6)
    # sft_nes_top %>% pheatmap(scale = "row", annotation_col = sft_cluster_mad_anno, cluster_cols = T,
    #                          cutree_rows = 4, show_colnames = F, show_rownames = F, border = NA)
    # dev.off()
    
}

cluster_volcano = function(res, comp, flt = T){
    
    if (flt){
        filtered_res = res %>% 
            filter(comparison == comp) %>%
            mutate(sig = ifelse((adj.P.Val < 0.05 & abs(logFC) > log2(1.5)), 
                                ifelse(logFC > 0, "up", "down"), "non-sig")) %>%
            rownames_to_column("ID") %>%
            separate(ID, sep = "-", into = c("geneID", "geneName"))
    } else {
        filtered_res = res %>% 
            mutate(sig = ifelse((adj.P.Val < 0.05 & abs(logFC) > log2(1.5)), 
                                ifelse(logFC > 0, "up", "down"), "non-sig")) %>%
            rownames_to_column("ID") %>%
            separate(ID, sep = "-", into = c("geneID", "geneName"))
    }
    
    
    top_up = filtered_res %>%
        filter(sig == "up") %>%
        arrange(adj.P.Val) %>%
        head(5)
    
    top_down = filtered_res %>%
        filter(sig == "down") %>%
        arrange(adj.P.Val) %>%
        head(5)
    
    ggplot(filtered_res, aes(x = logFC, y = -log10(adj.P.Val), 
                             size = -log10(adj.P.Val), color = sig)) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
        geom_vline(xintercept = c(-(log2(1.5)), log2(1.5)), linetype = "dashed", color = "#999999") +
        geom_point(alpha = 0.5) +
        scale_color_manual(values = c(rb_two[1], "grey", rb_two[2])) +
        scale_size_continuous(range = c(0.5, 2.5)) +
        geom_text_repel(data = top_up, aes(label = geneName), size = 3, box.padding = 0.3, point.padding = 0.3) +
        geom_text_repel(data = top_down, aes(label = geneName), size = 3, box.padding = 0.3, point.padding = 0.3) +
        labs(title = comp) +
        theme_bw() +
        mytheme +
        theme(panel.grid = element_blank())
}

cluster_enrich_out = function(res = cluster_res, comp){
    filtered_res = res %>% 
        rownames_to_column("ID") %>%
        separate(ID, sep = "-", into = c("gene_id", "geneName")) %>%
        separate(gene_id, sep = "\\.", into = c("geneID", "gene_version")) %>%
        dplyr::rename(log2FoldChange= logFC,
                      Compare = comparison) %>%
        filter(Compare == comp) %>%
        inner_join(bitr(.$geneID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db"),
                   by = c("geneID" = "ENSEMBL")) %>%
        distinct(ENTREZID, log2FoldChange) %>%
        .[!(duplicated(.$ENTREZID)), ] %>%
        arrange(desc(log2FoldChange))
    
    # Prepare named vector for GSEA
    gene_list = filtered_res$log2FoldChange
    names(gene_list) = filtered_res$ENTREZID
    
    print(head(gene_list))
    
    # Perform GSEA
    gsea_result = gsePathway(geneList = gene_list, organism = "human") %>%
        setReadable(OrgDb = "org.Hs.eg.db") %>%
        as.data.frame() %>%
        mutate(Comparison = comp)
    
}

AllMsigdbEnrich = function(comp, cat){
    message("comparison is ", comp)
    message("category is ", cat)
    cluster_res %>% filter(comparison == comp) %>%
        filter(abs(logFC) > log2(1.5), P.Value < 0.05)  %>%
        rownames_to_column("ID") %>%
        separate(ID, sep = "-", into = c("gene_id", "geneName")) %>%
        separate(gene_id, sep = "\\.", into = c("geneID", "gene_version")) %>%
        pull(geneID) %>%
        bitr(fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db") %>%
        pull(SYMBOL) %>%
        OneStepEnrich(enrich_type = "MSigDB", category = cat, res_type = "clean") %>%
        mutate(category = cat,
               comparison = comp)
}

# load OneStepEnrich
devtools::load_all("D://_Scripts/R/Packages/OneStepEnrich/")

# GSVA
OneStepGSEA = function(res_sig = combat_res_sig, df = F, item_num = 1:3, filter = T, compare_group = "G2vs1") {
    # Filter and prepare gene list
    if(!filter){
        tmp = combat_res_sig %>%
            inner_join(bitr(.$geneID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db"),
                       by = c("geneID" = "ENSEMBL")) %>%
            distinct(ENTREZID, log2FoldChange) %>%
            .[!(duplicated(.$ENTREZID)), ] %>%
            arrange(desc(log2FoldChange))
        
    } else {
        tmp = combat_res_sig %>%
            filter(Compare == compare_group) %>%
            inner_join(bitr(.$geneID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db"),
                       by = c("geneID" = "ENSEMBL")) %>%
            distinct(ENTREZID, log2FoldChange) %>%
            .[!(duplicated(.$ENTREZID)), ] %>%
            arrange(desc(log2FoldChange))
        
    }
    
    # Prepare named vector for GSEA
    message("*** STATR ***")
    gene_list = tmp$log2FoldChange
    names(gene_list) <- tmp$ENTREZID
    
    # Perform GSEA
    gsea_result = gsePathway(geneList = gene_list, organism = "human") %>%
        setReadable(OrgDb = "org.Hs.eg.db")
    
    if(!df){
        # Plot GSEA results
        gseaplot2(
            gsea_result, geneSetID = item_num,
            title = compare_group, 
            color = "firebrick",
            base_size = 11,
            rel_heights = c(1.5, 0.5, 1),
            subplots = 1:3,
            pvalue_table = TRUE,
            ES_geom = "line"
        )
        
    } else {
        return(gsea_result)
    }
}

# RNA-seq PCA
OneStepPCA = function(data, row_id, selected_columns, expr_cut, prefix){
    matrix = data %>%
        column_to_rownames(row_id) %>%
        dplyr::select(any_of(selected_columns))
    
    sample_rna_anno = sample_rna_filter %>%
        filter(PatientID %in% colnames(matrix))
    
    matrix = matrix[rowSums(matrix) > expr_cut,]
    
    message("row:", nrow(matrix), " col:", ncol(matrix))    
    
    pca = matrix %>% t() %>%
        PCA(., scale.unit = T, ncp = 5, graph = F)
    
    fviz_eig(pca, addlabels = TRUE, 
             barfill = "#ADA7B0", barcolor = "#ADA7B0", 
             main = "Percentage of variance")
    ggsave(paste0(prefix, "_pca_eig.pdf"), 
           units = "mm", width = 110, height = 110)
    
    fviz_pca_ind(pca, pointsize = "cos2", 
                 pointshape = 21, fill = "#E7B800",
                 repel = TRUE)
    
    fviz_pca_ind(pca,
                 geom.ind = "point", # show points only (nbut not "text")
                 col.ind = as.character(sample_rna_anno$Grade), # color by groups
                 palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                 addEllipses = TRUE, # Concentration ellipses
                 legend.title = "Grade")
    ggsave(paste0(prefix, "_pca_grade.pdf"), 
           units = "mm", width = 110, height = 110)
    
    fviz_pca_ind(pca,
                 geom.ind = "point", # show points only (nbut not "text")
                 col.ind = as.character(sample_rna_anno$Batch), # color by groups
                 palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                 addEllipses = TRUE, # Concentration ellipses
                 legend.title = "Grade")
    ggsave(paste0(prefix, "_pca_batch.pdf"), 
           units = "mm", width = 110, height = 110)
    
    return(pca)
}


# volcano
VolcanoPlot = function(res, p_cut, fc_cut){
    mutate(res, sig = ifelse(padj < p_cut & log2FoldChange > fc_cut, "up", 
                             ifelse(padj < p_cut & log2FoldChange < -(fc_cut), "down", "non-sig"))) %>%
        ggplot(aes(x = log2FoldChange, y = -log10(padj), size = log2FoldChange, color = sig))+
        geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "#999999")+
        geom_vline(xintercept = c(-(fc_cut),fc_cut), linetype = "dashed", color = "#999999")+
        geom_point() +
        scale_color_manual(values = c(rb_two[1], "grey", rb_two[2])) +
        scale_size_continuous(range = c(0.5,1.5)) +
        theme_bw() +
        mytheme +
        theme(panel.grid = element_blank())
}


# Function for differential gene expression analysis of limma-voom
PerformLimmaDEG = function(countdata, coldata, levels = c("CNS", "Periphery")) {
    
    group = factor(coldata[colnames(countdata),]$Location, levels = levels)
    
    # Create DGEList object and calculate normalization factors
    dge = DGEList(counts = countdata, group = group) %>%
        calcNormFactors()
    
    # Create design matrix
    design = model.matrix(~ group)
    
    # Apply voom transformation
    peri_idx = which(group == "Periphery")
    cns_idx = setdiff(seq_along(group), peri_idx)
    
    v = voom(dge[, c(peri_idx, cns_idx)], design)
    
    # Fit the linear model
    fit = lmFit(v, design)
    
    # Apply empirical Bayes moderation
    fit = eBayes(fit)
    
    # Extract all genes results from the linear model fit
    res = topTable(fit, coef = 2, number = Inf, sort.by = "none")
    
    return(res)
}

# network
SubNetViual = function(subnetwork = cns_subnetwork, top_number = 5, file_name = "cns.html"){
    edge_list = subnetwork %>%
        as.data.frame() %>%
        rownames_to_column("from") %>%
        pivot_longer(cols = starts_with("ENSG"), names_to = "to", values_to = "weight") %>%
        filter(weight != 1) %>%
        separate(from, sep = "-", into = c("from_id", "from")) %>%
        separate(to, sep = "-", into = c("to_id", "to")) %>%
        filter(!str_detect(to, "^ENSG")) %>%
        filter(!str_detect(from, "^ENSG")) %>%
        dplyr::select(-from_id, -to_id) %>%
        filter((from %in% location_enrich_gene)) %>%
        group_by(from) %>%
        arrange(desc(weight)) %>%
        dplyr::slice(1:top_number)
    
    node_names = data.frame(
        name = unique(c(edge_list$from, edge_list$to)),
        group = ifelse(unique(c(edge_list$from, edge_list$to)) %in% location_enrich_gene, "1", "0"),
        size = ifelse(unique(c(edge_list$from, edge_list$to)) %in% location_enrich_gene, 1, 0.8)
    ) 
    
    edge_list = edge_list %>%
        mutate(
            source = match(from, node_names$name) - 1,
            target = match(to, node_names$name) - 1
        )
    
    # 使用 forceNetwork 绘制图形
    p = forceNetwork(
        Links = edge_list,
        Nodes = node_names,
        Source = "source",
        Target = "target",
        Value = "weight",
        NodeID = "name",
        Group = "group",  # Use 'group' column for coloring nodes
        Nodesize = "size",
        radiusCalculation = "Math.sqrt(d.nodesize)+8", legend = TRUE,
        colourScale = JS("d3.scaleOrdinal().domain([0, 1]).range(['#003366','#990033'])"),
        charge = -40,
        fontSize = 12,
        fontFamily = "Arial",
        width = 800,
        height = 800,
        opacity = 1,
        opacityNoHover = 1,
        linkDistance = 75,
        linkWidth = JS("function(d) { return d.value * 15; }"),  # Set the link width based on the weight
        #linkColour = ifelse(edge_list$to %in% location_enrich_gene, "#bf3eff", "#666"),
        bounded = T,
        zoom = T
    )
    
    saveWidget(p, file_name)
}

CalculateScore = function(type, target, expr){

    if (type == "file") {
        target_genes = read_tsv(target, col_names = "gene_name") %>% pull(gene_name)    
    } else if (type == "gene") {
        target_genes = target
    } else if (type == "gs") {
        target_genes = msigdbr(species = "Homo sapiens", category = "H") %>%
            filter(gs_name == target) %>%
            pull(gene_symbol)
    } else {
        stop(" file (one column) or gene (vector) or gs (HALLMARK_HYPOXIA)")
    }
    
    # calculate scores
    expr_matrix  = expr %>%
        as.data.frame() %>%
        mutate(gene_name = str_split(rownames(expr), "-", simplify = T)[,2]) %>%
        .[!duplicated(.$gene_name),] %>%
        `rownames<-`(NULL) %>%
        column_to_rownames("gene_name") %>% 
        as.matrix()
    
    scores = gsva(expr_matrix, list(Scores = target_genes), method = "ssgsea")
    scores_df = scores %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("batch_id") %>%
        separate(batch_id, "_", into = c("Batch", "PatientID")) %>%
        right_join(sample_info) %>%
        mutate(Cluster = factor(Cluster, levels = c( 'Classical','Neural-like', 'Inflamed', 'Migratory')))
    
    print(PMCMRplus::kwAllPairsDunnTest(Scores ~ Cluster, data = scores_df, 
                                        p.adjust.method = "bonferroni") %>%
              summary()
    )
    return(scores_df)
}

# st
RenameCellType = function(obj){
    
    message("*** Rename cell types ***")
    
    # rename cell type
    obj$first_type = as.character(obj$first_type)
    obj$first_type[obj$first_type == "Oligo"] = "Fibroblast"
    obj$first_type[obj$first_type == "Neuron"] = "Unknown"
    
    obj$first_type[obj$first_type == "SFT0"] = "SFT_classical"
    obj$first_type[obj$first_type == "SFT1"] = "SFT_inflamed"
    obj$first_type[obj$first_type == "SFT2"] = "SFT_neural"
    obj$first_type[obj$first_type == "SFT3"] = "SFT_migratory"
    obj$first_type[obj$first_type == "SFT4"] = "SFT_angiogenic"
    
    DefaultAssay(obj) = "Spatial"
    Idents(obj) = "first_type"
    
    message("*** Recalculate PCA and umap ***")
    #    obj = FindVariableFeatures(obj)
    #    obj = ScaleData(obj)
    obj = RunPCA(obj, reduction.name = "pca.sc")
    #    obj = FindNeighbors(obj, reduction = "pca.sc", dims = 1:50)
    #    obj = FindClusters(obj, cluster.name = "first_type", resolution = 3)
    obj = RunUMAP(obj, reduction = "pca.sc", reduction.name = "umap.sc",
                  return.model = T, dims = 1:50)
    
    return(obj)
}


PlotCellPosition = function(obj, sample_prefix){
    
    message("*** Visual cell position ***")
    
    DefaultAssay(obj) = "Spatial"
    Idents(obj) = "first_type"
    cells = CellsByIdentities(obj)
    
    p1 = SpatialDimPlot(obj, 
                        image.alpha = 0,
                        cells.highlight = cells, 
                        cols.highlight = c(cell_colors, "grey10"),
                        label = F,
                        #    label.color = "black",
                        #    label.size = 1.8,
                        #    label.box = T,
                        #    repel = T,
                        facet.highlight = F, 
                        combine = T
    )
    
    p2 = SpatialDimPlot(obj, 
                        image.alpha = 0,
                        cells.highlight = cells, 
                        cols.highlight = c(cell_colors, "grey10"),
                        label = F,
                        facet.highlight = T, 
                        combine = T,
                        ncol = 4
    )
    
    #    sample_prefix = deparse(substitute(obj))
    #    print(sample_prefix)
    
    p3 = DimPlot(obj, reduction = "umap.sc", label = F, cols = cell_colors) +
        theme(legend.position = "bottom")
    
    ggsave(paste0(sample_prefix, "_cell_position.pdf"), plot = p1, unit = "mm", width = 120, height = 120)
    ggsave(paste0(sample_prefix, "_cell_position_seperate.pdf"), 
           plot = p2, unit = "mm", width = 400, height = 400)
    ggsave(paste0(sample_prefix, "_scumap.pdf"), plot = p3, unit = "mm", width = 100, height = 120)
    
    ggsave(paste0(sample_prefix, "_cell_position.png"), plot = p1, unit = "mm", width = 120, height = 120, dpi = 600)
    ggsave(paste0(sample_prefix, "_cell_position_seperate.png"), 
           plot = p2, unit = "mm", width = 400, height = 400, dpi = 600)
    ggsave(paste0(sample_prefix, "_scumap.png"), plot = p3, unit = "mm", width = 100, height = 120, dpi = 600)
    
    
}

PlotCellPositionSeperate = function(obj, sample_prefix){
    message("*** Visual cell position ***")
    
    DefaultAssay(obj) = "Spatial"
    Idents(obj) = "first_type"
    cells = CellsByIdentities(obj)
    cell_type = names(cells)
    
    plots = SpatialDimPlot(obj, 
                           image.alpha = 0,
                           cells.highlight = cells, 
                           cols.highlight = c("red", "grey10"),
                           label = F,
                           facet.highlight = T, 
                           combine = F
    )
    
    map2(plots, cell_type, ~ ggsave(filename = paste0(sample_prefix,"_", .y,"_cell_position.pdf"), 
                                    plot = .x, unit = "mm", width = 100, height = 100))
    
    map2(plots, cell_type, ~ ggsave(filename = paste0(sample_prefix,"_", .y,"_cell_position.png"), 
                                    plot = .x, unit = "mm", width = 100, height = 100, dpi = 600))
    
}

PlotSpecialCellPosition = function(obj, cell_pattern, cell_col = cell_colors, sample_prefix){
    
    message("*** Visual Special Cell position ***")
    
    DefaultAssay(obj) = "Spatial"
    Idents(obj) = "first_type"
    cells = CellsByIdentities(obj)
    
    selected_names = sort(grep(cell_pattern, names(cells), value = TRUE))
    
    p1 = SpatialDimPlot(obj, 
                        image.alpha = 0,
                        cells.highlight = cells[selected_names], 
                        cols.highlight = cell_col,
                        label = F,
                        pt.size.factor = 2,
                        alpha = rep(0.75, length(cell_col)),
                        facet.highlight = F, 
                        combine = T
    )
    
    ggsave(paste0(sample_prefix, "_speical_cell_position.pdf"), plot = p1, unit = "mm", width = 110, height = 110)
    ggsave(paste0(sample_prefix, "_speical_cell_position.png"), plot = p1, unit = "mm", width = 110, height = 110, dpi = 600)
    
}

PlotCellTypeMarker = function(obj, cell = "Endo", sample_prefix){
    
    DefaultAssay(obj) ="Spatial"
    Idents(obj) = "first_type"
    
    message("*** plot ", cell, " ***")
    marker_gene = cell_markers[[cell]]
    message("*** marker genes: ", marker_gene, " ***")
    
    p1 = FeaturePlot(obj, features = marker_gene,
                     reduction = "umap.sc",
                     keep.scale = "feature",
                     cols = c("lightgrey","#990033"), 
                     ncol = 3,
                     combine = T)
    
    ggsave(paste0(sample_prefix, "_", cell, "_featureplot.pdf"), 
           plot = p1, unit = "mm", width = 300, height = 100)
    ggsave(paste0(sample_prefix, "_", cell, "_featureplot.png"), 
           plot = p1, unit = "mm", width = 300, height = 100, dpi = 600)
    
    
    plots = SpatialFeaturePlot(obj, image.alpha = 0,
                               features = marker_gene, ncol = 3,
                               combine = F)
    
    plots = lapply(plots, function(p) {
        p + scale_fill_gradientn(colors = c("lightgrey", "#990033"))
    })
    
    
    combined_plot = wrap_plots(plots, ncol = 3)
    
    ggsave(paste0(sample_prefix, "_", cell, "_markergene.pdf"), 
           plot = combined_plot, unit = "mm", width = 300, height = 100)
    ggsave(paste0(sample_prefix, "_", cell, "_markergene.png"), 
           plot = combined_plot, unit = "mm", width = 300, height = 100, dpi = 600)
    
}

