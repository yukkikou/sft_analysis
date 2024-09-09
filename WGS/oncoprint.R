# environment
workdir = "D://_SFT/"
figdir = "D://_SFT/4_figures/"
setwd(workdir)

# packages
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# color and theme
source("../2_scripts/7_plot/theme_colors.R")

# function
source("../2_scripts/functions.R")

oncoDataPrecess = function(plot_df, plot_sample_info = plot_sample_wgs_info, plot_gene){
    
    # get all mutation type
    mutation_tp = unique(plot_df$mutation_type)
    
    result = map(mutation_tp, function(mutation_tp) {
        # build mutation matrix for each mutation type character by 1 and 0
        message("mutation type is ", mutation_tp)
        plot_df_part = plot_df %>%
            filter(mutation_type == mutation_tp) %>%
            dplyr::select(-mutation_type)
        
        plot_df_part = plot_df_part %>% 
            pivot_wider(names_from = sample, values_from = state) %>%
            #arrange(SYMBOL) %>%
            column_to_rownames("SYMBOL") %>%
            as.matrix()
        
        plot_df_part = plot_df_part[plot_gene, ]
        
        print(rownames(plot_df_part))
        return(plot_df_part)

    })
    
    names(result) = mutation_tp
    return(result)
}


# sample info
plot_sample_wgs_info = sample_info_keep %>%
    mutate(Grade = paste("WHO", Grade))

#write_excel_csv(plot_sample_wgs_info, "3_results/Summary/sample_infomation_wgs.csv")
#write_tsv(esnv_all_filter, "1_somaticMerge/esnv_all_filter.tsv")
#write_tsv(fCNV_res_filter, "4_CNV/fCNV_res_filter.tsv")

# mutation enrichment
unique(c(esnv_all_filter$SYMBOL, filter(fCNV_res_filter, CN>4|CN<0.5)$SYMBOL)) %>%
    OneStepEnrich(res_type = "clean") %>%
    mutate(count = as.numeric(str_split(GeneRatio, pattern = "\\/", simplify = T)[,1]),
           enrichfold = DOSE::parse_ratio(GeneRatio)/DOSE::parse_ratio(BgRatio),
           Dtype = "mutation", Etype = "mutation") %>%
    .[c(1:2,4,7,9,12,13),] %>%
    ggplot(aes(y = str_wrap(Description, 30), x = enrichfold,
                       color = p.adjust, size = count)) +
    geom_point() +
    scale_color_gradient2(high = "blue", mid = "red") +
    scale_size_continuous(range = c(1,4)) +
    labs(size = "Count", x = "Enrichment Fold", y = "Description") +
    theme_bw() +
    mytheme
ggsave("3_results/Summary/WGS/wgs_mutation_enrichement.pdf", width = 100, height = 80, units = "mm")

# mutation data
plot_snv_sum = esnv_all_filter %>%
    filter(SYMBOL %in% c(filter(merged_drivers, evi_tag > 2)$SYMBOL, human_rtk$geneName)) %>%
    distinct(SYMBOL, sample, ExonicFunc.refGene) %>%
    dplyr::rename(
        mutation_type = ExonicFunc.refGene
    ) %>%
    mutate(state = 1)

plot_cnv_sum = fCNV_res_filter %>%
    filter(SYMBOL %in% c(merged_drivers$SYMBOL, human_rtk$geneName), CN>4|CN<0.5) %>%
    distinct(SYMBOL, SampleID, CNV_TYPE) %>%
    dplyr::rename(
        mutation_type = CNV_TYPE,
        sample = SampleID
    ) %>%
    mutate(state = 1)

# statistics
rbind(plot_snv_sum, plot_cnv_sum) %>%
    #filter(SYMBOL %in% human_rtk$geneName) %>%
    filter(SYMBOL %in% (has_cancers_pwy %>%
                            filter(pathwayDescript %in% pull(filter(has_cancers_pwy, geneName %in% human_rtk$geneName), 
                                                             pathwayDescript)) %>%
                            pull(geneName) %>%
                            unique()
    )) %>%
    distinct(sample)

#write.csv(rbind(plot_snv_sum, plot_cnv_sum)[,1:3], "3_results/Summary/sample_mutation_wgs.csv")

# all gene 
plot_gene_wgs_info = unique(c(plot_snv_sum$SYMBOL, plot_cnv_sum$SYMBOL))

# all mutation
plot_mutation_wgs_info = unique(c(plot_snv_sum$mutation_type, plot_cnv_sum$mutation_type))

# mutation data frame
plot_wgs_df = data.frame(
    sample = rep(plot_sample_wgs_info$SampleID, each = length(plot_gene_wgs_info)*length(plot_mutation_wgs_info)),
    SYMBOL = plot_gene_wgs_info,
    mutation_type = plot_mutation_wgs_info) %>%
    left_join(plot_cnv_sum) %>%
    left_join(plot_snv_sum, by = c("sample" = "sample", "SYMBOL" = "SYMBOL", "mutation_type" = "mutation_type"),
              suffix = c(".cnv", ".snv")) %>%
    mutate(state = ifelse(mutation_type %in% unique(plot_cnv_sum$mutation_type), state.cnv, state.snv),
           state = ifelse(is.na(state), 0, state)) %>%
    group_by(SYMBOL) %>%
    add_count() %>%
    # filter(n > 1) %>%
    dplyr::select(-state.cnv, -state.snv, -n)

# check dataframe
identical(
    (plot_wgs_df %>%
         filter(state == 1) %>%
         pull(SYMBOL) %>%
         table() %>%
         as.data.frame() %>%
         arrange(desc(Freq))),
    (rbind(plot_snv_sum, plot_cnv_sum) %>%
         pull(SYMBOL) %>%
         table() %>%
         as.data.frame() %>%
         arrange(desc(Freq)))
)

# mutation list
oncoplot_rownames = as.character((plot_wgs_df %>%
                                      filter(state == 1) %>%
                                      pull(SYMBOL) %>%
                                      table() %>%
                                      as.data.frame() %>%
                                      arrange(desc(Freq)) %>% 
                                      filter(Freq > 1) %>%
                                      pull(".")))

plot_wgs_list = oncoDataPrecess(plot_df = plot_wgs_df, plot_sample_info = plot_sample_wgs_info, plot_gene = oncoplot_rownames)
str(plot_wgs_list)

# check list
plot_wgs_df %>%
    filter(state == 1) %>%
    pull(mutation_type) %>%
    table() %>%
    as.data.frame() %>%
    arrange(desc(Freq))

lapply(plot_wgs_list, sum) %>%
    as.data.frame() %>% t() %>%
    as.data.frame() %>%
    arrange(desc(V1))

# oncoprint configure
wgs_col = custom_colors$mutation_type

alter_fun = list(
    "background" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = wgs_col["background"], col = wgs_col["background"])),
    "nonsynonymous SNV" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = wgs_col["nonsynonymous SNV"], col = NA)),
    "frameshift deletion" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.75, gp = gpar(fill = wgs_col["frameshift deletion"], col = NA)),
    "frameshift insertion" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.6, gp = gpar(fill = wgs_col["frameshift insertion"], col = NA)),
    "stopgain" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, gp = gpar(fill = wgs_col["stopgain"], col = NA)),
    "stoploss" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = wgs_col["stoploss"], col = NA)),
    "startloss" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.3, gp = gpar(fill = wgs_col["startloss"], col = NA)),
    "duplication" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2, gp = gpar(fill = wgs_col["duplication"], col = NA)),
    "deletion" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.1, gp = gpar(fill = wgs_col["duplication"], col = NA))
)

# global size option
ht_opt(
    legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
    legend_labels_gp = gpar(fontsize = 7), 
    heatmap_column_names_gp = gpar(fontsize = 8),
    heatmap_column_title_gp = gpar(fontsize = 8),
    heatmap_row_title_gp = gpar(fontsize = 8),
    merge_legends = T,
    verbose = T,
    ROW_ANNO_PADDING = unit(1.5, 'mm')
)

plot_wgs_top_annotation = HeatmapAnnotation(
    cbar = anno_oncoprint_barplot(),
    #Recurrence = plot_sample_wgs_info[match(plot_sample_wgs_info$SampleID, colnames(plot_wgs_list[[1]])),"Recurrence"]$Recurrence,
    Age = plot_sample_wgs_info[match(plot_sample_wgs_info$SampleID, colnames(plot_wgs_list[[1]])),"Age"]$Age,
    Gender = plot_sample_wgs_info[match(plot_sample_wgs_info$SampleID, colnames(plot_wgs_list[[1]])),"Sex"]$Sex,
    Grade = plot_sample_wgs_info[match(plot_sample_wgs_info$SampleID, colnames(plot_wgs_list[[1]])),"Grade"]$Grade,
    col = custom_colors,
    annotation_name_side = "left",
    annotation_name_gp = gpar(col = "grey20", fontsize = 8, fontface = "bold"), # row legend title
    simple_anno_size = unit(0.2, "cm"), # row legend height
    gap = unit(c(6,rep(2,4),6), "points") # row legend gap
)

plot_wgs_right_annotation = rowAnnotation(
    rbar = anno_oncoprint_barplot(show_fraction = F)
)

lgd_list = list( #customized legend
    Legend(labels = c("RTK", "Oncogene"), 
           title = "Gene label", type = "points", pch = 16, 
           legend_gp = gpar(col = c("#990033", "black")),
           labels_gp = gpar(fontsize = 7, col = c("#990033", "black")), # legend size
           title_gp =  gpar(fontsize = 8, fontface = "bold")) # legend title size
)

ha = oncoPrint(plot_wgs_list, alter_fun = alter_fun, col = wgs_col,
               row_names_gp = gpar(col = ifelse(oncoplot_rownames %in% human_rtk$geneName, "#990033", "black"),
                                   fontsize = 7, fontface = "italic"), #gene name
               pct_gp = gpar(fontsize = 8, col = "grey50"), pct_digits = 0, # precent
               remove_empty_columns = FALSE, 
               alter_fun_is_vectorized = T,
               pct_side = "right", row_names_side = "left",
               top_annotation = plot_wgs_top_annotation,
               right_annotation = plot_wgs_right_annotation
)


# plot
pdf("wgs_oncoplot.pdf", width = 7.3, height = 6) # inch
draw(ha, annotation_legend_list = lgd_list, merge_legend = T, legend_grouping = "original")
dev.off()
