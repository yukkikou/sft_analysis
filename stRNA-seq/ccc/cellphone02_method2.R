# environment
workdir = "/media/data5/hxy/SFT/RNAseq/4_Spatial/CCC/cellphoneDB"
setwd(workdir)

# Get command line arguments
args = commandArgs(trailingOnly = TRUE)
# Assign arguments to variables
#sample_prefix = "A-2401158-04"
sample_prefix = args[1]

resdir = paste0(sample_prefix, "/results/method2/")
setwd(resdir)

# packages
library(tidyverse)

# functions
SigCellFrom = function(cell_type, mean_thred = 0, inter_sc = 50, pct = 0.1, pval = 0.05){

    message("*** processing means ***")
    means_matrix = means_results %>%
        column_to_rownames("id_cp_interaction") %>%
        dplyr::select(starts_with(paste0(cell_type, "|"))) 
    
    means_idx = rownames(means_matrix[rowSums(means_matrix) > mean_thred,])
   
    msig_matrix = means_sig %>%
        column_to_rownames("id_cp_interaction") %>%
        dplyr::select(starts_with(paste0(cell_type, "|"))) 
    
    msig_idx = rownames(msig_matrix[rowSums(msig_matrix, na.rm = T) > 0,])
 

    message("*** processing intersection ***")
    scores_matrix = intersection_scores %>%
        column_to_rownames("id_cp_interaction") %>%
        dplyr::select(starts_with(paste0(cell_type, "|"))) 
    
    scores_idx = rownames(scores_matrix[rowSums(scores_matrix > inter_sc) > 0,])

    pvalue_matrix = pvalue_results %>%
        column_to_rownames("id_cp_interaction") %>%
        dplyr::select(starts_with(paste0(cell_type, "|"))) 
    
    pvalue_idx = rownames(pvalue_matrix[rowSums(pvalue_matrix < pval) > 0,])
 
    keep_idx = Reduce(intersect, list(means_idx, msig_idx,scores_idx, pvalue_idx))

    message("*** keep id number: ", length(keep_idx))
    
    message("*** generating final results ***")

    filtered_interactions = means_results %>%
        filter(id_cp_interaction %in% keep_idx) %>%
        inner_join(intersection_scores, by = col_keep, suffix = c(".means", ".score")) %>%
        inner_join(pvalue_results, by = col_keep) %>%
        inner_join(means_sig, by = col_keep, suffix= c(".pvalue",".msig")) %>%
        dplyr::select(all_of(col_keep), starts_with(paste0(cell_type, "|")))
     
    final_interactions = deconvoluted_pct %>%
        dplyr::select(1:7, any_of(cell_type)) %>%
        filter(if_any(all_of(cell_type), ~ . > pct)) %>%
        inner_join(filtered_interactions) %>%
        mutate(key_cell = cell_type)
    
    return(final_interactions)
}



# data
intersection_scores = read_tsv(paste0("statistical_analysis_interaction_scores_", sample_prefix, ".txt"))
means_results = read_tsv(paste0("statistical_analysis_means_",sample_prefix,".txt"))
means_sig =  read_tsv(paste0("statistical_analysis_significant_means_",sample_prefix,".txt"))
pvalue_results =  read_tsv(paste0("statistical_analysis_pvalues_", sample_prefix, ".txt"))
deconvoluted_results = read_tsv(paste0("statistical_analysis_deconvoluted_",sample_prefix,".txt"))
deconvoluted_pct = read_tsv(paste0("statistical_analysis_deconvoluted_percents_",sample_prefix,".txt"))


col_keep = colnames(intersection_scores)[!str_detect(colnames(intersection_scores), "\\|")]

# config
cell_types = read_tsv("../../data/microenvironment.tsv")$cell_type
mean_thred = 0
inter_sc = 30
pct = 0.01
pval = 0.1

final_res = map_dfr(cell_types, ~ SigCellFrom(.x, mean_thred = mean_thred, inter_sc = inter_sc, pct = pct, pval = pval))

final_res = final_res %>% 
    pivot_longer(cols = contains("|"), names_to = "value_type", values_to = "value") %>%
    filter(!is.na(value)) %>%
    separate(value_type, sep = "\\.", into = c("interaction","value_type")) %>%
    pivot_wider(names_from = "value_type",values_from ="value") %>%
    filter(means > mean_thred, (score > inter_sc)|(pvalue < pval), !is.na(msig)) %>% 
    dplyr::select(key_cell,interaction,interacting_pair,means,score,pvalue,msig,classification,directionality,everything()) %>%
    distinct() 


# Save the filtered interactions to a file
write.csv(final_res, paste0(sample_prefix, "_significant_interactions.csv"), row.names = FALSE)

