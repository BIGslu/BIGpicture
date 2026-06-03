# library(kimma)
# library(tidyverse)
#
# dat <- kimma::example.voom
#
# model <- kmFit(dat = dat,
#                model = "~ virus*asthma + (1|ptID)",
#                run_lme = TRUE,
#                use_weights = TRUE,
#                metrics = TRUE,
#                run_contrast = TRUE)
#
# gene_list <- model$lme.contrast %>%
#   filter(variable == "virus*asthma" & FDR < 0.05) %>%
#   filter(contrast_ref == "none asthma" & contrast_lvl == "HRV asthma") %>%
#   pull(gene)
#
# # VARS:
# ### gene_list: vector of genes to plot in heatmap
# ### dat: voom object
# ### vars: variables to color columns by
#
# IL33_stim <- final_model$lme.contrast %>%
#   filter(variable == "treatment*stimulation" & FDR < 0.2) %>%
#   filter(grepl("IL-33",contrast_ref) & grepl("IL-33",contrast_lvl)) %>%
#   arrange(FDR, contrast_ref,contrast_lvl) %>%
#   left_join(ensembl_key, by = c("gene"="geneName")) %>%
#   mutate(higher_in=case_when(estimate > 0 ~ contrast_lvl,
#                              estimate < 0 ~ contrast_ref)) %>%
#   select(variable,contrast_ref,contrast_lvl,gene,hgnc_symbol,estimate,higher_in,FDR)
#
# IL33_stim_genes <- IL33_stim %>% pull(gene)
#
# IL33_meta <- as.data.frame(dat$targets) %>%
#   filter(stimulation == "IL-33") %>%
#   mutate(treat.stim=paste0(treatment,"_",stimulation)) %>%
#   select(libid,treatment,stimulation,donorId,AECdonor) %>%
#   arrange(treatment) %>%
#   mutate(treatment = factor(treatment,levels=c("None","AEC")))
#
# IL33_meta %>%
#   mutate(group = paste(donorId,treatment,stimulation,AECdonor,sep="_")) %>%
#   group_by(group) %>%
#   count() %>%
#   knitr::kable()
#
# IL33_dat <- dat$E[IL33_stim_genes,unique(IL33_meta$libid)]
# rowlabels <- (ensembl_key %>% column_to_rownames("geneName"))[IL33_stim_genes,]
#
# ann_colors <- list(treatment = c(None = "#B52AC7", AEC = "#FAA916"))
# column_ha = HeatmapAnnotation(treatment = IL33_meta$treatment,
#                               col = ann_colors)
#
# scaled_dat = t(scale(t(IL33_dat)))
# Heatmap(scaled_dat,show_row_names = T,show_column_names = F,
#         cluster_rows = T,cluster_columns = F,
#         clustering_distance_rows = "euclidean",
#         column_split = IL33_meta$treatment,#cluster_column_slices = T,
#         top_annotation = column_ha,
#         heatmap_legend_param = list(title = "row scaled expr"),
#         use_raster = FALSE) +
#   rowAnnotation(labels = anno_text(rowlabels,which="row"))
