library(kimma)
library(readr)
library(dplyr)
library(SEARchways)

#Linear models
example_model <- kmFit(dat = example.voom, kin = example.kin,
                      run.lmerel = TRUE, run.lme = TRUE, metrics=TRUE,
                      model = "~ virus*asthma + (1|ptID)", processors = 6)

usethis::use_data(example_model, overwrite = TRUE)

#### Done within function examples. Could make data if too slow when package grows ####
#Enrichment colors for STRING
# genes.OI <- example_model$lmerel %>%
#   filter(variable == "virus" & FDR < 0.05) %>%
#   distinct(variable, hgnc_symbol)
#
# example_enrich <- BIGprofiler(gene_df = genes.OI, category = "H")
#
# usethis::use_data(example_enrich, overwrite = TRUE)

# GSEA results
# genes.FC <- example_model$lmerel %>%
#   filter(variable == "virus") %>%
#   distinct(variable, hgnc_symbol, estimate)
#
# example_gsea <- BIGprofiler(gene_df = genes.OI, category = "H")
#
# usethis::use_data(example_enrich, overwrite = TRUE)
