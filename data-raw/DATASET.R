library(kimma)
library(readr)
library(dplyr)

#Linear models
model_result <- kmFit(dat = example.voom, kin = example.kin, patientID = "donorID",
                      run.lmekin = TRUE, run.lme = TRUE, metrics=TRUE,
                      model = "~ virus*asthma + (1|donorID)", processors = 6)

usethis::use_data(model_result, overwrite = TRUE)

#Enrichment colors for STRING
enrichment <- read_csv("data-raw/example.enrich.csv") %>%
  mutate(SYMBOLs = strsplit(as.character(SYMBOLs), "/"))
usethis::use_data(enrichment, overwrite = TRUE)
