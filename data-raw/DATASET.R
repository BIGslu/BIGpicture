library(kimma)

#Linear models
model_result <- kmFit(dat = example.voom, kin = example.kin,
                      run.lmerel = TRUE, run.lme = TRUE,
                      model = "~ virus*asthma + (1|ptID)", processors = 6)

usethis::use_data(model_result, overwrite = TRUE)
