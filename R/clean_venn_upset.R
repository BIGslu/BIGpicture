#' Format gene lists for venn or upset plots
#'
#' @return Data frame
#' @keywords internal

clean_venn_upset <- function(model_result, models,
                             variables, contrasts,
                             intercept, random){
  model <- variable <- contrast_ref <- contrast_lvl <- label <- dataset <- NULL

  #### Extract results ####
  dat_all <- data.frame()
  for(i in 1:length(model_result)){
    dat_temp <- model_result[[i]]

    #list all model df
    m <- names(dat_temp)[!grepl(".fit", names(dat_temp))]
    m <- m[!grepl(".error", m)]

    dat <- data.frame()
    for(j in m){
      dat <- dat_temp[[j]] %>%
        dplyr::mutate(dataset= names(model_result)[i]) %>%
        dplyr::bind_rows(dat)
    }

    dat_all <- dplyr::bind_rows(dat, dat_all)
  }

  #Subset to models of interest if provided
  if(!is.null(models)){
    dat_subset <- dat_all %>%
      dplyr::filter(model %in% models)
  } else {
    dat_subset <- dat_all
  }

  #### List variables of interest ####
  if(!is.null(variables)){
    var_all <- variables
  } else {
    var_all <- unique(dat_subset$variable)
  }

  #remove intercept if specified
  if(!intercept){
    var_all <- var_all[var_all != "(Intercept)"]
  }

  #remove random effects if specified
  if(!random){
    var_all <- var_all[!grepl("\\|", var_all)]
  }

  ##### List contrasts of interest ####
  if(any(grepl("contrast", unique(dat_subset$model)))){
    if(is.null(contrasts)){
      con_filter <- dplyr::distinct(dat_subset, variable, contrast_ref, contrast_lvl)
      if(!is.null(variables)){
        con_filter <- dplyr::filter(con_filter, variable %in% variables) %>%
          dplyr::select(-variable)
      }
    } else {
      con_filter <- strsplit(contrasts, split=" - ") %>%
        as.data.frame() %>% t() %>%  as.data.frame()
      colnames(con_filter) <- c("contrast_lvl", "contrast_ref")
      rownames(con_filter) <- NULL
    }
  } else {
    con_filter <- NULL
  }

  #### filter data to variables/contrasts of interest ####
  if(any(grepl("contrast", unique(dat_subset$model)))){
    dat_filter <- dat_subset %>%
      dplyr::filter(variable %in% var_all) %>%
      #add dataset name to labels
      dplyr::inner_join(con_filter, by = c("contrast_ref", "contrast_lvl")) %>%
      dplyr::mutate(label = paste(contrast_lvl, contrast_ref, sep="\n-\n"))
  } else if(length(model_result) > 1) {
    dat_filter <- dat_subset %>%
      dplyr::filter(variable %in% var_all) %>%
      #add dataset name to labels
      dplyr::mutate(label = variable)
  } else {
    dat_filter <- dat_subset %>%
      dplyr::filter(variable %in% var_all) %>%
      dplyr::mutate(label = variable)
  }

  #Add data set and/or model name to label
  if(length(unique(dat_filter$model)) > 1){
    dat_filter <- dat_filter %>%
      dplyr::mutate(label = paste(model, label, sep="\n"))
  }

  if(length(unique(dat_filter$dataset)) > 1 &
     !identical(sort(unique(dat_filter$dataset)),
                sort(unique(dat_filter$model)))){
    dat_filter <- dat_filter %>%
      dplyr::mutate(label = paste(dataset, label, sep="\n"))
  }

  dat_filter_all <- list()
  dat_filter_all[["dat"]] <- dat_filter
  dat_filter_all[["con"]] <- con_filter
  return(dat_filter_all)
}
