
#' Make Correlation Plot Matrix
#'
#' @param dat Dat object with sample-level metadata contained in dat$targets and expression data contained in dat$E
#' @param features_plot_square Vector of strings naming features to be plotted along both axes. Will create a square matrix of plots: points in the lower triangle, lines in the upper triangle. Boxplots of groupwise or overall median & IQR will be included along the vertical axis. Please provide only features_plot_square OR features_plot_x and features_plot_y.
#' @param features_plot_x Vector of strings naming features to be plotted along the horizontal axis. Will create two rectangular matrices of plots: points and lines. Boxplots of groupwise or overall median & IQR will be included along the horizontal axis. Please provide only features_plot_square OR features_plot_x and features_plot_y.
#' @param features_plot_y Vector of strings naming features to be plotted along the vertical axis. Will create two rectangular matrices of plots: points and lines. Boxplots of groupwise or overall median & IQR will be included along the vertical axis. Please provide only features_plot_square OR features_plot_x and features_plot_y.
#' @param lib String indicating name of column with unique library identifier in dat$targets (should match colnames in dat$E). Default is "libID"
#' @param group_col String indicating name of column with sample groups to be plotted from dat$targets. Can be left blank.
#' @param bp_groupname_rotate T/F: Do you want to rotate the group names on the x axis of the summary boxplots? Keeps them from running into each other if they are more than a few characters long.
#' @param color_palette_custom Custom color plot named vector, names are values in "group_col" column in dat$targets. otherwise ggplot default used.
#' @param fontsize_labels font size of feature labels along sides
#' @param labels_wrapwidth wrap width of feature labels along sides
#' @param point_size How big to make the points
#'
#' @returns list of ggplot2 objects
#' @export
#'
#' @examples
#' example.voom <- kimma::example.voom
#' make_corr_plot_matrix(dat = example.voom,
#'                       features_plot_square = rownames(example.voom$genes)[c(1:5)])
#'
make_corr_plot_matrix <- function(dat = NULL,

                                  # for square plot with empty diagonal, feature list x == feature list y
                                  # for features_plot_square, will return lines in top diagonal and dots in bottom
                                  features_plot_square = NULL,

                                  # for rectangular plot of feature list x vs feature list y
                                  # for features_plot_x and features_plot_y, will return list of plots
                                  features_plot_x = NULL,
                                  features_plot_y = NULL,

                                  # column names in dat$targets, grouping column optional
                                  lib = "libID",
                                  group_col = NULL,

                                  # format boxplots
                                  bp_groupname_rotate = FALSE,

                                  # custom color plot named vector, names are values in dat$targets[[group_col]]. otherwise ggplot default used.
                                  color_palette_custom = NULL,

                                  # size and wrap width of labels along sides
                                  fontsize_labels = 25,
                                  labels_wrapwidth = 10,

                                  # how big to make the points
                                  point_size = 1
){

  group <- NULL

  reslist <- list()
  # 1. Make input data frame
  if(is.null(features_plot_square) & (is.null(features_plot_x) | is.null(features_plot_y))){
    stop("You must provide features to be plotted. This can either be a vector of features provided under 'features_plot_square' which will result in a square matrix of plots OR two vectors of features to go along the x and y axes provided under 'features_plot_x' and 'features_plot_y' respectively.")
  }
  if(!is.null(features_plot_x) & (is.null(features_plot_y) | !is.null(features_plot_square))){
    stop("You must provide features to be plotted. This can either be a vector of features provided under 'features_plot_square' which will result in a square matrix of plots OR two vectors of features to go along the x and y axes provided under 'features_plot_x' and 'features_plot_y' respectively. Make sure you have not provided both.")
  }

  if(!is.null(features_plot_square)){
    allfeats <- features_plot_square
  } else if(!is.null(features_plot_x) & !is.null(features_plot_y)){
    allfeats <- unique(c(features_plot_x, features_plot_y))
  }
  allfeats <- allfeats[which(allfeats %in% rownames(dat$E))]

  ex_df <- as.data.frame(t(dat$E)) |>
    dplyr::select(dplyr::all_of(allfeats)) |>
    tibble::rownames_to_column(var = "lib")

  plot_df <- dat$targets
  plot_df[["lib"]] <- plot_df[[lib]]

  if(!is.null(group_col)){
    plot_df[["group"]] <- plot_df[[group_col]]
  } else{
    plot_df[["group"]] <- " "
  }

  plot_df <- plot_df |>
    dplyr::left_join(ex_df)

  # 2. Make boxplots and text labels
  plist_box <- list()
  plist_labs <- list()
  for(var in allfeats){
    ymin = min(plot_df[[var]]) *0.9
    ymax = max(plot_df[[var]]) *1.1

    plist_box[[var]] <- local({
      var <- var
      ymin <- ymin
      ymax <- ymax
      plot_df <- plot_df
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = group, y = get(var))) +
        ggplot2::geom_boxplot(ggplot2::aes(fill = group), outliers = F, show.legend = FALSE) +
        ggplot2::geom_jitter(shape=16, alpha = 0.25, position=ggplot2::position_jitter(0.1))+
        ggplot2::theme_bw()+
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       legend.position = "none",
                       axis.title = ggplot2::element_blank()) +
        ggplot2::ylim(ymin,ymax)

      if(!is.null(color_palette_custom)){
        p <- p +
          ggplot2::scale_fill_manual(values = color_palette_custom)
      }
      if(bp_groupname_rotate){
        p <- p +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
      }

      p
    })

    fs <- fontsize_labels
    plist_labs[[var]] <- local({
      var <- var
      fs <- fs
      labels_wrapwidth <- labels_wrapwidth

      var_wrap <- stringr::str_wrap(var, width = labels_wrapwidth)
      grid::textGrob(var_wrap,
                     gp=grid::gpar(fontsize=fs, col="black"),
                     just = "center")
    })
  }

  # 3. Make dotplots and line plots
  if(!is.null(features_plot_square)){
    xvars <- features_plot_square[which(features_plot_square %in% rownames(dat$E))]
    yvars <- features_plot_square[which(features_plot_square %in% rownames(dat$E))]
  } else{
    xvars <- features_plot_x[which(features_plot_x %in% rownames(dat$E))]
    yvars <- features_plot_y[which(features_plot_y %in% rownames(dat$E))]
  }

  plist_corr_points <- list()
  plist_corr_lines <- list()
  for(xvar in xvars){
    xmin = min(plot_df[[xvar]]) *0.9
    xmax = max(plot_df[[xvar]]) *1.1
    for(yvar in yvars){
      if(yvar==xvar){next}

      ymin = min(plot_df[[yvar]]) *0.9
      ymax = max(plot_df[[yvar]]) *1.1

      tc <- stats::lm(formula = paste0(xvar, " ~ ", yvar), data = plot_df)
      tc2 <- summary(tc)
      rsq <- round(tc2$r.squared,3)
      clab <- paste0("Rsq = ", rsq)
      grob <- grid::grobTree(grid::textGrob(clab, x=0.1,  y=0.9, hjust=0,
                                gp=grid::gpar(col="black", fontsize=10, fontface="italic")))



      plist_corr_lines[[paste0(xvar, "_",yvar)]] <- local({
        xvar <- xvar
        yvar <- yvar
        plot_df <- plot_df
        grob <- grob
        xmin <- xmin
        xmax <- xmax
        ymin <- ymin
        ymax <- ymax

        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = get(xvar), y = get(yvar))) +
          ggplot2::geom_smooth(method = "lm", color = "black", se = T)+
          ggplot2::labs(x = xvar,
               y = yvar) +
          ggplot2::theme_bw() +
          ggplot2::theme(panel.grid = ggplot2::element_blank(),
                legend.position = "none",
                axis.title = ggplot2::element_blank()) +
          ggplot2::annotation_custom(grob) +
          ggplot2::xlim(xmin,xmax) +
          ggplot2::ylim(ymin,ymax)
        if(!is.null(group_col)){
          p <- p +
            ggplot2::geom_smooth(method = "lm",  ggplot2::aes(color = group, fill = group), se = T)
        }
        if(!is.null(color_palette_custom)){
          p <- p +
            ggplot2::scale_color_manual(values = color_palette_custom) +
            ggplot2::scale_fill_manual(values = color_palette_custom)
        }
        p
      })

      plist_corr_points[[paste0(xvar, "_",yvar)]] <- local({
        xvar <- xvar
        yvar <- yvar
        plot_df <- plot_df
        grob <- grob
        xmin <- xmin
        xmax <- xmax
        ymin <- ymin
        ymax <- ymax

        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = get(xvar), y = get(yvar))) +
          ggplot2::labs(x = xvar,
               y = yvar) +
          ggplot2::theme_bw() +
          ggplot2::theme(panel.grid = ggplot2::element_blank(),
                legend.position = "none",
                axis.title = ggplot2::element_blank()) +
          ggplot2::annotation_custom(grob) +
          ggplot2::xlim(xmin,xmax) +
          ggplot2::ylim(ymin,ymax)
        if(!is.null(group_col)){
          p <- p +
            ggplot2::geom_point(ggplot2::aes(fill = as.factor(group)), color = "black", shape = 21, size = point_size, stroke = 0.1) +
            ggplot2::geom_smooth(method = "lm", color = "black", se = F)
        } else{
          p <- p +
            ggplot2::geom_point(fill = "grey70", color = "black", shape = 21, size = point_size, stroke = 0.1) +
            ggplot2::geom_smooth(method = "lm", color = "black", se = F)
        }
        if(!is.null(color_palette_custom)){
          p <- p +
            ggplot2::scale_fill_manual(values = color_palette_custom)
        }
        p
      })
    }
  }

  # make layout for square implementation
  if(!is.null(features_plot_square)){
    # 1. Define your module names (dynamic)
    mod_names <- unique(allfeats)
    n_mods <- length(mod_names)
    # 2. Initialize an empty list to store all grid elements
    plist_dynamic <- list()

    # 3. Build the grid row by row
    for(i in 1:n_mods) {
      row_mod <- mod_names[i]

      # Column 1: Row Label (e.g., "m1")
      plist_dynamic[[length(plist_dynamic) + 1]] <- plist_labs[[row_mod]]

      # Column 2: The Boxplot for that module
      plist_dynamic[[length(plist_dynamic) + 1]] <- plist_box[[row_mod]]

      # Columns 3 to (n_mods + 2): The Correlation plots
      for(j in 1:n_mods) {
        col_mod <- mod_names[j]

        if (j < i) {
          # Lower Triangle: Points
          plist_dynamic[[length(plist_dynamic) + 1]] <- plist_corr_points[[paste0(col_mod, "_", row_mod)]]
        } else if (j == i) {
          # Diagonal: Spacer
          plist_dynamic[[length(plist_dynamic) + 1]] <- patchwork::plot_spacer()
        } else {
          # Upper Triangle: Regression Lines
          plist_dynamic[[length(plist_dynamic) + 1]] <- plist_corr_lines[[paste0(col_mod, "_", row_mod)]]
        }
      }
    }
    # 4. Add the Bottom Label Row
    plist_dynamic[[length(plist_dynamic) + 1]] <- patchwork::plot_spacer() # Spacer for Row Label col
    plist_dynamic[[length(plist_dynamic) + 1]] <- patchwork::plot_spacer() # Spacer for Boxplot col

    for(j in 1:n_mods) {
      col_mod <- mod_names[j]
      plist_dynamic[[length(plist_dynamic) + 1]] <- plist_labs[[col_mod]]
    }
    # 5. Wrap it all up
    # Total columns = Number of modules + Label column + Boxplot column
    total_cols <- n_mods + 2

    final_plot <- patchwork::wrap_plots(plist_dynamic, ncol = total_cols) +
      patchwork::plot_layout(guides = "collect", axes = "collect")

    reslist[["square"]] <- final_plot
  }else{
    # 1. Define your module names (dynamic)
    mod_names_x <- features_plot_x[which(features_plot_x %in% rownames(dat$E))]
    mod_names_y <- features_plot_y[which(features_plot_y %in% rownames(dat$E))]
    n_mods_x <- length(mod_names_x)
    n_mods_y <- length(mod_names_y)
    # 2. Initialize an empty list to store all grid elements
    plist_dynamic_points <- list()
    plist_dynamic_lines <- list()

    # 3. Build the grid row by row
    for(i in 1:n_mods_y) {
      row_mod <- mod_names_y[i]

      # Column 1: Row Label (e.g., "m1")
      plist_dynamic_points[[length(plist_dynamic_points) + 1]] <- plist_labs[[row_mod]]
      plist_dynamic_lines[[length(plist_dynamic_lines) + 1]] <- plist_labs[[row_mod]]

      # Column 2: The Boxplot for that module
      plist_dynamic_points[[length(plist_dynamic_points) + 1]] <- plist_box[[row_mod]]
      plist_dynamic_lines[[length(plist_dynamic_lines) + 1]] <- plist_box[[row_mod]]

      # Columns 3 to (n_mods + 2): The Correlation plots
      for(j in 1:n_mods_x) {
        col_mod <- mod_names_x[j]

        plist_dynamic_points[[length(plist_dynamic_points) + 1]] <- plist_corr_points[[paste0(col_mod, "_", row_mod)]]
        plist_dynamic_lines[[length(plist_dynamic_lines) + 1]] <- plist_corr_lines[[paste0(col_mod, "_", row_mod)]]
      }
    }



    # 4. Add the Bottom module boxplot
    plist_dynamic_points[[length(plist_dynamic_points) + 1]] <- patchwork::plot_spacer() # Spacer for Row Label col
    plist_dynamic_points[[length(plist_dynamic_points) + 1]] <- patchwork::plot_spacer() # Spacer for Boxplot col
    plist_dynamic_lines[[length(plist_dynamic_lines) + 1]] <- patchwork::plot_spacer() # Spacer for Row Label col
    plist_dynamic_lines[[length(plist_dynamic_lines) + 1]] <- patchwork::plot_spacer() # Spacer for Boxplot col
    for(j in 1:n_mods_x){
      col_mod <- mod_names_x[j]
      plist_dynamic_points[[length(plist_dynamic_points) + 1]] <- plist_box[[col_mod]]
      plist_dynamic_lines[[length(plist_dynamic_lines) + 1]] <- plist_box[[col_mod]]
    }

    # 5. Add the Bottom Label Row
    plist_dynamic_points[[length(plist_dynamic_points) + 1]] <- patchwork::plot_spacer() # Spacer for Row Label col
    plist_dynamic_points[[length(plist_dynamic_points) + 1]] <- patchwork::plot_spacer() # Spacer for Boxplot col
    plist_dynamic_lines[[length(plist_dynamic_lines) + 1]] <- patchwork::plot_spacer() # Spacer for Row Label col
    plist_dynamic_lines[[length(plist_dynamic_lines) + 1]] <- patchwork::plot_spacer() # Spacer for Boxplot col
    for(j in 1:n_mods_x){
      col_mod <- mod_names_x[j]
      plist_dynamic_points[[length(plist_dynamic_points) + 1]] <- plist_labs[[col_mod]]
      plist_dynamic_lines[[length(plist_dynamic_lines) + 1]] <- plist_labs[[col_mod]]
    }

    # 6. Wrap it all up
    # Total columns = Number of x modules + Label column + Boxplot column
    total_cols <- n_mods_x + 2
    final_plot_points <- patchwork::wrap_plots(plist_dynamic_points, ncol = total_cols) +
      patchwork::plot_layout(guides = "collect", axes = "collect")
    final_plot_lines <- patchwork::wrap_plots(plist_dynamic_lines, ncol = total_cols) +
      patchwork::plot_layout(guides = "collect", axes = "collect")

    reslist[["points"]] <- final_plot_points
    reslist[["lines"]] <- final_plot_lines
  }
  return(reslist)
}
