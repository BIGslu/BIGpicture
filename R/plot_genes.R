#' Plot boxplots and dotplots for genes and variables of interest
#'
#' @param dat DGEList or EList object with expression data (counts, E), sample metadata (samples, targets), and gene annotation (genes)
#' @param counts If dat not provided. Data frame of gene expression. Genes are rows, samples are columns. geneID must be rownames or first column
#' @param meta If dat not provided. Data frame of sample meta data with samples as rows.
#' @param genes Optional if dat not provided. Data frame of gene annotation with genes as rows.
#' @param fdr Optional. Model results output by kmFit( ). Used to label significant differences
#' @param libraryID Character string of variable name to use when combining expression and sample data
#' @param geneID Character string of variable name to use when combining expression and gene data
#' @param subset_genes Optional. Character vector of genes to plot. Must match names in geneID column. If not provided, all genes are plotted.
#' @param variables Character vector of variable names to include in plot. Variables can be character, factor, or numeric. One two-variable interaction term allowed
#' @param colorID Optional. Character string for variable to color point by
#' @param processors Numeric processors to run in parallel. Default is 2 less than the total available
#'
#' @param subset.genes Deprecated form of subset_genes
#'
#' @importFrom foreach %dopar%
#' @return List of ggplot objects
#' @export
#'
#' @examples
#' #Data from kimma
#' example.voom <- kimma::example.voom
#'
#' subset_genes <- c("ENSG00000250479","ENSG00000250510","ENSG00000255823")
#'
#' plot_genes(dat = example.voom, fdr = example.model$lmerel,
#'      subset_genes = subset_genes, geneID="geneName",
#'      variables = c("virus*asthma", "lib.size"), colorID = "virus")

plot_genes <- function(dat=NULL, counts=NULL, meta=NULL, genes=NULL,
                       fdr=NULL,
                       libraryID="libID", geneID="geneName",
                       subset_genes=NULL,
                       variables, colorID=NULL,
                       processors=NULL,
                       subset.genes = NULL) {
  i <- gene <- E <- model <- variable <- contrast_ref <- contrast_lvl <- estimate <- pval <- FDR <- NULL

  # backwards compatibility
  if(!is.null(subset.genes)){subset_genes <- subset.genes}

  ###### Parallel ######
  #setup parallel processors
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    #Use 2 in CRAN/Travis/AppVeyor
    processors.to.use <- 2
  } else if (is.null(processors)){
    #Use 2 less than total if not user defined
    processors.to.use <- parallel::detectCores()-2
    if(processors.to.use == 0){
      stop("Error processors: Default resulted in 0. Please correct.")}
  } else {
    #Use user defined number
    processors.to.use <- processors
  }

  cl <- parallel::makeCluster(processors.to.use)

  #Set seed
  set.seed(4389)

  ########## Check data ##########
  if(!is.null(dat) & !is.null(counts)){ stop("Please provide only one of dat or counts.")}
  if(!is.null(counts) & is.null(meta)){ stop("When using counts, meta must also be provided.")}

  ########## Load data ##########
  if(!is.null(dat) & !is.null(counts)){ stop("Please provide only one of dat or counts.")}
  if(!is.null(counts) & is.null(meta)){ stop("When using counts, meta must also be provided.")}

  if(!is.null(dat)){
    if(inherits(dat, "DGEList")){
      dat.counts <- as.data.frame(dat$counts) %>%
        tibble::rownames_to_column(geneID)
      dat.meta <- as.data.frame(dat$samples)
      dat.genes <- as.data.frame(dat$genes)
    } else if(inherits(dat, "EList")){
      dat.counts <- as.data.frame(dat$E) %>%
        tibble::rownames_to_column(geneID)
      dat.meta <- as.data.frame(dat$targets)
      dat.genes <- as.data.frame(dat$genes)
    } else {
      stop("dat must be a DGEList or EList object.")
    }
  }
  if(!is.null(counts)){
    #Move rownames, if exist
    if(is.numeric(as.matrix(counts))){
      dat.counts <- as.data.frame(counts) %>%
        tibble::rownames_to_column(geneID)
    } else {
      dat.counts <- as.data.frame(counts)
      colnames(dat.counts)[1] <- geneID
    }}
  if(!is.null(meta)){
    dat.meta <- as.data.frame(meta)
  }
  if(!is.null(genes)){
    dat.genes <- as.data.frame(genes)
  }

  ########## Combine data ##########
  #find match variable
  match_var <- which(colSums(dat.counts[[geneID]][1] == dat.genes, na.rm = TRUE) > 0)
  match_var <- colnames(dat.genes)[match_var]
  colnames(dat.genes)[colnames(dat.genes)==match_var] <- geneID

  if(!is.null(subset_genes)){
    plot.dat <- dat.counts %>%
      dplyr::filter(get(geneID) %in% subset_genes) %>%
      tidyr::pivot_longer(-geneID, names_to = libraryID,
                          values_to = "E") %>%
      dplyr::left_join(dat.meta, by=libraryID) %>%
      dplyr::left_join(dat.genes, by=geneID)
  } else{
    plot.dat <- dat.counts %>%
      tidyr::pivot_longer(-geneID, names_to = libraryID,
                          values_to = "E") %>%
      dplyr::left_join(dat.meta, by=libraryID) %>%
      dplyr::left_join(dat.genes, by=geneID)
  }

  ########## Create interaction variable if selected ##########
  var.i <- variables[grepl(":|[*]", variables)]
  if(length(var.i)>0){
    vars.i <- strsplit(var.i, split=":|[*]")[[1]]
    plot.dat <- plot.dat %>%
      dplyr::mutate(interaction = paste(get(vars.i[1]), get(vars.i[2]), sep="\n"))
    variables.noI <- unlist(strsplit(variables, split=":|[*]"))
  } else {
    variables.noI <- variables
  }

  ########## Subset genes ##########
  #List all genes/modules
  if(!is.null(subset_genes)){
    to_plot <- subset_genes
  } else{
    #list all genes
    to_plot <- sort(unique(unlist(plot.dat[,geneID])))
  }

  # Setup parallel computing
  doParallel::registerDoParallel(cl)

  ##########  Loop through genes ##########
  plot.ls <- list()
  plot.ls <- foreach::foreach(i = 1:length(to_plot),
                              .packages = c("dplyr","magrittr","tibble","ggplot2","ggpubr",
                                            "forcats",
                                            "patchwork","foreach","doParallel")) %dopar% {

    #Subset data to gene/module of interest
    plot.dat.sub <- plot.dat %>%
      dplyr::filter(get(geneID) == to_plot[i])

    ########## Loop through variables ##########
    gene.plot.ls = list()

    #fdr table
    if(!is.null(fdr) & "contrast_ref" %in% colnames(fdr)){
      dat.fdr.sub <- dplyr::filter(fdr, gene==to_plot[i]) %>%
        dplyr::select(model, gene, variable, contrast_ref,
                      contrast_lvl, estimate, pval, FDR)
      fdr.plot <- ggpubr::ggtexttable(dat.fdr.sub, rows = NULL)
    } else if(!is.null(fdr)){
      dat.fdr.sub <- dplyr::filter(fdr, gene==to_plot[i]) %>%
        dplyr::select(model, gene, variable, estimate, pval, FDR)
      fdr.plot <- ggpubr::ggtexttable(dat.fdr.sub, rows = NULL)
    } else{ fdr.plot = NULL }

    for(j in 1:length(variables.noI)){
      #Type = factor or character variables
      if(is.factor(plot.dat.sub[[variables.noI[j]]]) |
         is.character(plot.dat.sub[[variables.noI[j]]])) {

        plot1 <- plot.dat.sub %>%
          ggplot2::ggplot(ggplot2::aes_string(x=variables.noI[j], y="E")) +
          ggplot2::geom_boxplot(outlier.shape = NA) +
          ggplot2::geom_jitter(ggplot2::aes_string(color=colorID), height=0, width=0.2) +
          ggplot2::theme_classic() +
          ggplot2::stat_summary(fun.y=mean, geom="point",
                                shape="square", color="black", size=3)

        #Color legend if not x variable
        if(!is.null(colorID)){
        if(variables.noI[j] != colorID){
          gene.plot.ls[[variables.noI[j]]] <- plot1
        }} else {
          gene.plot.ls[[variables.noI[j]]] <- plot1 +
            ggplot2::theme(legend.position = "none")
        }
      } else
        #Type = numeric
        if(is.numeric(plot.dat.sub[[variables.noI[j]]])){
          plot1 <- plot.dat.sub %>%
            ggplot2::ggplot(ggplot2::aes_string(x=variables.noI[j], y="E")) +
            ggplot2::geom_point(ggplot2::aes_string(color=colorID)) +
            ggplot2::theme_classic() +
            ggplot2::geom_smooth(method='lm', formula= y~x, color="black",
                                 se=FALSE)

          #Color legend if not x variable
          if(!is.null(colorID)){
            if(variables.noI[j] != colorID){
            gene.plot.ls[[variables.noI[j]]] <- plot1
          } } else {
            gene.plot.ls[[variables.noI[j]]] <- plot1 +
              ggplot2::theme(legend.position = "none")
          }

        } else{
          stop("Variables of interest must be numeric, character, or factor.")
        }
    }

    #Interaction plot
    if(length(var.i)>0){
      plot3 <- plot.dat.sub %>%
        ggplot2::ggplot(ggplot2::aes(x=interaction, y=E)) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::geom_jitter(ggplot2::aes_string(color=colorID), height=0, width=0.2) +
        ggplot2::theme_classic() +
        ggplot2::labs(x="", color=colorID) +
        ggplot2::theme(legend.position = "right") +
        ggplot2::stat_summary(fun.y=mean, geom="point",
                              shape="square", color="black", size=3)

    } else{
      plot3 <- NULL
    }

    #### Combine plots ####
    #Main plot title
    if (!is.null(dat.genes)){
      symbol <- dat.genes %>%
        dplyr::filter(get(geneID) == to_plot[i]) %>%
        dplyr::select(dplyr::ends_with("symbol")) %>%
        unlist()

      title <- paste(to_plot[i], unique(symbol), sep=" ", collapse=" ")
    } else{
      title <- to_plot[i]
    }

    plot_row1 <- patchwork::wrap_plots(gene.plot.ls, nrow = 1)

    if(!is.null(plot3) & !is.null(fdr.plot)){
      plot_final <- patchwork::wrap_plots(fdr.plot, plot_row1, plot3, nrow=3) +
        patchwork::plot_annotation(title = title)
    } else if(!is.null(plot3)){
      plot_final <- patchwork::wrap_plots(plot_row1, plot3, nrow=2) +
        patchwork::plot_annotation(title = title)
    } else if(!is.null(fdr.plot)){
      plot_final <- patchwork::wrap_plots(fdr.plot, plot_row1, nrow=2) +
        patchwork::plot_annotation(title = title)
    } else {
      plot_final <- plot_row1
    }

    #### Save to list
    plot.ls[[to_plot[i]]] <- plot_final
  }
  parallel::stopCluster(cl)

  names(plot.ls) <- to_plot
  return(plot.ls)
}
