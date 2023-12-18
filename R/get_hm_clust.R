#' Get data frame of genes per cluster from ComplexHeatmap output
#'
#' @param dat Ordered and scaled matrix of gene counts
#' @param hm ComplexHeatMap output specified with draw(h)
#' @param dimension Row or column clusters
#' @return Data frame of genes per row/column, which row or column they owe to, and which cluster number
#'
#' @examples
#' heatmap_cluster_genes <- get_hm_clust(dat = dat, hm = draw(hm), dimension = "row")

get_hm_clust <- function(dat, hm, dimension){

  cluster.result <- data.frame()

  # Rows of heatmap
  if(dimension == "row"){
    #Deal with single cluster results
    if(is.list(row_order(hm))){
      clust.tot <- length(row_order(hm))
    } else {
      clust.tot <- 1
    }

    for (i in 1:clust.tot){
      #Get row indices
      if(is.list(row_order(hm))){
        # version of hm row order needs to be numeric (CHM v2.15.1)
        cluster.index <- row_order(hm)[[(i)]]
      } else {
        cluster.index <- row_order(hm)
      }

      #Clusters with >1 element
      if(length(cluster.index) > 1){
        #Pull clusters
        cluster.result <- t(t(row.names(dat[cluster.index,]))) %>%
          #Convert to data frame
          as.data.frame() %>%
          #row order within cluster
          mutate(row_within_cluster = 1:length(cluster.index)) %>%
          #add cluster name
          mutate(cluster = paste0("cluster", i)) %>%
          #Rename default column
          dplyr::rename(row=V1) %>%
          #concatenate results
          bind_rows(cluster.result)
      } else {
        #Clusters with only 1 element
        cluster.result <- data.frame(row = rownames(dat)[cluster.index],
                                     row_within_cluster = 1,
                                     cluster = paste0("cluster", i)) %>%
          bind_rows(cluster.result)
      }
    }
  } else if(dimension == "col"){
    # Columns of heatmap

    #Deal with single cluster results
    if(is.list(column_order(hm))){
      clust.tot <- length(column_order(hm))
    } else {
      clust.tot <- 1
    }

    for (i in 1:clust.tot){
      #Get column indices
      if(is.list(column_order(hm))){
        cluster.index <- column_order(hm)[[(i)]]
      } else {
        cluster.index <- column_order(hm)
      }

      #Clusters with >1 element
      if(length(cluster.index) > 1){
        cluster.result <- t(t(colnames(dat[,cluster.index]))) %>%
          as.data.frame() %>%
          mutate(col_within_cluster = 1:length(cluster.index)) %>%
          mutate(cluster = paste0("cluster", i)) %>%
          dplyr::rename(col=V1) %>%
          bind_rows(cluster.result)
      } else {
        #Clusters with only 1 element
        cluster.result <- data.frame(col = colnames(dat)[cluster.index],
                                     col_within_cluster = 1,
                                     cluster = paste0("cluster", i)) %>%
          bind_rows(cluster.result)
      }
    }
  } else{ stop("dimension must be one of row or col.") }

  cluster.result <- cluster.result[order(cluster.result$cluster),]

  return(cluster.result)
}
