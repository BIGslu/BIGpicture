#' Map genes to string network
#'
#' To be used in conjunction with plot_string
#'
#' @param genes Character vector of HGNC symbols or ENSEMBL ID of genes
#' @param version Numeric version of STRING data base. Default 11.5 which is current as of 2022-03-30
#' @param score_threshold Numeric combined score threshold for edge inclusion in network. Default it 400 which is STRING defined "medium". 0 to 1000 allowed
#' @param species Character string of "human" or "mouse"
#'
#' @return List containing "subgraph" and "map"
#' @export
#'
#' @examples
#' #See plot_string for full example

map_string <- function(genes, version=11.5, score_threshold=400,
                       species="human"){
  #Set species ID
  if(species == "human"){
    speciesID <- 9606
  } else if(species == "mouse"){
    speciesID <- 10090
  } else{
   stop("Only human and mouse species supported.")
  }

  #STRING database
  string_db <- STRINGdb::STRINGdb$new(version=as.character(version),
                                      species=speciesID,
                                      score_threshold=score_threshold,
                                      input_directory="")

  #Format gene vector to matrix
  genes.mat <- as.matrix(genes)
  colnames(genes.mat) <- "gene"

  #Map genes to STRING
  map <- string_db$map(genes.mat, "gene", removeUnmappedRows = TRUE)

  #Fill in NA with 0
  map2 <- dplyr::mutate_if(map, is.numeric, ~replace_na(., 0))

  # Create igraph object
  subgraph <- string_db$get_subnetwork(map2$STRING_id)

  result.ls <- list()
  result.ls[["map"]] <- map2
  result.ls[["subgraph"]] <- subgraph
  return(result.ls)
}
