calc_com_louv <- function(df_sim) {
  df_sim <- dplyr::rename(df_sim, weight= 3 )
  G <- igraph::graph_from_data_frame(df_sim, directed=FALSE)
  louv_com <- igraph::cluster_louvain(G)
  output <- list(community = louv_com, G = G)
}

create_reduced_bipartite_network <- function(freq_df, comms_1, comms_2){
  '%>%' <- magrittr::'%>%'
  all_of <- tidyselect::all_of
  as_tibble <- tibble::as_tibble
  pull <- dplyr::pull
  members1 <- igraph::membership(comms_1)
  members2 <- igraph::membership(comms_2)
  t1 <- paste(names(freq_df)[1],"Community",sep="_")
  t2 <- paste(names(freq_df)[2],"Community",sep="_")
  freq_df[,t1] <- unname(members1[pull(freq_df[all_of(names(freq_df[1]))])])
  freq_df[,t2] <- unname(members2[pull(freq_df[all_of(names(freq_df[2]))])])
  gr <- c(all_of(t1),all_of(t2))
  dots <- lapply(gr, as.symbol)
  freq_df <- freq_df %>% dplyr::group_by(.dots=dots) %>% dplyr::summarise(weight = sum(weight))
  freq_df[,t1] = paste(t1,pull(freq_df[,t1]))
  freq_df[,t2] = paste(t2,pull(freq_df[,t2]))
  G <- igraph::graph_from_data_frame(freq_df, directed=FALSE)
  igraph::V(G)$type <- igraph::bipartite_mapping(G)$type
  output <- list(df = freq_df, graph = G)
}
