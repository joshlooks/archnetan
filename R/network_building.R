create_graph_single_file <- function(input_string, col1, col2) {
  '%>%' <- magrittr::'%>%'
  all_of <- tidyselect::all_of
  df_graph <- readr::read_csv(input_string)
  df_graph <- df_graph %>%
    tidyr::unite("Tup", all_of(col1):all_of(col2)) %>%
    dplyr::count(Tup) %>%
    tidyr::separate(Tup, sep = "_", into = c(col1,col2)) %>%
    dplyr::rename(weight = n)
  G <- igraph::graph_from_data_frame(df_graph, directed=FALSE)
  igraph::V(G)$type <- igraph::bipartite_mapping(G)$type
  output <- list(df = df_graph, graph = G)
}

create_graphs_two_files <- function(input_string1, input_string2, col1, col2, col3, col4, join_string){
  '%>%' <- magrittr::'%>%'
  all_of <- tidyselect::all_of
  as_tibble <- tibble::as_tibble
  df_f <- readr::read_csv(input_string1)
  df_f2 <- readr::read_csv(input_string2)
  df_full <- merge(df_f, df_f2, by=join_string, all=TRUE)

  df_graph_12 <- df_full[,c(col1, col2)] %>%
    tidyr::unite("Tup", all_of(col1):all_of(col2)) %>%
    dplyr::count(Tup) %>%
    tidyr::separate(Tup, sep = "_", into = c(col1, col2)) %>%
    dplyr::rename(weight = n)
  G_12 <- igraph::graph_from_data_frame(df_graph_12, directed=FALSE)
  igraph::V(G_12)$type <- igraph::bipartite_mapping(G_12)$type

  df_graph_34 <- df_full[,c(col3, col4)] %>%
    tidyr::unite("Tup", all_of(col3):all_of(col4)) %>%
    dplyr::count(Tup) %>%
    tidyr::separate(Tup, sep = "_", into = c(col3, col4)) %>%
    dplyr::rename(weight = n)
  G_34 <- igraph::graph_from_data_frame(df_graph_34, directed=FALSE)
  igraph::V(G_34)$type <- igraph::bipartite_mapping(G_34)$type

  output <- list(df_12 = as_tibble(df_graph_12), graph_12 = G_12,
                 df_34 = as_tibble(df_graph_34), graph_34 = G_34,
                 df_full = as_tibble(df_full))
}
