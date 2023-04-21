#' Creates a bipartite network from a single input file
#'
#' @param input_string A path to the input file (assumed to be a csv)
#' @param col1 A string of the first column to create a network from
#' @param col2 A string of the second column to create a network from
#' @returns An output list containing two components:
#' \item{df}{Tibble containing the edge weights of the bipartite network}
#' \item{graph}{igraph::graph object, the bipartite network}
#' @export
#' @examples
#' create_graph_single_file(archnetan_example("originalDataset.csv"),"Study Site","Source")
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

#' Creates a bipartite network from a single input file with weighted frequencies
#'
#' @param input_string A path to the input file (assumed to be a csv)
#' @param col1 A string of the first column to create a network from
#' @param col2 A string of the second column to create a network from
#' @param colWt A string of the column to weight the frequencies by
#' @returns An output list containing two components:
#' \item{df}{Tibble containing the edge weights of the bipartite network}
#' \item{graph}{igraph::graph object, the bipartite network}
#' @export
#' @examples
#' create_weighted_graph_single_file(archnetan_example("originalDataset.csv"),
#' "Study Site","Source","Weight (g)")
create_weighted_graph_single_file <- function(input_string, col1, col2, colWt) {
  '%>%' <- magrittr::'%>%'
  df_full <- readr::read_csv(input_string)
  df_temp <- dplyr::rename(df_full, wt = {{colWt}})
  df_graph <- df_temp %>%
    tidyr::unite("Tup", {{col1}}:{{col2}}) %>%
    dplyr::group_by(Tup) %>%
    dplyr::summarise(weight = sum(wt,na.rm=TRUE)) %>%
    tidyr::separate(Tup, sep = "_", into = c(col1,col2)) %>%
    dplyr::filter(weight > 0)
  G <- igraph::graph_from_data_frame(df_graph, directed=FALSE)
  igraph::V(G)$type <- igraph::bipartite_mapping(G)$type
  output <- list(df = df_graph, graph = G, df_full = df_full)
}

#' Creates one or two bipartite networks from two input files
#'
#' @param input_string1 A path to the first input file (assumed to be a csv)
#' @param input_string2 A path to the second input file (assumed to be a csv)
#' @param col1 First column of first network
#' @param col2 Second column of first network
#' @param col3 First column of second network (Optional)
#' @param col4 Second column of second network (Optional)
#' @param join_string Column to merge the two input csvs on
#' @returns An output list containing five components:
#' \item{df_12}{Tibble containing the edge weights of the first bipartite network}
#' \item{graph_12}{igraph::graph object, the first bipartite network}
#' \item{df_34}{Tibble containing the edge weights of the second bipartite network - requires col3,col4}
#' \item{graph_34}{igraph::graph object, the second bipartite network - requires col3,col4}
#' \item{df_full}{Raw, merged dataframe of the two input files}
#' @export
#' @examples
#' \dontrun{create_graphs_two_files("trial.csv", "trial2.csv", "Source",
#' "Litho", "Meta_Assemblage", "Meta_Litho", "Source")}
create_graphs_two_files <- function(input_string1, input_string2, col1, col2, col3=NULL, col4=NULL, join_string){
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

  if(is.null(col3)){
    output <- list(df = as_tibble(df_graph_12), graph = G_12, df_full = as_tibble(df_full))
  } else{
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
}

#' Creates a bipartite network from an inputted dataframe
#'
#' @param df A dataframe to create the network from (assumed to be a tibble)
#' @param col1 First column of network
#' @param col2 Second column network
#' @returns An output list containing two components:
#' \item{df}{Tibble containing the edge weights of the bipartite network}
#' \item{graph}{igraph::graph object, the bipartite network}
#' @export
#' @examples
#' #create_graph_df(df_full, "Source", "Combined")
create_graph_df <- function(df, col1, col2) {
  '%>%' <- magrittr::'%>%'
  all_of <- tidyselect::all_of
  df_graph <- df[,c(col1, col2)] %>%
    tidyr::unite("Tup", all_of(col1):all_of(col2)) %>%
    dplyr::count(Tup) %>%
    tidyr::separate(Tup, sep = "_", into = c(col1, col2)) %>%
    dplyr::rename(weight = n)
  G <- igraph::graph_from_data_frame(df_graph, directed=FALSE)
  igraph::V(G)$type <- igraph::bipartite_mapping(G)$type
  output <- list(df = df_graph, graph = G)
}

#' Creates a bipartite network from an inputted dataframe with weighted frequencies
#'
#' @param df A dataframe to create the network from (assumed to be a tibble)
#' @param col1 First column of network
#' @param col2 Second column network
#' @param colWt Column of weights in df
#' @returns An output list containing two components:
#' \item{df}{Tibble containing the edge weights of the bipartite network}
#' \item{graph}{igraph::graph object, the bipartite network}
#' @export
#' @examples
#' #create_weighted_graph_df(df_full, "Source", "Combined", "Weight (g)")
create_weighted_graph_df <- function(df, col1, col2, colWt) {
  '%>%' <- magrittr::'%>%'
  df_temp <- dplyr::rename(df, wt = {{colWt}})
  df_graph <- df_temp %>%
    tidyr::unite("Tup", {{col1}}:{{col2}}) %>%
    dplyr::group_by(Tup) %>%
    dplyr::summarise(weight = sum(wt,na.rm=TRUE)) %>%
    tidyr::separate(Tup, sep = "_", into = c(col1,col2)) %>%
    dplyr::filter(weight > 0)
  G <- igraph::graph_from_data_frame(df_graph, directed=FALSE)
  igraph::V(G)$type <- igraph::bipartite_mapping(G)$type
  output <- list(df = df_graph, graph = G, df_full = df_full)
}
