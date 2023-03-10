#' Bootstraps communities from unweighted dataset
#'
#' @param df A tibble of the unweighted dataset
#' @param col1 First column of network
#' @param col2 Second column of network
#' @param nBoot Number of bootstrap replications to create
#' @returns An output list containing three components:
#' \item{comm1}{Tibble containing bootstrapped communities of col1}
#' \item{comm2}{Tibble containing bootstrapped communities of col2}
#' \item{df}{Nested list of bootstrap replicant weighted dataframes}
#' @export
#' @examples
#' #create_boostrap_communities(base_df, "Source", "Combined")
create_bootstrap_communities <- function(df, col1, col2, nBoot=100){
  bootRepComm1 <- list()
  bootRepComm2 <- list()
  bootRepdf <- list()
  as_tibble <- tibble::as_tibble
  left_join <- dplyr::left_join
  tibble <- tibble::tibble
  bootN <- dim(df)[1]
  for (i in 1:nBoot){
    bootInd <- sample(1:bootN,bootN,replace=TRUE)
    thisBoot <- df[bootInd,]
    dfBoot <- create_graph_df(thisBoot, col1, col2)
    dfBoot <- as_tibble(dfBoot$df)
    bootRepdf <- append(bootRepdf, dfBoot)
    bootRCA <- calc_rca_sim(dfBoot, col1, col2)
    col1_rca_boot = calc_com_louv(bootRCA$RCA_df_col)
    col2_rca_boot = calc_com_louv(bootRCA$RCA_df_row)
    col1mem <- membership(col1_rca_boot$community)
    col2mem <- membership(col2_rca_boot$community)
    if (i==1){
      bootRepComm1 <- tibble(node_id = names(col1mem), community=col1mem)
      bootRepComm2 <- tibble(node_id = names(col2mem), community=col2mem)
    } else{
      bootTemp <- tibble(node_id = names(col1mem), community=col1mem)
      bootTemp2 <- tibble(node_id = names(col2mem), community=col2mem)
      bootRepComm1 <- left_join(bootRepComm1, bootTemp, by="node_id")
      bootRepComm2 <- left_join(bootRepComm2, bootTemp2, by="node_id")
    }
  }
  names(bootRepComm1) <- c(col1,paste0("Rep",1:nBoot))
  names(bootRepComm2) <- c(col2,paste0("Rep",1:nBoot))
  output <- list(comm1 = bootRepComm1, comm2 = bootRepComm2, df = bootRepdf)
}

#' Creates communities by switching edges within a weighted edgelist
#'
#' @param df A tibble of the weighted edgelist
#' @param col1 First column of network
#' @param col2 Second column of network
#' @returns An output list containing three components:
#' \item{comm1}{Tibble containing permuted communities of col1}
#' \item{comm2}{Tibble containing permuted communities of col2}
#' \item{df}{Nested list of permuted weighted dataframes}
#' @export
#' @examples
#' #create_permutation_communities(base_df, "Source", "Combined")
create_permutation_communities <- function(df, col1, col2){
  as_tibble <- tibble::as_tibble
  tibble <- tibble::tibble
  df <- as_tibble(df)
  perms <- t(combn(dim(df)[1],2))
  permRepdf <- list()
  left_join <- dplyr::left_join
  RCA <- calc_rca_sim(df, col1, col2)
  col1_rca_perm = calc_com_louv(RCA$RCA_df_col)
  col2_rca_perm = calc_com_louv(RCA$RCA_df_row)
  col1mem <- membership(col1_rca_perm$community)
  col2mem <- membership(col2_rca_perm$community)
  permRepComm1 <- tibble(node_id = names(col1mem), community=col1mem)
  permRepComm2 <- tibble(node_id = names(col2mem), community=col2mem)
  for (i in 1:dim(perms)[1]){
    temp <- df[perms[i,1],3]
    temp2 <- df[perms[i,2],3]
    thisPerm <- df
    thisPerm[perms[i,1],3] <- temp2
    thisPerm[perms[i,2],3] <- temp
    permRepdf <- append(permRepdf, thisPerm)
    permRCA <- calc_rca_sim(thisPerm, col1, col2)
    col1_rca_perm = calc_com_louv(permRCA$RCA_df_col)
    col2_rca_perm = calc_com_louv(permRCA$RCA_df_row)
    col1mem <- membership(col1_rca_perm$community)
    col2mem <- membership(col2_rca_perm$community)
    permTemp <- tibble(node_id = names(col1mem), community=col1mem)
    permTemp2 <- tibble(node_id = names(col2mem), community=col2mem)
    permRepComm1 <- left_join(permRepComm1, permTemp, by="node_id")
    permRepComm2 <- left_join(permRepComm2, permTemp2, by="node_id")
  }
  names(permRepComm1) <- c(col1,paste0("Rep",1:(dim(permRepComm1)[2]-1)))
  names(permRepComm2) <- c(col2,paste0("Rep",1:(dim(permRepComm2)[2]-1)))
  output <- list(comm1 = permRepComm1, comm2 = permRepComm2, df = permRepdf)
}
