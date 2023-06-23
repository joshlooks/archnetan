#' Calculates the cosine similarities of a bipartite network
#'
#' @param df_graph A dataframe of edge weights of a bipartite network
#' @returns An output list containing two components:
#' \item{col1_cossim}{Tibble containing the cosine similarities of first node type}
#' \item{col2_cossim}{Tibble containing the cosine similarities of the second node type}
#' @export
#' @examples
#' #calc_cos_sim(df_graph_main)
calc_cos_sim <- function(df_graph) {
  '%>%' <- magrittr::'%>%'
  col1 <- colnames(df_graph)[1]
  col2 <- colnames(df_graph)[2]
  pull <- dplyr::pull
  df_c1vec <- df_graph %>% tidyr::spread(key = col2, value = weight, fill = 0)
  df_c2vec <- df_graph %>% tidyr::spread(key = col1, value = weight, fill = 0)
  temp <- combinat::combn(unique(pull(df_c1vec[col1])), 2)
  temp2 <- combinat::combn(unique(pull(df_c2vec[col2])), 2)
  df_c1edg <- tidyr::tibble(temp[1,],temp[2,]) %>%
    dplyr::rename(p1 = 1, p2 = 2) %>%
    tibble::add_column(sim = 0)
  df_c2edg <- tidyr::tibble(temp2[1,],temp2[2,]) %>%
    dplyr::rename(p1 = 1, p2 = 2) %>%
    tibble::add_column(sim = 0)
  nc1 = ncol(df_c1vec)
  nc2 = ncol(df_c2vec)
  for(i in 1:nrow(df_c1edg)){
    temp <- as.numeric(df_c1vec[pull(df_c1vec[col1]) == df_c1edg$p1[i], 2:nc1])
    temp2 <- as.numeric(df_c1vec[pull(df_c1vec[col1]) == df_c1edg$p2[i], 2:nc1])
    df_c1edg$sim[i]=(temp %*% temp2) / (norm(temp,type = "2")*norm(temp2, type = "2"))
  }
  for(i in 1:nrow(df_c2edg)){
    temp <- as.numeric(df_c2vec[pull(df_c2vec[col2]) == df_c2edg$p1[i], 2:nc2])
    temp2 <- as.numeric(df_c2vec[pull(df_c2vec[col2]) == df_c2edg$p2[i], 2:nc2])
    df_c2edg$sim[i] <- (temp %*% temp2) / (norm(temp,type = "2")*norm(temp2, type = "2"))
  }
  output <- list(col1_cossim = df_c1edg, col2_cossim = df_c2edg)
}

#' Calculates the cosine similarities of a bipartite network using RCA normalisation
#'
#' @param df_graph A dataframe of edge weights of a bipartite network
#' @returns An output list containing three components:
#' \item{col1_cossim}{Tibble containing the cosine similarities of first node type/column after RCA normalisation}
#' \item{col2_cossim}{Tibble containing the cosine similarities of the second node type/column after RCA normalisation}
#' \item{RCA_full}{Tibble containing the RCA values of both columns/node types}
#' @export
#' @examples
#' #calc_rca_sim(df_graph_main)
calc_rca_sim <- function(df_graph){
  calc_rca <- function(m, cols, rows) {
    sum_cols <- colSums(m)
    sum_rows <- rowSums(m)
    Rn <- m / matrix(rep(sum_rows, cols), ncol = cols)
    Rd <- t(matrix(rep(sum_cols, rows), ncol = rows)) / matrix(sum(m), rows, cols)
    R <- Rn / Rd
  }

  calc_phi <- function(m) {
    m[which(m >= 1)] <- 1
    m[which(m < 1)] <- 0
    sum_cols <- colSums(m)
    mm <- t(m) %*% m
    phi <- mm / outer(sum_cols,sum_cols,pmax)
    diag(phi) <- 0
    phi
  }

  create_phi_df <- function(phi_mat, c1, cnames) {
    '%>%' <- magrittr::'%>%'
    cols <- ncol(phi_mat)
    df_temp <- matrix(, nrow=cols, ncol=(cols+1))
    colnames(df_temp) <- c(c1,cnames)
    df_final <- tibble::as_tibble(df_temp)
    df_final[,2:(cols+1)] <- phi_mat
    df_final[,1] <- cnames
    df_final
  }
  type1 <- colnames(df_graph)[1]
  type2 <- colnames(df_graph)[2]
  '%>%' <- magrittr::'%>%'
  as_tibble <- tibble::as_tibble
  df_wide <- df_graph %>% tidyr::spread(key = type1, value = weight, fill = 0)
  df_rca <- df_wide
  rows <- nrow(df_rca)
  cols <- ncol(df_rca)
  df_mat <- as.matrix(df_rca[,2:cols])
  dimnames(df_mat) <- NULL
  vals <- calc_rca(df_mat, (cols-1), rows)
  phi_1 <- calc_phi(vals)
  phi_2 <- calc_phi(t(vals))
  df_phi_1 <- create_phi_df(phi_1, type1, colnames(df_wide)[2:cols])
  df_phi_2 <- create_phi_df(phi_2, type2, dplyr::pull(df_wide[,1]))
  df_phi_1 <- tidyr::pivot_longer(df_phi_1,-1,names_to="Secondary",values_to="RCA")
  df_phi_1 <- df_phi_1 %>% dplyr::rowwise() %>%
    dplyr::mutate(temp = paste(sort(c(get(type1),Secondary)), collapse = " "))
  df_phi_1 <- df_phi_1[!duplicated(df_phi_1$temp),1:3]
  df_phi_1 <- df_phi_1[df_phi_1$RCA > 0,]
  df_phi_2 <- tidyr::pivot_longer(df_phi_2,-1,names_to="Secondary",values_to="RCA")
  df_phi_2 <- df_phi_2 %>% dplyr::rowwise() %>%
    dplyr::mutate(temp = paste(sort(c(get(type2),Secondary)), collapse = " "))
  df_phi_2 <- df_phi_2[!duplicated(df_phi_2$temp),1:3]
  df_phi_2 <- df_phi_2[df_phi_2$RCA > 0,]
  df_rca <- as_tibble(vals)
  colnames(df_rca) <- colnames(df_wide)[2:cols]
  df_rca[,paste(type2)] <- df_wide[,type2]
  df_rca <- df_rca[,c(colnames(df_rca)[cols],colnames(df_rca)[1:(cols-1)])]
  output <- list(col1_cossim = df_phi_1, col2_cossim = df_phi_2, RCA_full = df_rca)
}

#' Calculates the similarities of a bipartite network using the Brainerd-Robinson coefficient
#'
#' @param df_graph A dataframe of edge weights of a bipartite network
#' @returns An output list containing three components:
#' \item{col1_brsim}{Tibble containing the Brainerd-Robinson coefficients of first node type/columns}
#' \item{col2_brsim}{Tibble containing the Brainerd-Robinson coefficients of second node type/columns}
#' \item{col1_brmat}{Tibble containing the Brainerd-Robinson coefficient matrix for the first column}
#' \item{col2_brmat}{Tibble containing the Brainerd-Robinson coefficient matrix for the second column}
#' @export
#' @examples
#' #calc_br_sim(df_graph_main)
calc_br_sim <- function(df_graph){
  brainerd_robinson <- function(df) {
    ncols <- dim(df)[1]
    br_coeffs <- matrix(0,ncols,ncols)
    for (i in 1:ncols) {
      for (j in 1:ncols) {
        temp1 <- df[i,]
        temp2 <- df[j,]
        br_coeffs[i,j] <- 1 - (sum(abs(temp1 - temp2)))/2}}
    return(br_coeffs)
  }

  create_br_df <- function(br_mat, c1, cnames) {
    '%>%' <- magrittr::'%>%'
    ncols <- ncol(br_mat)
    df_temp <- matrix(, nrow=ncols, ncol=(ncols+1))
    colnames(df_temp) <- c(c1,cnames)
    df_final <- tibble::as_tibble(df_temp)
    df_final[,2:(ncols+1)] <- br_mat
    df_final[,1] <- cnames
    df_final
  }

  '%>%' <- magrittr::'%>%'
  as_tibble <- tibble::as_tibble
  col1 <- colnames(df_graph)[1]
  col2 <- colnames(df_graph)[2]
  df_wide <- df_graph %>% tidyr::spread(key = {{col1}}, value = weight, fill = 0)
  df_wide_2 <- df_graph %>% tidyr::spread(key = {{col2}}, value = weight, fill = 0)
  rows <- nrow(df_wide)
  cols <- ncol(df_wide)
  df_mat <- prop.table(as.matrix(df_wide[,2:cols]),1)
  df_mat_2 <- prop.table(as.matrix(df_wide_2[,2:rows]),1)
  vals <- brainerd_robinson(df_mat)
  vals2 <- brainerd_robinson(df_mat_2)
  df_br_col1 <- create_br_df(vals, col1, colnames(df_wide_2)[2:(rows+1)])
  df_br_col2 <- create_br_df(vals2, col2, colnames(df_wide)[2:(cols)])
  BR_df_col1 <- tidyr::pivot_longer(df_br_col1,-1,names_to="Secondary",values_to="BR")
  BR_df_col1 <- BR_df_col1 %>% dplyr::rowwise() %>%
    dplyr::mutate(temp = paste(sort(c(get(col1),Secondary)),collapse = " "))
  BR_df_col1 <- BR_df_col1[!duplicated(BR_df_col1$temp),1:3]
  BR_df_col1 <- BR_df_col1[BR_df_col1[,1] != BR_df_col1[,2],]
  BR_df_col2 <- tidyr::pivot_longer(df_br_col2,-1,names_to="Secondary",values_to="BR")
  BR_df_col2 <- BR_df_col2 %>% dplyr::rowwise() %>%
    dplyr::mutate(temp = paste(sort(c(get(col2),Secondary)),collapse = " "))
  BR_df_col2 <- BR_df_col2[!duplicated(BR_df_col2$temp),1:3]
  BR_df_col2 <- BR_df_col2[BR_df_col2[,1] != BR_df_col2[,2],]
  output <- list(col1_brsim = BR_df_col2, col2_brsim = BR_df_col1, col1_brmat = df_br_col2, col2_brmat = df_br_col1)
}
