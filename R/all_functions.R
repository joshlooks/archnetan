create_graph_single_file <- function(input_string, col1, col2) {
  '%>%' <- magrittr::'%>%'
  df_full <- readr::read_csv(input_string)
  df_graph <- df_full %>%
    tidyr::unite("Tup", {{col1}}:{{col2}}) %>%
    dplyr::count(Tup) %>%
    tidyr::separate(Tup, sep = "_", into = c(col1,col2)) %>%
    dplyr::rename(weight = n) %>%
    dplyr::filter(weight > 0)
  G <- igraph::graph_from_data_frame(df_graph, directed=FALSE)
  igraph::V(G)$type <- igraph::bipartite_mapping(G)$type
  output <- list(df = df_graph, graph = G, df_full = df_full)
}

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

create_graph_df <- function(df, col1, col2) {
  '%>%' <- magrittr::'%>%'
  df_graph <- df %>%
    dplyr::select(c({{col1}},{{col2}})) %>%
    tidyr::unite("Tup", {{col1}}:{{col2}}) %>%
    dplyr::count(Tup) %>%
    tidyr::separate(Tup, sep = "_", into = c(col1, col2)) %>%
    dplyr::rename(weight = n) %>%
    dplyr::filter(weight > 0)
  G <- igraph::graph_from_data_frame(df_graph, directed=FALSE)
  igraph::V(G)$type <- igraph::bipartite_mapping(G)$type
  output <- list(df = df_graph, graph = G)
}

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

create_graphs_two_files <- function(input_string1, input_string2, col1, col2, col3=NULL, col4=NULL, join_string){
  '%>%' <- magrittr::'%>%'
  as_tibble <- tibble::as_tibble
  df_f <- readr::read_csv(input_string1)
  df_f2 <- readr::read_csv(input_string2)
  df_full <- merge(df_f, df_f2, by=join_string, all=TRUE)

  df_graph_12 <- df_full %>%
    dplyr::select(c({{col1}},{{col2}})) %>%
    tidyr::unite("Tup", {{col1}}:{{col2}}) %>%
    dplyr::count(Tup) %>%
    tidyr::separate(Tup, sep = "_", into = c(col1, col2)) %>%
    dplyr::rename(weight = n) %>%
    dplyr::filter(weight > 0)
  G_12 <- igraph::graph_from_data_frame(df_graph_12, directed=FALSE)
  igraph::V(G_12)$type <- igraph::bipartite_mapping(G_12)$type

  if(is.null(col3)){
    output <- list(df = as_tibble(df_graph_12), graph = G_12, df_full = as_tibble(df_full))
  } else{
    df_graph_34 <- df_full %>%
      dplyr::select(c({{col3}},{{col4}})) %>%
      tidyr::unite("Tup", {{col3}}:{{col4}}) %>%
      dplyr::count(Tup) %>%
      tidyr::separate(Tup, sep = "_", into = c(col3, col4)) %>%
      dplyr::rename(weight = n) %>%
      dplyr::filter(weight > 0)
    G_34 <- igraph::graph_from_data_frame(df_graph_34, directed=FALSE)
    igraph::V(G_34)$type <- igraph::bipartite_mapping(G_34)$type

    output <- list(df_12 = as_tibble(df_graph_12), graph_12 = G_12,
                   df_34 = as_tibble(df_graph_34), graph_34 = G_34,
                   df_full = as_tibble(df_full))
  }
}

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

calc_rca_sim <- function(df_graph, type1="Source", type2="POI"){
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
  '%>%' <- magrittr::'%>%'
  as_tibble <- tibble::as_tibble
  df_wide <- df_graph %>% tidyr::spread(key = type1, value = weight, fill = 0)
  df_rca <- df_wide
  rows <- nrow(df_rca)
  cols <- ncol(df_rca)
  df_mat <- as.matrix(df_rca[,2:cols])
  dimnames(df_mat) <- NULL
  vals <- calc_rca(df_mat, (cols-1), rows)
  phi_source <- calc_phi(vals)
  phi_place <- calc_phi(t(vals))
  df_phi_source <- create_phi_df(phi_source, type1, colnames(df_wide)[2:cols])
  df_phi_place <- create_phi_df(phi_place, type2, dplyr::pull(df_wide[,1]))
  df_phi_source <- tidyr::pivot_longer(df_phi_source,-1,names_to="Secondary",values_to="RCA")
  df_phi_source <- df_phi_source %>% dplyr::rowwise() %>%
    dplyr::mutate(temp = paste(sort(c(get(type1),Secondary)), collapse = " "))
  df_phi_source <- df_phi_source[!duplicated(df_phi_source$temp),1:3]
  df_phi_source <- df_phi_source[df_phi_source$RCA > 0,]
  df_phi_place <- tidyr::pivot_longer(df_phi_place,-1,names_to="Secondary",values_to="RCA")
  df_phi_place <- df_phi_place %>% dplyr::rowwise() %>%
    dplyr::mutate(temp = paste(sort(c(get(type2),Secondary)), collapse = " "))
  df_phi_place <- df_phi_place[!duplicated(df_phi_place$temp),1:3]
  df_phi_place <- df_phi_place[df_phi_place$RCA > 0,]
  df_rca <- as_tibble(vals)
  colnames(df_rca) <- colnames(df_wide)[2:cols]
  df_rca[,paste(type2)] <- df_wide[,type2]
  df_rca <- df_rca[,c(colnames(df_rca)[cols],colnames(df_rca)[1:(cols-1)])]
  output <- list(RCA_df_row = df_phi_place, RCA_df_col = df_phi_source, RCA_full = df_rca)
}

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
  output <- list(BR_df_col1 = BR_df_col2, BR_df_col2 = BR_df_col1, BR_mat_col1 = df_br_col2, BR_mat_col2 = df_br_col1)
}

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

RCA_communities_single <- function(df, col1, col2, proportion=0.10, numReps=100){
  as_tibble <- tibble::as_tibble
  tibble <- tibble::tibble
  df <- as_tibble(df)
  nRows <- dim(df)[1]
  numNoise <- ceiling(nRows*proportion)
  left_join <- dplyr::left_join
  noiseRepdf <- list()
  origRCA <- calc_rca_sim(df, col1, col2)
  origMem1 <- membership(calc_com_louv(origRCA$RCA_df_col)$community)
  origMem2 <- membership(calc_com_louv(origRCA$RCA_df_row)$community)
  col1Diff <- rep(0,numReps)
  col2Diff <- rep(0,numReps)
  col1Change <- matrix(0, length(origMem1), (numReps+1))
  col1Change <- as_tibble(col1Change)
  col1Change[,1] <- names(origMem1)
  names(col1Change) <- c(col1,paste0("Rep",1:numReps))
  col2Change <- matrix(0, length(origMem2), (numReps+1))
  col2Change <- as_tibble(col2Change)
  col2Change[,1] <- names(origMem2)
  names(col2Change) <- c(col2,paste0("Rep",1:numReps))
  wts <- dplyr::pull(df[,3])
  for (i in 1:numReps){
    dfRep <- df
    noiseInd <- sample(1:nRows,numNoise,replace=TRUE,prob=wts)
    for (j in noiseInd){
      dfRep[j,3] <- dfRep[j,3]+1
    }
    noiseRepdf <- append(noiseRepdf, dfRep)
    noiseRCA <- calc_rca_sim(dfRep, col1, col2)
    col1_rca_noise = calc_com_louv(noiseRCA$RCA_df_col)
    col2_rca_noise = calc_com_louv(noiseRCA$RCA_df_row)
    col1mem <- membership(col1_rca_noise$community)
    col2mem <- membership(col2_rca_noise$community)
    repTemp1 <- tibble(node_id = names(col1mem), community=col1mem)
    repTemp2 <- tibble(node_id = names(col2mem), community=col2mem)
    swaps <- combinat::permn(max(col1mem))
    swaps2 <- combinat::permn(max(col2mem))
    dif <- rep(0,length(swaps))
    dif2 <- rep(0,length(swaps2))
    for (k in 1:length(swaps)){
      temp <- swaps[[k]]
      col1tempMem <- temp[col1mem]
      dif[k] <- sum(origMem1 != col1tempMem)
    }
    for (k in 1:length(swaps2)){
      temp <- swaps2[[k]]
      col2tempMem <- temp[col2mem]
      dif2[k] <- sum(origMem2 != col2tempMem)
    }
    col1Diff[i] <- min(dif)
    col2Diff[i] <- min(dif2)
    col1Change[,i+1] <- as.integer(swaps[[which.min(dif)]][col1mem] != origMem1)
    col2Change[,i+1] <- as.integer(swaps2[[which.min(dif2)]][col2mem] != origMem2)
    if (i==1){
      noiseRepComm1 <- tibble(node_id = names(col1mem), community=swaps[[which.min(dif)]][col1mem])
      noiseRepComm2 <- tibble(node_id = names(col2mem), community=swaps2[[which.min(dif2)]][col2mem])
    } else{
      repTemp <- tibble(node_id = names(col1mem), community=swaps[[which.min(dif)]][col1mem])
      repTemp2 <- tibble(node_id = names(col2mem), community=swaps2[[which.min(dif2)]][col2mem])
      noiseRepComm1 <- left_join(noiseRepComm1, repTemp, by="node_id")
      noiseRepComm2 <- left_join(noiseRepComm2, repTemp2, by="node_id")
    }
  }
  names(noiseRepComm1) <- c(col1,paste0("Rep",1:(dim(noiseRepComm1)[2]-1)))
  names(noiseRepComm2) <- c(col2,paste0("Rep",1:(dim(noiseRepComm2)[2]-1)))
  output <- list(comm1 = noiseRepComm1, comm2 = noiseRepComm2,
                 hamming1 = col1Diff, hamming2 = col2Diff,
                 col1Changes = col1Change, col2Changes = col2Change,
                 df = noiseRepdf)
}

cos_communities_single <- function(df, col1, col2, proportion=0.10, numReps=100){
  as_tibble <- tibble::as_tibble
  tibble <- tibble::tibble
  df <- as_tibble(df)
  nRows <- dim(df)[1]
  numNoise <- ceiling(nRows*proportion)
  left_join <- dplyr::left_join
  noiseRepdf <- list()
  origCos <- calc_cos_sim(df)
  origMem1 <- membership(calc_com_louv(origCos$col1_cossim)$community)
  origMem2 <- membership(calc_com_louv(origCos$col2_cossim)$community)
  col1Diff <- rep(0,numReps)
  col2Diff <- rep(0,numReps)
  col1Change <- matrix(0, length(origMem1), (numReps+1))
  col1Change <- as_tibble(col1Change)
  col1Change[,1] <- names(origMem1)
  names(col1Change) <- c(col1,paste0("Rep",1:numReps))
  col2Change <- matrix(0, length(origMem2), (numReps+1))
  col2Change <- as_tibble(col2Change)
  col2Change[,1] <- names(origMem2)
  names(col2Change) <- c(col2,paste0("Rep",1:numReps))
  wts <- dplyr::pull(df[,3])
  for (i in 1:numReps){
    dfRep <- df
    noiseInd <- sample(1:nRows,numNoise,replace=TRUE,prob=wts)
    for (j in noiseInd){
      dfRep[j,3] <- dfRep[j,3]+1
    }
    noiseRepdf <- append(noiseRepdf, dfRep)
    noiseCos <- calc_cos_sim(dfRep)
    col1_cos_noise = calc_com_louv(noiseCos$col1_cossim)
    col2_cos_noise = calc_com_louv(noiseCos$col2_cossim)
    col1mem <- membership(col1_cos_noise$community)
    col2mem <- membership(col2_cos_noise$community)
    repTemp1 <- tibble(node_id = names(col1mem), community=col1mem)
    repTemp2 <- tibble(node_id = names(col2mem), community=col2mem)
    swaps <- combinat::permn(max(col1mem))
    swaps2 <- combinat::permn(max(col2mem))
    dif <- rep(0,length(swaps))
    dif2 <- rep(0,length(swaps2))
    for (k in 1:length(swaps)){
      temp <- swaps[[k]]
      col1tempMem <- temp[col1mem]
      dif[k] <- sum(origMem1 != col1tempMem)
    }
    for (k in 1:length(swaps2)){
      temp <- swaps2[[k]]
      col2tempMem <- temp[col2mem]
      dif2[k] <- sum(origMem2 != col2tempMem)
    }
    col1Diff[i] <- min(dif)
    col2Diff[i] <- min(dif2)
    col1Change[,i+1] <- as.integer(swaps[[which.min(dif)]][col1mem] != origMem1)
    col2Change[,i+1] <- as.integer(swaps2[[which.min(dif2)]][col2mem] != origMem2)
    if (i==1){
      noiseRepComm1 <- tibble(node_id = names(col1mem), community=swaps[[which.min(dif)]][col1mem])
      noiseRepComm2 <- tibble(node_id = names(col2mem), community=swaps2[[which.min(dif2)]][col2mem])
    } else{
      repTemp <- tibble(node_id = names(col1mem), community=swaps[[which.min(dif)]][col1mem])
      repTemp2 <- tibble(node_id = names(col2mem), community=swaps2[[which.min(dif2)]][col2mem])
      noiseRepComm1 <- left_join(noiseRepComm1, repTemp, by="node_id")
      noiseRepComm2 <- left_join(noiseRepComm2, repTemp2, by="node_id")
    }
  }
  names(noiseRepComm1) <- c(col1,paste0("Rep",1:(dim(noiseRepComm1)[2]-1)))
  names(noiseRepComm2) <- c(col2,paste0("Rep",1:(dim(noiseRepComm2)[2]-1)))
  output <- list(comm1 = noiseRepComm1, comm2 = noiseRepComm2,
                 hamming1 = col1Diff, hamming2 = col2Diff,
                 col1Changes = col1Change, col2Changes = col2Change,
                 df = noiseRepdf)
}

cos_communities_multi <- function(df, col1, col2, proportion=0.10, numReps=100){
  as_tibble <- tibble::as_tibble
  tibble <- tibble::tibble
  df <- as_tibble(df)
  nRows <- dim(df)[1]
  numNoise <- ceiling(nRows*proportion)
  left_join <- dplyr::left_join
  noiseRepdf <- list()
  origCos <- calc_cos_sim(df)
  origMem1 <- membership(calc_com_louv(origCos$col1_cossim)$community)
  origMem2 <- membership(calc_com_louv(origCos$col2_cossim)$community)
  col1Diff <- rep(0,numReps)
  col2Diff <- rep(0,numReps)
  col1Change <- matrix(0, length(origMem1), (numReps+1))
  col1Change <- as_tibble(col1Change)
  col1Change[,1] <- names(origMem1)
  names(col1Change) <- c(col1,paste0("Rep",1:numReps))
  col2Change <- matrix(0, length(origMem2), (numReps+1))
  col2Change <- as_tibble(col2Change)
  col2Change[,1] <- names(origMem2)
  names(col2Change) <- c(col2,paste0("Rep",1:numReps))
  col3Name <- names(df)[3]
  wtsCol1 <- dplyr::group_by(df, !!! rlang::syms(col1)) %>%
    dplyr::summarise(weight = sum(!!! rlang::syms(col3Name)))
  wtsCol2 <- dplyr::group_by(df, !!! rlang::syms(col2)) %>%
    dplyr::summarise(weight = sum(!!! rlang::syms(col3Name)))
  for (i in 1:numReps){
    dfRep <- df
    col1Ind <- sample(dplyr::pull(wtsCol1[,1]),numNoise,replace=TRUE,prob=wtsCol1$weight)
    col2Ind <- sample(dplyr::pull(wtsCol2[,1]),numNoise,replace=TRUE,prob=wtsCol2$weight)
    for (j in 1:numNoise){
      c1 <- col1Ind[j]
      c2 <- col2Ind[j]
      if (!any((dfRep[,1]==c1)&(dfRep[,2]==c2))){
        dfRep <- dfRep %>% dplyr::add_row(!!!setNames(list(c1,c2,1), names(.)))
      } else{
        dfRep[((dfRep[,1]==c1)&(dfRep[,2]==c2)),3] <-
          dfRep[((dfRep[,1]==c1)&(dfRep[,2]==c2)),3] + 1
      }
    }
    noiseRepdf <- append(noiseRepdf, dfRep)
    noiseCos <- calc_cos_sim(dfRep)
    col1_cos_noise = calc_com_louv(noiseCos$col1_cossim)
    col2_cos_noise = calc_com_louv(noiseCos$col2_cossim)
    col1mem <- membership(col1_cos_noise$community)
    col2mem <- membership(col2_cos_noise$community)
    repTemp1 <- tibble(node_id = names(col1mem), community=col1mem)
    repTemp2 <- tibble(node_id = names(col2mem), community=col2mem)
    swaps <- combinat::permn(max(col1mem))
    swaps2 <- combinat::permn(max(col2mem))
    dif <- rep(0,length(swaps))
    dif2 <- rep(0,length(swaps2))
    for (k in 1:length(swaps)){
      temp <- swaps[[k]]
      col1tempMem <- temp[col1mem]
      dif[k] <- sum(origMem1 != col1tempMem)
    }
    for (k in 1:length(swaps2)){
      temp <- swaps2[[k]]
      col2tempMem <- temp[col2mem]
      dif2[k] <- sum(origMem2 != col2tempMem)
    }
    col1Diff[i] <- min(dif)
    col2Diff[i] <- min(dif2)
    col1Change[,i+1] <- as.integer(swaps[[which.min(dif)]][col1mem] != origMem1)
    col2Change[,i+1] <- as.integer(swaps2[[which.min(dif2)]][col2mem] != origMem2)
    if (i==1){
      noiseRepComm1 <- tibble(node_id = names(col1mem), community=swaps[[which.min(dif)]][col1mem])
      noiseRepComm2 <- tibble(node_id = names(col2mem), community=swaps2[[which.min(dif2)]][col2mem])
    } else{
      repTemp <- tibble(node_id = names(col1mem), community=swaps[[which.min(dif)]][col1mem])
      repTemp2 <- tibble(node_id = names(col2mem), community=swaps2[[which.min(dif2)]][col2mem])
      noiseRepComm1 <- left_join(noiseRepComm1, repTemp, by="node_id")
      noiseRepComm2 <- left_join(noiseRepComm2, repTemp2, by="node_id")
    }
  }
  names(noiseRepComm1) <- c(col1,paste0("Rep",1:(dim(noiseRepComm1)[2]-1)))
  names(noiseRepComm2) <- c(col2,paste0("Rep",1:(dim(noiseRepComm2)[2]-1)))
  output <- list(comm1 = noiseRepComm1, comm2 = noiseRepComm2,
                 hamming1 = col1Diff, hamming2 = col2Diff,
                 col1Changes = col1Change, col2Changes = col2Change,
                 df = noiseRepdf)
}

RCA_communities_multi <- function(df, col1, col2, proportion=0.10, numReps=100){
  as_tibble <- tibble::as_tibble
  tibble <- tibble::tibble
  df <- as_tibble(df)
  nRows <- dim(df)[1]
  numNoise <- ceiling(nRows*proportion)
  left_join <- dplyr::left_join
  noiseRepdf <- list()
  origRCA <- calc_rca_sim(df, col1, col2)
  origMem1 <- membership(calc_com_louv(origRCA$RCA_df_col)$community)
  origMem2 <- membership(calc_com_louv(origRCA$RCA_df_row)$community)
  col1Diff <- rep(0,numReps)
  col2Diff <- rep(0,numReps)
  col1Change <- matrix(0, length(origMem1), (numReps+1))
  col1Change <- as_tibble(col1Change)
  col1Change[,1] <- names(origMem1)
  names(col1Change) <- c(col1,paste0("Rep",1:numReps))
  col2Change <- matrix(0, length(origMem2), (numReps+1))
  col2Change <- as_tibble(col2Change)
  col2Change[,1] <- names(origMem2)
  names(col2Change) <- c(col2,paste0("Rep",1:numReps))
  col3Name <- names(df)[3]
  wtsCol1 <- dplyr::group_by(df, !!! rlang::syms(col1)) %>%
    dplyr::summarise(weight = sum(!!! rlang::syms(col3Name)))
  wtsCol2 <- dplyr::group_by(df, !!! rlang::syms(col2)) %>%
    dplyr::summarise(weight = sum(!!! rlang::syms(col3Name)))
  for (i in 1:numReps){
    dfRep <- df
    col1Ind <- sample(dplyr::pull(wtsCol1[,1]),numNoise,replace=TRUE,prob=wtsCol1$weight)
    col2Ind <- sample(dplyr::pull(wtsCol2[,1]),numNoise,replace=TRUE,prob=wtsCol2$weight)
    for (j in 1:numNoise){
      c1 <- col1Ind[j]
      c2 <- col2Ind[j]
      if (!any((dfRep[,1]==c1)&(dfRep[,2]==c2))){
        dfRep <- dfRep %>% dplyr::add_row(!!!setNames(list(c1,c2,1), names(.)))
      } else{
        dfRep[((dfRep[,1]==c1)&(dfRep[,2]==c2)),3] <-
          dfRep[((dfRep[,1]==c1)&(dfRep[,2]==c2)),3] + 1
      }
    }
    noiseRepdf <- append(noiseRepdf, dfRep)
    noiseRCA <- calc_rca_sim(dfRep, col1, col2)
    col1_rca_noise = calc_com_louv(noiseRCA$RCA_df_col)
    col2_rca_noise = calc_com_louv(noiseRCA$RCA_df_row)
    col1mem <- membership(col1_rca_noise$community)
    col2mem <- membership(col2_rca_noise$community)
    repTemp1 <- tibble(node_id = names(col1mem), community=col1mem)
    repTemp2 <- tibble(node_id = names(col2mem), community=col2mem)
    swaps <- combinat::permn(max(col1mem))
    swaps2 <- combinat::permn(max(col2mem))
    dif <- rep(0,length(swaps))
    dif2 <- rep(0,length(swaps2))
    for (k in 1:length(swaps)){
      temp <- swaps[[k]]
      col1tempMem <- temp[col1mem]
      dif[k] <- sum(origMem1 != col1tempMem)
    }
    for (k in 1:length(swaps2)){
      temp <- swaps2[[k]]
      col2tempMem <- temp[col2mem]
      dif2[k] <- sum(origMem2 != col2tempMem)
    }
    col1Diff[i] <- min(dif)
    col2Diff[i] <- min(dif2)
    col1Change[,i+1] <- as.integer(swaps[[which.min(dif)]][col1mem] != origMem1)
    col2Change[,i+1] <- as.integer(swaps2[[which.min(dif2)]][col2mem] != origMem2)
    if (i==1){
      noiseRepComm1 <- tibble(node_id = names(col1mem), community=swaps[[which.min(dif)]][col1mem])
      noiseRepComm2 <- tibble(node_id = names(col2mem), community=swaps2[[which.min(dif2)]][col2mem])
    } else{
      repTemp <- tibble(node_id = names(col1mem), community=swaps[[which.min(dif)]][col1mem])
      repTemp2 <- tibble(node_id = names(col2mem), community=swaps2[[which.min(dif2)]][col2mem])
      noiseRepComm1 <- left_join(noiseRepComm1, repTemp, by="node_id")
      noiseRepComm2 <- left_join(noiseRepComm2, repTemp2, by="node_id")
    }
  }
  names(noiseRepComm1) <- c(col1,paste0("Rep",1:(dim(noiseRepComm1)[2]-1)))
  names(noiseRepComm2) <- c(col2,paste0("Rep",1:(dim(noiseRepComm2)[2]-1)))
  output <- list(comm1 = noiseRepComm1, comm2 = noiseRepComm2,
                 hamming1 = col1Diff, hamming2 = col2Diff,
                 col1Changes = col1Change, col2Changes = col2Change,
                 df = noiseRepdf)
}

create_added_noise_RCA_communities <- function(df, col1, col2, proportion=0.10, numReps=100, method='single', propSave=NULL){
  if (method=='single'){
    RCAc <- RCA_communities_single
  } else {
    RCAc <- RCA_communities_multi
  }
  if (length(proportion) == 1){
    output <- RCAc(df, col1, col2, proportion, numReps)
  } else{
    as_tibble <- tibble::as_tibble
    tibble <- tibble::tibble
    col1Stats <- matrix(0, length(proportion), 3)
    col1Stats <- as_tibble(col1Stats)
    col1Stats[,1] <- proportion
    names(col1Stats) <- c("Proportion","Mean","Std")
    col2Stats <- matrix(0, length(proportion), 3)
    col2Stats <- as_tibble(col2Stats)
    col2Stats[,1] <- proportion
    names(col2Stats) <- c("Proportion","Mean","Std")
    outSave <- NULL
    for (l in 1:length(proportion)){
      prop <- proportion[l]
      outProp <- RCAc(df, col1, col2, proportion=prop, numReps)
      if (prop == propSave){
        outSave <- outProp
      }
      col1Stats[l,2] <- mean(outProp$hamming1)
      col1Stats[l,3] <- sd(outProp$hamming1)
      col2Stats[l,2] <- mean(outProp$hamming2)
      col2Stats[l,3] <- sd(outProp$hamming2)
    }
    output <- list(col1Stats = col1Stats, col2Stats = col2Stats, propSave <- outSave)
  }
}

create_added_noise_cos_communities <- function(df, col1, col2, proportion=0.10, numReps=100, method='single', propSave=NULL){
  if (method=='single'){
    cc <- cos_communities_single
  } else{
    cc <- cos_communities_multi
  }
  if (length(proportion) == 1){
    output <- cc(df, col1, col2, proportion, numReps)
  } else{
    as_tibble <- tibble::as_tibble
    tibble <- tibble::tibble
    col1Stats <- matrix(0, length(proportion), 3)
    col1Stats <- as_tibble(col1Stats)
    col1Stats[,1] <- proportion
    names(col1Stats) <- c("Proportion","Mean","Std")
    col2Stats <- matrix(0, length(proportion), 3)
    col2Stats <- as_tibble(col2Stats)
    col2Stats[,1] <- proportion
    names(col2Stats) <- c("Proportion","Mean","Std")
    outSave <- NULL
    for (l in 1:length(proportion)){
      prop <- proportion[l]
      outProp <- cc(df, col1, col2, proportion=prop, numReps)
      if (prop == propSave){
        outSave <- outProp
      }
      col1Stats[l,2] <- mean(outProp$hamming1)
      col1Stats[l,3] <- sd(outProp$hamming1)
      col2Stats[l,2] <- mean(outProp$hamming2)
      col2Stats[l,3] <- sd(outProp$hamming2)
    }
    output <- list(col1Stats = col1Stats, col2Stats = col2Stats, propSave <- outSave)
  }
}

create_bootstrap_cos_communities <- function(df, df_full, col1, col2, nBoot=100, method='single'){
  # Function to deal with missing nodes during resampling
  copy_and_add_list <- function(original, replication){
    copy <- original
    copy[names(copy)] <- 0
    copy[names(replication)] <- replication
    copy
  }
  # Bootstrap variables
  bootRepdf <- list()
  as_tibble <- tibble::as_tibble
  left_join <- dplyr::left_join
  tibble <- tibble::tibble
  bootN <- dim(df_full)[1]
  # Hamming and output variables
  origCos <- calc_cos_sim(df)
  origMem1 <- membership(calc_com_louv(origCos$col1_cossim)$community)
  origMem2 <- membership(calc_com_louv(origCos$col2_cossim)$community)
  col1Diff <- rep(0,nBoot)
  col2Diff <- rep(0,nBoot)
  col1Change <- matrix(0, length(origMem1), (nBoot+1))
  col1Change <- as_tibble(col1Change)
  col1Change[,1] <- names(origMem1)
  names(col1Change) <- c(col1,paste0("Rep",1:nBoot))
  col2Change <- matrix(0, length(origMem2), (nBoot+1))
  col2Change <- as_tibble(col2Change)
  col2Change[,1] <- names(origMem2)
  names(col2Change) <- c(col2,paste0("Rep",1:nBoot))
  # Bootstrapping
  for (i in 1:nBoot){
    if (method=='single'){
      bootInd <- sample(1:bootN,bootN,replace=TRUE)
      thisBoot <- df_full[bootInd,]
    } else{
      bootInd1 <- sample(1:bootN,bootN,replace=TRUE)
      bootInd2 <- sample(1:bootN,bootN,replace=TRUE)
      thisBoot1 <- df_full[bootInd1,names(df_full)==col1]
      thisBoot2 <- df_full[bootInd2,names(df_full)==col2]
      thisBoot <- dplyr::bind_cols(thisBoot1,thisBoot2)
    }
    dfBoot <- create_graph_df(thisBoot, col1, col2)
    dfBoot <- as_tibble(dfBoot$df)
    bootRepdf <- append(bootRepdf, dfBoot)
    noiseCos <- calc_cos_sim(dfBoot)
    col1_cos_noise = calc_com_louv(noiseCos$col1_cossim)
    col2_cos_noise = calc_com_louv(noiseCos$col2_cossim)
    col1mem <- membership(col1_cos_noise$community)
    col2mem <- membership(col2_cos_noise$community)
    swaps <- combinat::permn(max(col1mem))
    swaps2 <- combinat::permn(max(col2mem))
    dif <- rep(0,length(swaps))
    dif2 <- rep(0,length(swaps2))
    flag1 <- length(col1mem) != length(origMem1)
    flag2 <- length(col2mem) != length(origMem2)
    for (k in 1:length(swaps)){
      temp <- swaps[[k]]
      col1tempMem <- temp[col1mem]
      names(col1tempMem) <- names(col1mem)
      if(flag1){
        col1tempMem <- copy_and_add_list(origMem1, col1tempMem)
      }
      dif[k] <- sum(origMem1 != col1tempMem)
    }
    for (k in 1:length(swaps2)){
      temp <- swaps2[[k]]
      col2tempMem <- temp[col2mem]
      names(col2tempMem) <- names(col2mem)
      if(flag2){
        col2tempMem <- copy_and_add_list(origMem2, col2tempMem)
      }
      dif2[k] <- sum(origMem2 != col2tempMem)
    }
    col1Diff[i] <- min(dif)
    col2Diff[i] <- min(dif2)
    t1 <- swaps[[which.min(dif)]][col1mem]
    names(t1) <- names(col1mem)
    t2 <- swaps2[[which.min(dif2)]][col2mem]
    names(t2) <- names(col2mem)
    if (flag1){
      t1 <- copy_and_add_list(origMem1, t1)
    }
    if (flag2){
      t2 <- copy_and_add_list(origMem2, t2)
    }
    col1Change[,i+1] <- as.integer(t1 != origMem1)
    col2Change[,i+1] <- as.integer(t2 != origMem2)
    if (i==1){
      bootRepComm1 <- tibble(node_id = names(t1), community=t1)
      bootRepComm2 <- tibble(node_id = names(t2), community=t2)
    } else{
      repTemp <- tibble(node_id = names(t1), community=t1)
      repTemp2 <- tibble(node_id = names(t2), community=t2)
      bootRepComm1 <- left_join(bootRepComm1, repTemp, by="node_id")
      bootRepComm2 <- left_join(bootRepComm2, repTemp2, by="node_id")
    }
  }
  names(bootRepComm1) <- c(col1,paste0("Rep",1:nBoot))
  names(bootRepComm2) <- c(col2,paste0("Rep",1:nBoot))
  output <- list(comm1 = bootRepComm1, comm2 = bootRepComm2,
                 hamming1 = col1Diff, hamming2 = col2Diff,
                 col1Changes = col1Change, col2Changes = col2Change,
                 col1Stats = list(Mean=mean(col1Diff), STD=sd(col1Diff)),
                 col2Stats = list(Mean=mean(col2Diff), STD=sd(col2Diff)),
                 df = bootRepdf)
}

create_bootstrap_RCA_communities <- function(df, df_full, col1, col2, nBoot=100, method='single'){
  # Function to deal with missing nodes during resampling
  copy_and_add_list <- function(original, replication){
    copy <- original
    copy[names(copy)] <- 0
    copy[names(replication)] <- replication
    copy
  }
  # Bootstrap variables
  bootRepdf <- list()
  as_tibble <- tibble::as_tibble
  left_join <- dplyr::left_join
  tibble <- tibble::tibble
  bootN <- dim(df_full)[1]
  # Hamming and output variables
  origRCA <- calc_rca_sim(df, col1, col2)
  origMem1 <- membership(calc_com_louv(origRCA$RCA_df_col)$community)
  origMem2 <- membership(calc_com_louv(origRCA$RCA_df_row)$community)
  col1Diff <- rep(0,nBoot)
  col2Diff <- rep(0,nBoot)
  col1Change <- matrix(0, length(origMem1), (nBoot+1))
  col1Change <- as_tibble(col1Change)
  col1Change[,1] <- names(origMem1)
  names(col1Change) <- c(col1,paste0("Rep",1:nBoot))
  col2Change <- matrix(0, length(origMem2), (nBoot+1))
  col2Change <- as_tibble(col2Change)
  col2Change[,1] <- names(origMem2)
  names(col2Change) <- c(col2,paste0("Rep",1:nBoot))
  # Bootstrapping
  for (i in 1:nBoot){
    if (method=='single'){
      bootInd <- sample(1:bootN,bootN,replace=TRUE)
      thisBoot <- df_full[bootInd,]
    } else{
      bootInd1 <- sample(1:bootN,bootN,replace=TRUE)
      bootInd2 <- sample(1:bootN,bootN,replace=TRUE)
      thisBoot1 <- df_full[bootInd1,names(df_full)==col1]
      thisBoot2 <- df_full[bootInd2,names(df_full)==col2]
      thisBoot <- dplyr::bind_cols(thisBoot1,thisBoot2)
    }
    dfBoot <- create_graph_df(thisBoot, col1, col2)
    dfBoot <- as_tibble(dfBoot$df)
    bootRepdf <- append(bootRepdf, dfBoot)
    noiseRCA <- calc_rca_sim(dfBoot, col1, col2)
    col1_rca_noise = calc_com_louv(noiseRCA$RCA_df_col)
    col2_rca_noise = calc_com_louv(noiseRCA$RCA_df_row)
    col1mem <- membership(col1_rca_noise$community)
    col2mem <- membership(col2_rca_noise$community)
    swaps <- combinat::permn(max(col1mem))
    swaps2 <- combinat::permn(max(col2mem))
    dif <- rep(0,length(swaps))
    dif2 <- rep(0,length(swaps2))
    flag1 <- length(col1mem) != length(origMem1)
    flag2 <- length(col2mem) != length(origMem2)
    for (k in 1:length(swaps)){
      temp <- swaps[[k]]
      col1tempMem <- temp[col1mem]
      names(col1tempMem) <- names(col1mem)
      if(flag1){
        col1tempMem <- copy_and_add_list(origMem1, col1tempMem)
      }
      dif[k] <- sum(origMem1 != col1tempMem)
    }
    for (k in 1:length(swaps2)){
      temp <- swaps2[[k]]
      col2tempMem <- temp[col2mem]
      names(col2tempMem) <- names(col2mem)
      if(flag2){
        col2tempMem <- copy_and_add_list(origMem2, col2tempMem)
      }
      dif2[k] <- sum(origMem2 != col2tempMem)
    }
    col1Diff[i] <- min(dif)
    col2Diff[i] <- min(dif2)
    t1 <- swaps[[which.min(dif)]][col1mem]
    names(t1) <- names(col1mem)
    t2 <- swaps2[[which.min(dif2)]][col2mem]
    names(t2) <- names(col2mem)
    if (flag1){
      t1 <- copy_and_add_list(origMem1, t1)
    }
    if (flag2){
      t2 <- copy_and_add_list(origMem2, t2)
    }
    col1Change[,i+1] <- as.integer(t1 != origMem1)
    col2Change[,i+1] <- as.integer(t2 != origMem2)
    if (i==1){
      bootRepComm1 <- tibble(node_id = names(t1), community=t1)
      bootRepComm2 <- tibble(node_id = names(t2), community=t2)
    } else{
      repTemp <- tibble(node_id = names(t1), community=t1)
      repTemp2 <- tibble(node_id = names(t2), community=t2)
      bootRepComm1 <- left_join(bootRepComm1, repTemp, by="node_id")
      bootRepComm2 <- left_join(bootRepComm2, repTemp2, by="node_id")
    }
  }
  names(bootRepComm1) <- c(col1,paste0("Rep",1:nBoot))
  names(bootRepComm2) <- c(col2,paste0("Rep",1:nBoot))
  output <- list(comm1 = bootRepComm1, comm2 = bootRepComm2,
                 hamming1 = col1Diff, hamming2 = col2Diff,
                 col1Changes = col1Change, col2Changes = col2Change,
                 col1Stats = list(Mean=mean(col1Diff), STD=sd(col1Diff)),
                 col2Stats = list(Mean=mean(col2Diff), STD=sd(col2Diff)),
                 df = bootRepdf)
}
