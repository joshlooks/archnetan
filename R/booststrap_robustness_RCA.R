#' Title
#'
#' @param df
#' @param df_full
#' @param col1
#' @param col2
#' @param nBoot
#' @param method
#'
#' @return
#' @export
#'
#' @examples
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
