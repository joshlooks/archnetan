#' Adds noise (in the form of full data points), applies RCA and then create cosine-similarity communities from an inputted dataset
#'
#' @param df Dataframe original frequency graph
#' @param col1 A string of the first column used to make graphs
#' @param col2 A string of the second column used to make graphs
#' @param proportion A double corresponding to the proportion of the original dataset's size to add as noise
#' @param numReps Number of runs noise addition runs to complete
#' @returns An output list containing nine components:
#' \itemize{
#' \item{comm1}{Tibble containing the assigned communities of column 1 in each run}
#' \item{comm2}{Tibble containing the assigned communities of column 2 in each run}
#' \item{hamming1}{Vector containing the Hamming distance to the original column 1 communities in each run}
#' \item{hamming2}{Vector containing the Hamming distance to the original column 1 communities in each run}
#' \item{col1Changes}{Tibble of binary values for whether a column 1 node changes communities in each run}
#' \item{col2Changes}{Tibble of binary values for whether a column 1 node changes communities in each run}
#' \item{df}{List of lists corresponding to the df_graph like tibbles produced in each bootstrap}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' RCA_communities_single(df_graph, "Source", "Site", proportion=0.2, numReps=50)
#' }
RCA_communities_single <- function(df, col1, col2, proportion=0.10, numReps=100){
  as_tibble <- tibble::as_tibble
  tibble <- tibble::tibble
  # Bootstrap variables
  df <- as_tibble(df)
  nRows <- dim(df)[1]
  numNoise <- ceiling(nRows*proportion)
  left_join <- dplyr::left_join
  noiseRepdf <- list()
  # Hamming and output variables
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
  # Bootstrapping added noise
  for (i in 1:numReps){
    # Drawing with variable concordance/draw a full datum
    dfRep <- df
    noiseInd <- sample(1:nRows,numNoise,replace=TRUE,prob=wts)
    for (j in noiseInd){
      dfRep[j,3] <- dfRep[j,3]+1
    }
    # create the graph and communities
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
    # checking for relabelling of original communities
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
    # actual change should be minimum of re-ordering changes (unless all resemblance is broken)
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
  # construct output list
  names(noiseRepComm1) <- c(col1,paste0("Rep",1:(dim(noiseRepComm1)[2]-1)))
  names(noiseRepComm2) <- c(col2,paste0("Rep",1:(dim(noiseRepComm2)[2]-1)))
  output <- list(comm1 = noiseRepComm1, comm2 = noiseRepComm2,
                 hamming1 = col1Diff, hamming2 = col2Diff,
                 col1Changes = col1Change, col2Changes = col2Change,
                 df = noiseRepdf)
}

#' Adds noise (by drawing each variable independently), applies RCA and then create cosine-similarity communities from an inputted dataset
#'
#' @param df Dataframe original frequency graph
#' @param col1 A string of the first column used to make graphs
#' @param col2 A string of the second column used to make graphs
#' @param proportion A double corresponding to the proportion of the original dataset's size to add as noise
#' @param numReps Number of runs noise addition runs to complete
#' @returns An output list containing nine components:
#' \item{comm1}{Tibble containing the assigned communities of column 1 in each run}
#' \item{comm2}{Tibble containing the assigned communities of column 2 in each run}
#' \item{hamming1}{Vector containing the Hamming distance to the original column 1 communities in each run}
#' \item{hamming2}{Vector containing the Hamming distance to the original column 1 communities in each run}
#' \item{col1Changes}{Tibble of binary values for whether a column 1 node changes communities in each run}
#' \item{col2Changes}{Tibble of binary values for whether a column 1 node changes communities in each run}
#' \item{df}{List of lists corresponding to the df_graph like tibbles produced in each bootstrap}
#' @export
#'
#' @examples
#' \dontrun{
#' RCA_communities_multi(df_graph, "Source", "Site", proportion=0.2, numReps=50)
#' }
RCA_communities_multi <- function(df, col1, col2, proportion=0.10, numReps=100){
  as_tibble <- tibble::as_tibble
  tibble <- tibble::tibble
  df <- as_tibble(df)
  nRows <- dim(df)[1]
  numNoise <- ceiling(nRows*proportion)
  left_join <- dplyr::left_join
  noiseRepdf <- list()
  # Hamming and output variables
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
  # Bootstrapping added noise
  for (i in 1:numReps){
    # Drawing each variable individually
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
    # create the graph and communities
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
    # checking for relabelling of original communities
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
    # actual change should be minimum of re-ordering changes (unless all resemblance is broken)
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

#' Runs RCA_communities_single or RCA_communities multiple over a range of proportions and aggregates Hamming distances
#'
#' @param df Dataframe original frequency graph
#' @param col1 A string of the first column used to make graphs
#' @param col2 A string of the second column used to make graphs
#' @param proportion A double corresponding to the proportion of the original dataset's size to add as noise
#' @param numReps Number of runs noise addition runs to complete
#' @param method A string determining the method of drawing added noise data points:
#' \itemize{
#' \item{single (default)}{Draws a full artifact from the df dataframe}
#' \item{multi}{Draws each characteristic from the df dataframe}
#' }
#' @param propSave Proportion of noise added for which to save the cos_communities_single/multi output for
#'
#' @returns An output list containing three components:
#' \item{col1Stats}{Tibble containing the mean and standard deviation in column 1 in Hamming distance for each proportion}
#' \item{col1Stats}{Tibble containing the mean and standard deviation in column 2 in Hamming distance for each proportion}
#' \item{propSave}{List in same structure as outputs from cos_communities_single/multi for chosen proportion value}
#' @export
#'
#' @examples
#' \dontrun{
#' create_added_noise_RCA_communities(df_graph, "Source", "Site",
#' proportion=c(0.2,0.3,0.4,0.5), numReps=50, method="multi", propSave=0.5)
#' }
create_added_noise_RCA_communities <- function(df, col1, col2, proportion=0.10, numReps=100, method='single', propSave=0.10){
  # Choose method of drawing random noise
  if (method=='single'){
    RCAc <- RCA_communities_single
  } else {
    RCAc <- RCA_communities_multi
  }
  # Choose output based on whether a vector or single proportion is chosen
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
    # Run for range of proportions and save communities dataframes for selected proportion
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
