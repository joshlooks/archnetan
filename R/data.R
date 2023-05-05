#' obsidian_graph_df
#'
#' Frequency weighted bipartite graph of New Zealand obsidian data collected by collaborators
#' and package co-authors
#'
#' @format A data frame with three variables:
#' \describe{
#' \item{\code{Study Site}}{Assemblage location where the obsidian artifact was found}
#' \item{\code{Source}}{Geographical origin of the obsidian artifact}
#' \item{\code{weight}}{Edgeweight in the corresponding bipartite graph of Study Site and Source}
#' }
#'
#' @source <https://doi.org/10.5334/joad.52>
"obsidian_graph_df"

#' obsidian_full_df
#'
#' Dataframe of New Zealand obsidian data collected by collaborators
#' and package co-authors
#'
#' @format A data frame with seven variables:
#' \describe{
#' \item{\code{Study Site}}{Assemblage location where the obsidian artifact was found}
#' \item{\code{Source}}{Geographical origin of the obsidian artifact}
#' \item{\code{Period Start}}{Start of the likely time period where the artifact was deposited}
#' \item{\code}{Period End}}{End of the likely time period where the artifact was deposited}
#' \item{\code}{Probable Period}{Whether the likely time period is 'early' or 'late'}
#' \item{\code}{Weight (g)}{Weight of the artifact in grams}
#' \item{\code}{Reference}{Citation for the original reference of the artifact}
#' }
#'
#' @source <https://doi.org/10.5334/joad.52>
"obsidian_full_df"
