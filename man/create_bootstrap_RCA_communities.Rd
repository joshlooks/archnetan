% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/booststrap_robustness_RCA.R
\name{create_bootstrap_RCA_communities}
\alias{create_bootstrap_RCA_communities}
\title{Apply RCA and then create cosine-similarity communities from bootstrapped data}
\usage{
create_bootstrap_RCA_communities(
  df,
  df_full,
  col1,
  col2,
  nBoot = 100,
  method = "single"
)
}
\arguments{
\item{df}{Dataframe original frequency graph}

\item{df_full}{Dataframe used to build original frequency graph/of raw data}

\item{col1}{A string of the first column used to make graphs}

\item{col2}{A string of the second column used to make graphs}

\item{nBoot}{Number of bootstrap runs to complete}

\item{method}{A string determining the method of drawing bootstrap replicants: single for full artifact,
multi for drawing each artifact variably separately}
}
\value{
An output list containing nine components:
\item{comm1}{Tibble containing the assigned communities of column 1 in each bootstrap run}
\item{comm2}{Tibble containing the assigned communities of column 2 in each bootstrap run}
\item{hamming1}{Vector containing the Hamming distance to the original column 1 communities in each bootstrap run}
\item{hamming2}{Vector containing the Hamming distance to the original column 1 communities in each bootstrap run}
\item{col1Changes}{Tibble of binary values for whether a column 1 node changes communities in each bootstrap run}
\item{col2Changes}{Tibble of binary values for whether a column 1 node changes communities in each bootstrap run}
\item{col1Stats}{List of the mean and standard deviation of the Hamming distance for column 1 across bootstrap runs}
\item{col2Stats}{List of the mean and standard deviation of the Hamming distance for column 2 across bootstrap runs}
\item{df}{List of lists corresponding to the df_graph like tibbles produced in each bootstrap}
}
\description{
Apply RCA and then create cosine-similarity communities from bootstrapped data
}
\examples{
\dontrun{
create_bootstrap_RCA_communities(df_graph, df_full, "Source", "Site", nboot=200, method="multi")
}
}
