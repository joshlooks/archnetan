# ======================================================
# New Zealand Obsidian data
# This code prepares the dataset found in DOI: <https://doi.org/10.5334/joad.52>, collected
# by collaborators and co-authors of this package, for usage in demonstrating
# package functionality and code-tests
# Created: 05/05/2023

# ======================================================

graph_info <- archnetan::create_graph_single_file(archnetan::archnetan_example("originalDataset.csv"),
                                                  col1="Study Site",col2="Source")
obsidian_graph_df <- graph_info$df
obsidian_full_df <- readr::read_csv(archnetan::archnetan_example("originalDataset.csv"))

usethis::use_data(obsidian_graph_df, obsidian_full_df, overwrite = TRUE)
