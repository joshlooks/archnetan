## Prepares rda files containing functional examples run on the included dataset.

example_data <- archnetan::create_graph_single_file(archnetan::archnetan_example("originalDataset.csv"),
                                                    "Study Site", "Source")
usethis::use_data(DATASET, overwrite = TRUE)
