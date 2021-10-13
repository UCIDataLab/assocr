This R package contains code for assessing the likelihood two temporal event series were generated by the same source. Two techniques are implemented
1. Score-based likelihood ratio
2. Coincidental match probability

Both methods rely on score functions for quantifying the similarity between two temparal event series. Currently this package supports the coefficient of segregation, mingling index, and mean and median inter-event time as defined as the time to the nearest neighbor from one series to another.

More details can be found in the following paper:

[Galbraith, C., Smyth, P. and Stern, H.S., 2020. Quantifying the association between discrete event time series with applications to digital forensics. *Journal of the Royal Statistical Society: Series A (Statistics in Society).*](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssa.12549)

To re-create the experiments in the paper, download the [zip file from the RSS](https://rss.onlinelibrary.wiley.com/pb-assets/hub-assets/rss/Datasets/RSSA%20183.3/A1549Galbraith-1591284720010.zip). This folder contains the input data and scripts used to generate all results in the paper. The scripts rely on having the `assocr` package installed, with can be accomplished with the [`devtools::install_github()` function](https://github.com/r-lib/devtools).

This work was supported by the [Center for Statistics and Applications in Forensic Evidence (CSAFE)](https://forensicstats.org/).
