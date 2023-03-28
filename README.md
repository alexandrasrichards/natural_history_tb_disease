# Natural History of TB Disease

Code used for the analysis in the Lancet Global Health paper

[Quantifying progression and regression across the spectrum of pulmonary tuberculosis: a data synthesis study](https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(23)00082-7/fulltext)


## Libraries
The calibration uses [LibBi](https://libbi.org/) (which needs an external download) and packages rbi and rbi.helpers. In order for the rbi.helpers package to run correctly on the data used here, the following line needs to be run after initial installation of rbi and rbi.helpers:
`remotes::install_github("sbfnk/rbi@sparse_input")`. This may not work on Windows computers, but the data from an example run has been saved for use in the rest of the analysis.

Other libraries required to run this analysis are:

* tidyverse
* patchwork
* png
* binom
* bayesplot
* ggpubr
* reshape2

## Folder structure and code
A script for the calibration code is contained in the Calibration folder, with the files contained in the calibration sub-folders all sourced within the script. Numerical outputs from the calibration are stored in the Results folder, and the figure published in the paper is stored in the Figures folder. Other figures are created within the script but not saved.

The post-calibration cohort model analysis code is contained in the Analysis folder, with no sub-folders. The file `analysis_script.R` sources both other files. The `median_trajectories` function takes a long time with `nreps = 1000` and `npeople = 10000`, change these at the top of the file for a quicker analysis. The results will be saved in the Analysis subfolder of the Results folder.

