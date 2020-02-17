# ABLS2020-Poster

The poster for the Applied Bioinformatics in Life Sciences 2020 conference. The poster was presented on the February 14, 2020. It is available as a PDF file (`Poster.pdf`) in the main directory of this repository, it can be downloaded from [here](https://doi.org/10.5281/zenodo.3670061) or the poster can be generated from the latex file.

## Installation

In order to reproduce the poster you will need:

* A working `R` installation with the required packages also installed (see start of the scripts). Note that the package `Features` and `MsCoreUtils` are currenly (as of 14/02/2020) only available from GitHub (see [here](https://www.rformassspectrometry.org/)). The easiest way is to install it using `devtools::install_github("rformassspectrometry/Features")`.
* A working `TeX` distribution

## Compiling poster

The repository contains the required files to generate the poster. The poster is generated by compiling the `Poster.tex` file. The figures used in the poster are generated from the R scripts `poster_figs.R` and `R/utils.R`.

## Reference

You can cite the poster as:

> Vanderaa, Christophe, & Gatto, Laurent. (2020). Towards a standardized workflow for mass spectrometry-based single-cell proteomics. Zenodo. http://doi.org/10.5281/zenodo.3670061
