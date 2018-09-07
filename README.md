# README #

### What is this repository for? ###

This repository contains a workflow for the analysis fo single-cell RNASeq data using R/bioconductor. 
The main steps included are:

* Quality control
* Data normalization
* Feature selection
* Clustering
* Differential expression analysis
* Visualization
* Identification of rare cell subtypes with CellSIUS

### Software requirements ###

This workflow is designed for and tested on Linux systems, but should run on other UNIX flavors as well. The workflow was built and tested using R 3.4.1 and Bioconductor 3.5. Because some of the packages it uses (especially scater) changed between Bioconductor 3.5 and 3.6, currently, the only supported version is R3.4.1 with Bioconductor 3.5.

### How do I get set up? ###

Check the vignette (vignettes/workflow_vignette.html) for a detailed description. In brief, here's what you need:

1. This repository, cloned somewhere you have both read and write permissions.
2. R version 3.4.1 and (optionally but highly recommended) RStudio.
3. A series of R/bioconductor packages. Please refer to the vignette for details. Downloading and installing these packages might take some time. If you are starting from scracth without anything pre-installed, this will take ~ 4 hours.
4. An installation of the markov clustering algirthm (MCL, Stijn van Dongen, Graph Clustering by Flow Simulation. PhD thesis, University of Utrecht, May 2000), which can be found [here](https://micans.org/mcl/).

All code can be run directly form the R-Markdown document (vignettes/workflow_vignette.Rmd), using the "vignettes" folder as the working directory. We also provide a small test dataset which was used to generate the vignette. Once you are set up, running the full workflow on this test dataset takes ~5-10 minutes on a standard desktop computer.

### Need help? ###

Please contact Rebekka Wegmann [wegmann@imsb.biol.ethz.ch](mailto:wegmann@imsb.biol.ethz.ch) or Marilisa Neri [marilisa.neri@novartis.com](mailto:marilisa.neri@novartis.com)