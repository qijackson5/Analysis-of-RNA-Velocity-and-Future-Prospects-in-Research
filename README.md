# Analysis of RNA Velocity and Future Prospects in Research

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/UCB-stat-159-s22/hw07-hw07-group12/main?labpath=finalproject%2Ffig1.ipynb)

**Analysis Jupyter Book Link:** https://ucb-stat-159-s22.github.io/hw07-hw07-group12/intro.html

Note: We also did analysis on figure 2 and utilized the utils package we created. This is not shown in the jupyterbook and the binder, but it is in the main.ipynb. Feel free to give it a look! 

To run figure 2, clone this repo, open main.ipynb and feel free to run the second half of the main.ipynb. 

## Purpose statement: 

This project is a conceptual overview of RNA velocity and an exploratory analysis of of a dentate gyrus dataset. 

RNA velocity is a computational method that, given a single cell RNA-seq dataset, predicts, for each cell, the future state of that cell. Analogous to a photograph with motion blur, it is a method that leverages a static “snapshot” of a cell to infer its trajectory, by comparing relative levels of spliced and unspliced mRNA. 

The analyses here explore a dentate gyrus dataset by evaluating the overall velocities of cells within annotated cell clusters, as well as exploring the transcriptional dynamics across all cells for specific marker genes. Additionally, multiple methods of inferring RNA velocity are compared. The sometimes divergent biological conclusions across models underscore the importance of model selection.

## Data accessibility: 

All of the data can be accessed through (https://scvelo.readthedocs.io/) and (https://www.nature.com/articles/s41587-020-0591-3), including a tutorial on how to use the tool kit.  

## Installation Instructions:

- To create the environment `dev_env`, run `make env` 
    - env: sets up the environment 
    - clean: removes the figures already in the folder
    - all: runs 3 analyses notebooks and 1 main notebook
- Run `make clean` to remove any already loaded figures
- Run `make all` to run the 3 analysis notebooks in the analysis_notebook directory and the main.ipynb notebook
- In order to run the environment successfully and to test the functions, after running `conda activate dev_env`, run `conda install -c anaconda pytest`.

## Further clarifications of directories:

- The `finalproject` directory is for the JupyterBook
- Tests fall under the `scvtools` directory
- `analysis_notebook` directory contains all notebooks with analyses

Authorship of the original built-in data rests with https://scvelo.readthedocs.io/.
