[![DOI](https://zenodo.org/badge/190090979.svg)](https://zenodo.org/badge/latestdoi/190090979)

# Introduction

This is the public data repository for:

[David A. Koweek, Clara García-Sánchez, Philip G. Brodrick, Parker Gassett, and Ken Caldeira.
Evaluating hypoxia alleviation through induced downwelling. Science of the Total Environment,
719:137334, 2020.](https://doi.org/10.1016/j.scitotenv.2020.137334)    

You can use this repository to reproduce any of the analyses in the paper. The analyses necessary to produce each figure are contained in the script under that figure's name. For instance, to reproduce figure 6 in the paper, simply run:

`source("build_figure_6.R")`

The resulting figure will be produced and deposited in the `figures` folder as `figure_6.pdf`.

Necessary packages for each figure are defined at the beginning of each script. If you don't have the necessary packages, install the packages before executing the script.

The `data` folder contains data sets from the field experiment and the publicly available hydrocasts used for the oxygen transfer efficiency modeling. Calls to the data sets necessary to produce a given figure are found at the beginning of the script to produce that figure. 

This repository contains several supporting scripts, many of which contain custom functions necessary to produce the analyses and figures. Again, these scripts are sourced at the beginning of the appropriate figure script.

# Details

Users interested in re-creating the experimental results (O<sub>2</sub>') in `Fig. 5` of the paper can run the `Python` script `generate_temporal_interpolations.py`

The schematic figures of the paper (e.g. Figure 1) have been added to the `figures` folder for interested users.

Figure S1, a representation of Figure 3 on a log-scale, can be re-created by sourcing the script to build Figure 3 (`source("build_figure_3.R")`).

If you have any questions about accessing the material in this repository and/or find an error, please do not hesitate to contact me.

David Koweek
