This is the public data repository for:

*Alleviating hypoxia through induced downwelling* (final citation to be updated after manuscript acceptance)

You can use this repository to reproduce any of the analyses in the paper. The analyses necessary to produce each figure are contained in the script under that figure's name. For instance, to reproduce figure 6 in the paper, simply run:

`source("build_figure_6.R")`

The resulting figure will be produced and deposited in the `figures` folder as `figure_6.pdf`.

Necessary packages for each figure are defined at the beginning of each script. If you don't have the necessary packages, install the packages before executing the script.

The `data` folder contains data sets from the field experiment and some of the publicly available hydrocasts used for the oxygen transfer efficiency modeling. Calls to the data sets necessary to produce a given figure are found at the beginning of the script to produce that figure. 

This repository contains several supporting scripts, many of which contain custom functions necessary to produce the analyses and figures. Again, these scripts are sourced at the beginning of the appropriate figure script.

Users interested in re-creating the experimental results ($O_2'$) in `Fig. 5` of the paper can run the `Python` script `generate_temporal_interpolations.py`

Some of the schematic figures of the paper (e.g. Figure 1) have been added to the `figures` folder for interested users.

If you have any questions about accessing the material in this repository and/or find an error, please do not hesitate to contact me.

David Koweek