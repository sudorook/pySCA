# pySCA

![Website Build Status](https://gitlab.com/sudorook/pysca/badges/master/pipeline.svg)

> 10.2020
>
> Copyright (C) 2020 Olivier Rivoire, Rama Ranganathan, and Kimberly Reynolds,
> and Ansel George
>
> This program is free software distributed under the BSD 3-clause license,
> please see the file LICENSE for details.

The current version of the Statistical Coupling Analysis (SCA) analysis is
implemented in Python. This directory contains the necessary code for running
the SCA calculations, as well examples/tutorials for the dihydrofolate
reductase (DHFR) enzyme family, the S1A serine proteases, the small G-protein
family and the Beta-lactamase enzyme family. The tutorials are distributed as
Jupyter notebooks; for details please see:
[https://jupyter.org/](https://jupyter.org/).

For installation instructions, and an introduction to using the toolbox, please
[refer to the website](https://sudorook.gitlab.io/pysca).

or look through the [RST files](docs/source) included with the pySCA
distribution.

## Contents of `/`

|            |                                                         |
| :---       | :---                                                    |
| bin/       | Executables for running SCA analysis functions          |
| data/      | Input data (including those needed for the tutorials)   |
| docs/      | HTML documentation (generated by Sphinx)                |
| figs/      | Figures used for the notebooks and documentation        |
| notebooks/ | Example SCA notebooks                                   |
| output/    | Output files (empty at install, use `runAllNBCalcs.sh`) |
| pysca/     | Python code for SCA                                     |
| scripts/   | Utility scripts used to generate example data           |

## Contents of `bin/`

|               |                                                                  |
| :---          | :---                                                             |
| annotateMSA   | Annotates alignments with phylogenetic/taxonomic information     |
| scaProcessMSA | Conducts some initial processing of the sequence alignment       |
| scaCore       | Runs the core SCA calculations                                   |
| scaSectorID   | Defines sectors given the results of the calculations in scaCore |

## Contents of `pysca/`

|             |                                                      |
| :---        | :---                                                 |
| scaTools.py | The SCA toolbox - functions for the SCA calculations |
| settings.py | Global configuration settings for the analysis       |

## Contents of `notebooks/`

|                         |                                     |
| :---                    | :---                                |
| SCA_DHFR.ipynb          | Example for DHFR                    |
| SCA_G.ipynb             | Example for the small G proteins    |
| SCA_betalactamase.ipynb | Example for the beta-lactamases     |
| SCA_S1A.ipynb           | Example for the S1A serine protease |
