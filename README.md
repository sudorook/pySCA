# pySCA

![Website Build Status](https://gitlab.com/sudorook/pySCA/badges/master/build.svg)

> 03.2015 
>
> Copyright (C) 2015 Olivier Rivoire, Rama Ranganathan, Kimberly Reynolds
>
> This program is free software distributed under the BSD 3-clause license,
> please see the file LICENSE for details.

The current version of the Statistical Coupling Analysis (SCA) analysis is
implemented in python. This directory contains the necessary code for running
the SCA calculations, as well examples/tutorials for the dihydrofolate
reductase (DHFR) enzyme family, the S1A serine proteases, the small G-protein
family and the Beta-lactamase enzyme family. The tutorials are distributed as
Jupyter notebooks; for details please see:
[https://jupyter.org/](https://jupyter.org/).

For installation instructions, and an introduction to using the toolbox, please
see:

[https://gitlab.com/sudorook/pySCA](https://gitlab.com/sudorook/pySCA)

or open the html files included with the pySCA distribution in a browser:

[docs/index.html](docs/index.html)

## Contents of `/`

|            |                                                         |
| :---       | :---                                                    |
| data/      | Input data (including those needed for the tutorials)   |
| output/    | Output files (empty at install, use `runAllNBCalcs.sh`) |
| docs/      | HTML documentation (generated by Sphinx)                |
| scripts/   | Utility scripts used to generate example data           |
| pysca/     | Python code for SCA                                     |
| notebooks/ | Example SCA notebooks                                   |

## Contents of `pysca/`

|                  |                                                                          |
| :---             | :---                                                                     |
| annotate_MSA.py  | Annotates alignments with phylogenetic/taxonomic information             |
| scaProcessMSA.py | Conducts some initial processing of the sequence alignment               |
| scaCore.py       | Runs the core SCA calculations                                           |
| scaSectorID.py   | Defines sectors given the results of the calculations in scaCore         |
| scaTools.py      | The SCA toolbox - contains all functions needed for the SCA calculations |

## Contents of `notebooks/`

|                         |                                     |
| :---                    | :---                                |
| SCA_DHFR.ipynb          | Example for DHFR                    |
| SCA_G.ipynb             | Example for the small G proteins    |
| SCA_betalactamase.ipynb | Example for the beta-lactamases     |
| SCA_S1A.ipynb           | Example for the S1A serine protease |
