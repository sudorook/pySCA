.. SPDX-FileCopyrightText: 2019 sudorook <daemon@nullcodon.com>
..
.. SPDX-License-Identifier: BSD-3-Clause

=====
Usage
=====

We provide tutorials that walk through the process of sector identification for
three protein families: the Ras-like small G-proteins, the metabolic enzyme
Dihydrofolate Reductase (DHFR), and the antibiotic resistance enzyme
Beta-lactamase.

To run the SCA calculations for all three examples, you can execute the
following shell script from the scripts/ directory::

  ./runAllNBCalcs.sh

For each example, this will generate the following outputs in the output/
directory:

  1.  A pickle database (\*.db file) that contains the results of the
      calculations (these are then read in and analyzed in the IPython
      notebooks - \*.ipynb)
  2.  A \*.log file that provides some information about the analysis
  3.  A processed alignment (\*.fasta file) resulting from the
      scaProcessMSA script.

Following this step, you can begin the tutorial as an interactive Jupyter
notebook from the command line as follows::

  jupyter notebook SCA_G.ipynb

This should open the notebook in a browser window, where you can run the code,
and examine the SCA results.
