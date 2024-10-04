.. SPDX-FileCopyrightText: 2019 sudorook <daemon@nullcodon.com>
..
.. SPDX-License-Identifier: BSD-3-Clause

===============
Getting Started
===============

Running a complete SCA analysis consists of five steps:

  1) Constructing an alignment
  2) Alignment pre-processing and conditioning
  3) Calculation of the conservation and co-evolution statistics
  4) Identifying statistically significant correlations
  5) Interpretation of the results

The core SCA calculations (steps 2,3, and 4) are each associated with a
particular Python analysis script (:code:`scaProcessMSA`, :code:`scaCore`, and
:code:`scaSectorID`, respectively). Sequential execution of each Python
analysis script stores the results in a pickle database. This means the core
SCA calculations can be run from the command line, or multiple proteins can be
analyzed using a shell script (for an example, see `runAllNBCalcs.sh`).
Following execution of the scripts, the pickle database can be loaded in an
Jupyter notebook for visualizing the results and interpreting the data.
Alternatively, the output of the analysis scripts can be saved as a MATLAB
workspace, and results plotted/analyzed in MATLAB. Below we describe the five
main analysis steps in more detail.


File and Directory Structure
============================

The pySCA repository contains the following files and directories:

Base Directory
--------------

bin/
  Contains the analysis scripts that use functions defined in the `scaTools.py`
  module.
data/
  Git submodule that contains the input sequence alignments (\*.fasta, \*.an)
  and structures (\*.pdb) for the analysis. The \*.an files correspond to
  fasta-formatted sequence files with taxonomic annotations. The inputs
  needed for all tutorials are here.
output/
  Contains the output of the analysis. Accordingly, it is an empty directory
  in a newly installed pySCA distribution. Running the scripts below will
  output a processed alignment (\*.fasta or \*.an file) and pickle database
  (\*.db file) to Outputs/. Similarly, if you choose to output results to a
  MATLAB workspace, the resulting \*.mat file will write to this directory.
figs/
  Contains a few figures that are loaded into the tutorials for illustration
  purposes.
docs/
  Contains this documentation.
LICENSE
  This work is distributed under the standard BSD 3-clause open source
  software license.
README.md
  Very basic introduction to the toolbox.
scripts/
  Contains scripts used to generate the input data from the example analyses.
notebooks/
  Contains a set of pySCA examples as Jupyter (formerly IPython) notebooks.
pysca/
  Contains the Python source code for the SCA implementation.

`bin` Directory
-----------------

alnFilterSeqSize, alnParseID, alnReplaceHeaders, alnChangeDelim, alnConvertGI
  These aren't essential to the main SCA utilities/package, but are little
  scripts that we often find useful in alignment construction.
annotateMSA
  A script for adding taxonomic annotations to fasta-formatted sequence
  alignments
scaProcessMSA
  The script that conducts alignment pre-processing and conditioning. This
  constitutes trimming the alignment for gaps, and removing low identity
  sequences.
scaCore
  The script that computes SCA conservation and co-evolution values.
scaSectorID
  The script that defines positions that show a statistically significant
  correlation.

`scripts` Directory
-------------------

runAllNBCalcs.sh
  A shell script that runs all of the calculations needed for the tutorials.
  This script also serves as an example for how to call the pySCA scripts
  from the command line.

`notebooks` Directory
---------------------

SCA_DHFR.ipynb
  Jupyter (formerly IPython) notebook tutorial for the Dihydrolate reductase
  enzyme family.
SCA_G.ipynb
  Jupyter notebook tutorial for the small G proteins.
SCA_betalactamase.ipynb
  Jupyter notebook tutorial for the Beta-lactamase enzyme family.
SCA_S1A.ipynb
  Jupyter notebook tutorial for the S1A serine protease enzyme family.

`pysca` Directory
-----------------

scaTools.py
  Contains the pySCA library - the functions that implement all of the SCA
  calculations.
settings.py
  Optional configuration file useful for specifying paths instead of having to
  so do on the command line.


1. Constructing and annotating a multiple sequence alignment
============================================================

The SCA method operates on a multiple sequence alignment of homologous protein
sequences. You can begin the analysis by obtaining an alignment for your
protein of interest from a curated database (for example PFAM:
http://pfam.xfam.org/ ) or by constructing your own alignment. The details of
alignment construction aren't covered here, but we may add a tutorial in future
versions of this documentation. The critical thing is that the alignment
contain on the order of 100 or more effective sequences.

Once you have an alignment, it is helpful to add taxonomic annotations to the
headers. These annotations are used in SCA to examine the relationship between
sector positions and phylogenetic divergence (i.e. in the mapping between
independent components and sequence space). The annotateMSA script contains
two utilities to automate sequence annotation: one which uses the NCBI Entrez
tools in BioPython, and one which uses PFAM database annotations (PFAM
alignment specific). Please note that the annotation step can be slow (on the
order of hours), but only needs to be done once per alignment. For further
details please see the :doc:`/annotateMSA` documentation.

2. Alignment pre-processing and conditioning
============================================

Following alignment construction and annotation, the alignment is processed to:
(1) remove highly gapped or low homology sequences, (2) remove highly gapped
positions, (3) calculate sequence weights and (4) to create a mapping of
alignment positions to a reference structure or sequence numbering system. This
process is handled by the script :doc:`/scaProcessMSA`. Please see the script
documentation for a complete list of optional arguments and notes on usage, and
for a full description of computations 1-4, see the Rivoire et al 2016 methods
paper (Box 1). [#Rivoire2016]_ The resulting output can be stored as either a
Python pickle database or MATLAB workspace for further analysis.

3. Calculation of the conservation and co-evolution statistics
==============================================================

The processed alignment and sequence weights computed in step 2 are then used
in the calculation of evolutionary statistics by the script :doc:`scaCore`.
This script handles the core calculations for:

    1. Pairwise sequence correlations/sequence similarity
    2. Single-site positional conservation from the Kullback-Leibler relative
       entropy, :math:`D_i^a`, and position weights from the gradient of the KL
       entropy, :math:`\frac{\partial{D_i^a}}{\partial{f_i^a}}`. See eqs. 1-2
       in Rivoire, 2016. [#Rivoire2016]_
    3. The SCA matrix :math:`\tilde{C_{ij}}`. See eq. 3 in Rivoire, 2016.
       [#Rivoire2016]_
    4. The projected alignment (eq. 10-11), and the projector (supplemental
       section 1H) [#Rivoire2016]_.
    5. N trials (default N=10) of the randomized SCA matrix and associated
       eigenvectors and eigenvalues; used to choose the number of significant
       eigenmodes.

The calculations and optional execution flags are further described in the
script documentation. As for :doc:`scaProcessMSA`, the output can be stored as
either a Python pickle database or MATLAB workspace for further analysis.

4. Identifying significant evolutionary correlations
====================================================

After the core calculations are complete, the next step is to define the
significant number of eigenmodes/independent components for analysis
(:math:`k_{max}`) and to select sector positions by their contributions to the
top :math:`k_{max}` independent components. This is handled by the script
:doc:`scaSectorID`. This script also computes the sequence-to-position space
mapping as in eq.10-11 and fig. 7. As for :doc:`scaProcessMSA` and
:doc:`scaCore`, the output can be stored as either a Python shelve database or
MATLAB workspace for further analysis.

5. Interpretation of the results and sector definition
======================================================

Execution of annotateMSA, scaProcessMSA, scaCore, and scaSectorID completes
the calculation of SCA terms and results in a single pickle database (\*.db
file, and optionally, a MATLAB workspace) containing the collected results. The
final step is to interpret these calculations and evaluate the
(non-)independence of the amino acid positions associated with each independent
component (as in Fig. 4).

The :doc:`tutorials <usage>` are designed to provide examples of this process,
and to illustrate different aspects of SCA usage (please see the individual
tutorial headers for more information).


**Further Reading/References:**

.. [#Halabi2009] Halabi N, Rivoire O, Leibler S, and Ranganathan R. "Protein
   sectors: evolutionary unis of three-dimensional structure." *Cell.* 2009
   v.138 p.774

.. [#Smock2010] Smock RG, Rivoire O, Russ WP, Swain JF, Leibler S, Ranganathan
   R, Gierasch LM. "An interdomain sector mediating allostery in Hsp70
   molecular chaperones." *MSB.* 2010 v.6 p.414

.. [#Reynolds2013] Reynolds KA, Russ WP, Socolich M, Ranganathan R.
   "Evolution-based design of proteins." *Methods Enzymol.* 2013 v.523 p.213

.. [#Rivoire2016] Rivoire, O., Reynolds, K. A., and Ranganathan, R.
   Evolution-Based Functional Decomposition of Proteins. *PLOS Computational
   Biology* 12, e1004817 (2016).
