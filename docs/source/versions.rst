.. SPDX-FileCopyrightText: 2019 - 2020 sudorook <daemon@nullcodon.com>
..
.. SPDX-License-Identifier: BSD-3-Clause

=============
Distributions
=============

Previous versions of SCA were implemented as MATLAB toolboxes and contain
various accessory codes for data formatting, display, and analysis.

:SCA Toolbox 1.5: 
  The original SCA method as specified in Lockless and Ranganathan (4) with one
  modification that was used in all subsequent papers: the division of binomial
  probabilities by the mean probability of amino acids in the alignment is
  removed. This version is longer in active use.

:SCA Toolbox 2.5:
  The bootstrap-based approach for SCA. Position-specific conservation
  calculated as in Eq. (4) and correlations calculated as in Eq.  (9). Matrix
  reduction per Eq. (32).

:SCA Toolbox 3.0:
  The analytical calculation of correlations weighted by gradients of relative
  entropy. Position-specific conservation calculated as in Eq. (4) and
  correlations calculated as in Eq. (9)-(33). For non-binarized alignments,
  matrix reduction is per Eq. (32).

:SCA Toolbox 4.0:
  Analytical calculations as in v3.0, but now including sector identification
  methods as described in Ref. (2).

:SCA Toolbox 5.0:
  Calculation of positional and sequence correlations matrices by the alignment
  projection method as per Eq. (19) and Eq. (20), and calculation of the
  mapping between them Eq. (21). Includes methods for sector identification and
  exploring relationships between positional and sequence correlations. 

:SCA Toolbox 6.0:
  Calculation of first-order and second-order statistics for positional amino
  acid frequencies using sequences weighted by similarity in the multiple
  sequence alignment. 

:SCA Toolbox 6.1:
  Port from Python 2 to Python 3. Also includes updated annotation scripts,
  changes to the command-line interface, and the option to install analysis
  scripts as system-wide executables.

:SCA Toolbox 6.2:
  Bugfix release.
