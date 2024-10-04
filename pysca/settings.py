#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2019 - 2020 sudorook <daemon@nullcodon.com>
#
# SPDX-License-Identifier: BSD-3-Clause

""" Global settings for pySCA. """

#
# PATHS
#
# These have to be changed to be consistent with user-defined paths. This
# script is tested against the `runAllNBCalcs.sh` scripts, and because the
# script includes a `cd ../` command before running any of the python scripts,
# the base directory is the root of the repository.
#

# Enter absolute path (e.g. /home/<user>/pfamseq.txt) to the file 'pfamseq.txt'
# from
# ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/database_files/
# and/or the SQLite database `pfamseq.db` if it exists.
path2pfamseq = "pfamseq.txt"  # replace with absolute path to pfamseq.txt
path2pfamseqdb = (
    "pfamseq.db"  # replace with absolute path to pfamseq.db (if present)
)

# the location of your PDB structures
path2structures = (
    "."  # replace with absolute path to directory of PDB structures
)

# Also assumes that a folder named 'output/' is in the path. Change to '.' if
# you want results printed in the current working directory by default.
path2output = "output/"

# Used for pulling species, taxonomy annotations from ncbi database. PLEASE
# change to your own your email!!!
entrezemail = "your.email@youruniversity.edu"

# If you are using a version of PyMOL not intalled in your system PATH, you can
# add the path here. Use an absolute path to the PyMOL executable, or leave
# empty to use PyMOL on the system PATH.
path2pymol = ""
