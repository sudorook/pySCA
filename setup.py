#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2019 - 2022 sudorook <daemon@nullcodon.com>
#
# SPDX-License-Identifier: BSD-3-Clause

"""Installation and setup for pySCA"""

from setuptools import setup

setup(
    name="pySCA",
    version="6.2",
    author="Olivier Rivoire, Rama Ranganathan, and Kimberly Reynolds",
    maintainer="Ansel George",
    packages=["pysca"],
    package_data={"pysca": ["settings.py"]},
    description="Python 3 implementation of Statistical Coupling Analysis (SCA)",
    url="https://sudorook.gitlab.io/pySCA",
    download_url="https://github.com/sudorook/pySCA",
    long_description=open("README.md", "r", encoding="utf-8").read(),
    install_requires=[
        "biopython",
        "matplotlib",
        "numpy",
        "scipy",
        "wheel",
    ],
    scripts=[
        "bin/alnChangeDelim",
        "bin/alnConvertGI",
        "bin/alnFilterSeqSize",
        "bin/alnParseID",
        "bin/alnReplaceHeaders",
        "bin/annotateMSA",
        "bin/scaCore",
        "bin/scaProcessMSA",
        "bin/scaSectorID",
    ],
)
