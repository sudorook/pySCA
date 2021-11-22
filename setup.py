#! /usr/bin/env python3
""" Installation and setup for pySCA"""

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
        "numpy",
        "scipy",
        "argparse",
        "wheel",
        "matplotlib",
    ],
    scripts=[
        "bin/alnChangeDelim",
        "bin/alnFilterSeqSize",
        "bin/alnParseID",
        "bin/annotateMSA",
        "bin/scaProcessMSA",
        "bin/alnConvertGI",
        "bin/alnReplaceHeaders",
        "bin/scaCore",
        "bin/scaSectorID",
    ],
)
