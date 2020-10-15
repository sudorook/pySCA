#! /bin/bash
set -eu

# Globals

datadir=data
outputdir=output

datarepo="https://github.com/sudorook/pySCA-data"
version=6.2

# Download the data

cd ../

# In the event git is not installed, just directly download the data from
# GitHub using wget or curl (in order of preference). Also, check to see if tar
# is installed. If not, download the zipped archive.
if [ -x "$(command -v git)" ] && [ -d ".git/" ]; then
  git submodule init
  git submodule update --force
elif [ -x "$(command -v wget)" ]; then
  echo "git not installed --- trying wget"
  if [ -x "$(command -v tar)" ]; then
    wget -nc ${datarepo}/archive/v${version}.tar.gz
    tar xf v${version}.tar.gz
  elif [ -x "$(command -v unzip)" ]; then
    wget -nc ${datarepo}/archive/v${version}.zip
    unzip v${version}.zip
  else
    echo "'unzip' or 'tar' (with gzip) is required for decompressing data."
    exit 3
  fi
  mkdir -p ${datadir}
  mv -v pySCA-data-${version}/* ${datadir}/
  rm -rvf pySCA-data-${version}
elif [ -x "$(command -v curl)" ]; then
  echo "git not installed --- trying curl"
  if [ -x "$(command -v tar)" ]; then
    curl -L -O -C - ${datarepo}/archive/v${version}.tar.gz
    tar xf v${version}.tar.gz
  elif [ -x "$(command -v unzip)" ]; then
    curl -L -O -C - ${datarepo}/archive/v${version}.zip
    unzip v${version}.zip
  else
    echo "'unzip' or 'tar' (with gzip) is required for decompressing data."
    exit 3
  fi
  mkdir -p ${datadir}
  mv -v pySCA-data-${version}/* ${datadir}/
  rm -rvf pySCA-data-${version}
fi

# Generate the output files

mkdir -vp ${outputdir}

# The S1A serine proteases
echo "S1A serine protease Calculations:" | tee ${outputdir}/s1A_halabi.log
scaProcessMSA \
  -a ${datadir}/s1Ahalabi_1470_nosnakes.an \
  -b ${datadir} \
  -s 3TGI \
  -c E \
  -d ${outputdir} \
  -t -n 2>&1 | tee -a ${outputdir}/s1A_halabi.log
scaCore -i ${outputdir}/s1Ahalabi_1470_nosnakes.db 2>&1 | \
  tee -a ${outputdir}/s1A_halabi.log
scaSectorID -i ${outputdir}/s1Ahalabi_1470_nosnakes.db 2>&1 | \
  tee -a ${outputdir}/s1A_halabi.log
echo

# Beta-lactamase
echo "Beta-lactamase Calculations:" | tee ${outputdir}/PF13354.log
scaProcessMSA \
  -a ${datadir}/PF13354_full.an \
  -b ${datadir} \
  -s 1FQG \
  -c A \
  -d ${outputdir} \
  -f 'Escherichia coli' \
  -t -n 2>&1 | tee -a ${outputdir}/PF13354.log
scaCore -i ${outputdir}/PF13354_full.db 2>&1 | \
  tee -a ${outputdir}/PF13354.log
scaSectorID -i ${outputdir}/PF13354_full.db 2>&1 | \
  tee -a ${outputdir}/PF13354.log
echo

# G-protein - this analysis is run with two alignments - the full Pfam
# alignment (PF00071_full) and the Pfam alignment filtered to remove several
# N-terminal truncation mutants. PF00071_rd2 is the aligment discussed in the
# manuscript.
echo "G-protein calculations:" | tee ${outputdir}/PF00071.log
scaProcessMSA \
  -a ${datadir}/PF00071_full.an \
  -b ${datadir} \
  -s 5P21 \
  -c A \
  -d ${outputdir} \
  -f 'Homo sapiens' \
  -t -n 2>&1 | tee -a ${outputdir}/PF00071.log
scaCore -i ${outputdir}/PF00071_full.db 2>&1 | \
  tee -a ${outputdir}/PF00071.log
scaSectorID -i ${outputdir}/PF00071_full.db 2>&1 | \
  tee -a ${outputdir}/PF00071.log
echo

echo "G-protein calculations:" | tee ${outputdir}/PF00071_rd2.log
scaProcessMSA \
  -a ${datadir}/PF00071_rd2.an \
  -b ${datadir} \
  -s 5P21 \
  -c A \
  -d ${outputdir} \
  -f 'Homo sapiens' \
  -t -n 2>&1 | tee -a ${outputdir}/PF00071_rd2.log
scaCore -i ${outputdir}/PF00071_rd2.db 2>&1 | \
  tee -a ${outputdir}/PF00071_rd2.log
scaSectorID -i ${outputdir}/PF00071_rd2.db 2>&1 | \
  tee -a ${outputdir}/PF00071_rd2.log
echo

# DHFR - this analysis is also run with two alignments for comparison -
# the full PFAM alignment (PF00186_full.an) and a manually curated alignment
# (DHFR_PEPM3.an)
echo "DHFR Calculations:" | tee ${outputdir}/PF00186.log
scaProcessMSA \
  -a ${datadir}/PF00186_full.an \
  -b ${datadir} \
  -s 1RX2 \
  -c A \
  -d ${outputdir} \
  -f 'Escherichia coli' \
  -t -n 2>&1 | tee -a ${outputdir}/PF00186.log
scaCore -i ${outputdir}/PF00186_full.db 2>&1 | \
  tee -a ${outputdir}/PF00186.log
scaSectorID -i ${outputdir}/PF00186_full.db 2>&1 | \
  tee -a ${outputdir}/PF00186.log
echo

echo "DHFR Calculations:" | tee ${outputdir}/DHFR_PEPM3.log
scaProcessMSA \
  -a ${datadir}/DHFR_PEPM3.an \
  -b ${datadir} \
  -s 1RX2 \
  -c A \
  -d ${outputdir} \
  -t -n 2>&1 | tee -a ${outputdir}/DHFR_PEPM3.log
scaCore -i ${outputdir}/DHFR_PEPM3.db 2>&1 | \
  tee -a ${outputdir}/DHFR_PEPM3.log
scaSectorID -i ${outputdir}/DHFR_PEPM3.db 2>&1 | \
  tee -a ${outputdir}/DHFR_PEPM3.log
