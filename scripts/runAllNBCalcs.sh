#! /bin/bash
set -eu

# Globals

DATA_DIR=data
OUTPUT_DIR=output

DATA_REPO="https://github.com/sudorook/pySCA-data"
VERSION=6.2

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
    wget -nc ${DATA_REPO}/archive/v${VERSION}.tar.gz
    tar xf v${VERSION}.tar.gz
  elif [ -x "$(command -v unzip)" ]; then
    wget -nc ${DATA_REPO}/archive/v${VERSION}.zip
    unzip v${VERSION}.zip
  else
    echo "'unzip' or 'tar' (with gzip) is required for decompressing data."
    exit 3
  fi
  mkdir -p ${DATA_DIR}
  mv -v pySCA-data-${VERSION}/* ${DATA_DIR}/
  rm -rvf pySCA-data-${VERSION}
elif [ -x "$(command -v curl)" ]; then
  echo "git not installed --- trying curl"
  if [ -x "$(command -v tar)" ]; then
    curl -L -O -C - ${DATA_REPO}/archive/v${VERSION}.tar.gz
    tar xf v${VERSION}.tar.gz
  elif [ -x "$(command -v unzip)" ]; then
    curl -L -O -C - ${DATA_REPO}/archive/v${VERSION}.zip
    unzip v${VERSION}.zip
  else
    echo "'unzip' or 'tar' (with gzip) is required for decompressing data."
    exit 3
  fi
  mkdir -p ${DATA_DIR}
  mv -v pySCA-data-${VERSION}/* ${DATA_DIR}/
  rm -rvf pySCA-data-${VERSION}
fi

# Generate the output files

mkdir -vp ${OUTPUT_DIR}

# The S1A serine proteases
echo "S1A serine protease Calculations:" | tee ${OUTPUT_DIR}/s1A_halabi.log
scaProcessMSA \
  -a ${DATA_DIR}/s1Ahalabi_1470_nosnakes.an \
  -b ${DATA_DIR} \
  -s 3TGI \
  -c E \
  -d ${OUTPUT_DIR} \
  -t -n 2>&1 | tee -a ${OUTPUT_DIR}/s1A_halabi.log
scaCore -i ${OUTPUT_DIR}/s1Ahalabi_1470_nosnakes.db 2>&1 |
  tee -a ${OUTPUT_DIR}/s1A_halabi.log
scaSectorID -i ${OUTPUT_DIR}/s1Ahalabi_1470_nosnakes.db 2>&1 |
  tee -a ${OUTPUT_DIR}/s1A_halabi.log
echo

# Beta-lactamase
echo "Beta-lactamase Calculations:" | tee ${OUTPUT_DIR}/PF13354.log
scaProcessMSA \
  -a ${DATA_DIR}/PF13354_full.an \
  -b ${DATA_DIR} \
  -s 1FQG \
  -c A \
  -d ${OUTPUT_DIR} \
  -f 'Escherichia coli' \
  -t -n 2>&1 | tee -a ${OUTPUT_DIR}/PF13354.log
scaCore -i ${OUTPUT_DIR}/PF13354_full.db 2>&1 |
  tee -a ${OUTPUT_DIR}/PF13354.log
scaSectorID -i ${OUTPUT_DIR}/PF13354_full.db 2>&1 |
  tee -a ${OUTPUT_DIR}/PF13354.log
echo

# G-protein - this analysis is run with two alignments - the full Pfam
# alignment (PF00071_full) and the Pfam alignment filtered to remove several
# N-terminal truncation mutants. PF00071_rd2 is the aligment discussed in the
# manuscript.
echo "G-protein calculations:" | tee ${OUTPUT_DIR}/PF00071.log
scaProcessMSA \
  -a ${DATA_DIR}/PF00071_full.an \
  -b ${DATA_DIR} \
  -s 5P21 \
  -c A \
  -d ${OUTPUT_DIR} \
  -f 'Homo sapiens' \
  -t -n 2>&1 | tee -a ${OUTPUT_DIR}/PF00071.log
scaCore -i ${OUTPUT_DIR}/PF00071_full.db 2>&1 |
  tee -a ${OUTPUT_DIR}/PF00071.log
scaSectorID -i ${OUTPUT_DIR}/PF00071_full.db 2>&1 |
  tee -a ${OUTPUT_DIR}/PF00071.log
echo

echo "G-protein calculations:" | tee ${OUTPUT_DIR}/PF00071_rd2.log
scaProcessMSA \
  -a ${DATA_DIR}/PF00071_rd2.an \
  -b ${DATA_DIR} \
  -s 5P21 \
  -c A \
  -d ${OUTPUT_DIR} \
  -f 'Homo sapiens' \
  -t -n 2>&1 | tee -a ${OUTPUT_DIR}/PF00071_rd2.log
scaCore -i ${OUTPUT_DIR}/PF00071_rd2.db 2>&1 |
  tee -a ${OUTPUT_DIR}/PF00071_rd2.log
scaSectorID -i ${OUTPUT_DIR}/PF00071_rd2.db 2>&1 |
  tee -a ${OUTPUT_DIR}/PF00071_rd2.log
echo

# DHFR - this analysis is also run with two alignments for comparison -
# the full PFAM alignment (PF00186_full.an) and a manually curated alignment
# (DHFR_PEPM3.an)
echo "DHFR Calculations:" | tee ${OUTPUT_DIR}/PF00186.log
scaProcessMSA \
  -a ${DATA_DIR}/PF00186_full.an \
  -b ${DATA_DIR} \
  -s 1RX2 \
  -c A \
  -d ${OUTPUT_DIR} \
  -f 'Escherichia coli' \
  -t -n 2>&1 | tee -a ${OUTPUT_DIR}/PF00186.log
scaCore -i ${OUTPUT_DIR}/PF00186_full.db 2>&1 |
  tee -a ${OUTPUT_DIR}/PF00186.log
scaSectorID -i ${OUTPUT_DIR}/PF00186_full.db 2>&1 |
  tee -a ${OUTPUT_DIR}/PF00186.log
echo

echo "DHFR Calculations:" | tee ${OUTPUT_DIR}/DHFR_PEPM3.log
scaProcessMSA \
  -a ${DATA_DIR}/DHFR_PEPM3.an \
  -b ${DATA_DIR} \
  -s 1RX2 \
  -c A \
  -d ${OUTPUT_DIR} \
  -t -n 2>&1 | tee -a ${OUTPUT_DIR}/DHFR_PEPM3.log
scaCore -i ${OUTPUT_DIR}/DHFR_PEPM3.db 2>&1 |
  tee -a ${OUTPUT_DIR}/DHFR_PEPM3.log
scaSectorID -i ${OUTPUT_DIR}/DHFR_PEPM3.db 2>&1 |
  tee -a ${OUTPUT_DIR}/DHFR_PEPM3.log
