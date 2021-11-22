#! /bin/bash
set -eu

#
# This script is intended to change the URLs and file names of the zipped
# output from a Jupyter notebook (see 'Download as rst' option). By default,
# the images names are 'output_\*.png', with corresponding URLs in the RST
# file, which will cause naming clashes when including multiple notebooks, each
# with different images, to the _static folder.
#
# This script will take the FILENAME of the zip file, extract its contents,
# rename the images from 'output' to '<FILENAME>', and update the URLS in the
# RST files.
#
# Input:
#   - zip file containing RST and image from Jupyter notebook
# Output:
#   - directory containing RST file and images with updated URLs
#
# Usage:
#   ./rstZipFixUrl.sh <path to zip>
#

DOCS_DIR="../docs/source"
DOCS_STATIC_DIR="_static"

FILENAME=$(basename "${1%.*}")
EXTENSION=${1##*.}

if [[ "${EXTENSION}" != "zip" ]]; then
  echo "ERROR: Input is not a zipped archive."
  exit 3
fi

TMPDIR=tmp_${FILENAME}

mkdir -p "${TMPDIR}"
cd "${TMPDIR}"

unzip ../"${1}"

sed -i "s,output_\([0-9_]\+\).png,${DOCS_STATIC_DIR}/${FILENAME}_\1.png,g" "${FILENAME}.rst"
sed -i "s,^\.\. code:: ipython3,\.\. code:: python3,g" "${FILENAME}.rst"
for PNG in *.png; do
  mv "${PNG}" "${PNG//output/${FILENAME}}"
done

cd ../

mv "${TMPDIR}/${FILENAME}.rst" ${DOCS_DIR}/
mv "${TMPDIR}/${FILENAME}"_*.png ${DOCS_DIR}/${DOCS_STATIC_DIR}/

rmdir "${TMPDIR}"
