#! /bin/bash
set -eu

#
# This script is intended to change the URLs and file names of the zipped
# output from a Jupyter notebook (see 'Download as rst' option). By default,
# the images names are 'output_\*.png', with corresponding URLs in the RST
# file, which will cause naming clashes when including multiple notebooks, each
# with different images, to the _static folder.
#
# This script will take the filename of the zip file, extract its contents,
# rename the images from 'output' to '<filename>', and update the URLS in the
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

docsdir="../docs/source"
docsstaticdir="_static"

filename=$(basename "${1%.*}")
extension=${1##*.}

if [[ "${extension}" != "zip" ]]; then
  echo "ERROR: Input is not a zipped archive."
  exit 3
fi

tmpdir=tmp_${filename}

mkdir -p "${tmpdir}"
cd "${tmpdir}"

unzip ../"${1}"

sed -i "s,output_\([0-9_]\+\).png,${docsstaticdir}/${filename}_\1.png,g" "${filename}.rst"
sed -i "s,^\.\. code:: ipython3,\.\. code:: python3,g" "${filename}.rst"
for png in *.png; do
  newpng=$(echo "${png}" | sed -e "s/output_\([0-9_]\+\).png/${filename}_\1.png/g")
  mv "${png}" "${newpng}"
done

cd ../

mv "${tmpdir}/${filename}.rst" ${docsdir}/
mv "${tmpdir}/${filename}"_*.png ${docsdir}/${docsstaticdir}/

rmdir "${tmpdir}"
