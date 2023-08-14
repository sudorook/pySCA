#! /bin/bash
set -euo pipefail

#
# The Pfam annotation script is much, much faster when using a database instead
# of iterating over a 20 GB text file line by line. This script is intended to
# download the text file and convert it into a SQLite3 database.
#
# I recommend running this overnight when you aren't using your computer.
# SQLite3 has to create key-value pairs for over 40 million sequences, and it
# is VERY, VERY slow.
#
# Dependencies:
# - wget
# - sqlite3
# - awk
# - gzip or pigz
#

#
# Globals
#

PFAM_URL="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files"
PFAM_HEADERS="pfamseq.sql"
PFAM_DATA="pfamseq.txt"
PFAM_MD5="md5_checksums"
PFAM_DB="pfamseq.db"

GZIP=pigz  # replace this value with whatever GZIP compression tool you use

#
# Download and extract the data
#

echo "Downloading the Pfam database files and generate a SQLite3 database."
echo "Requires ~90 GB of free storage and could take several hours."

echo "Downloading the Pfam annotated sequence data."

wget -Nc "${PFAM_URL}/${PFAM_HEADERS}.gz"
wget -Nc "${PFAM_URL}/${PFAM_DATA}.gz"
wget -Nc "${PFAM_URL}/${PFAM_MD5}"
! grep "\s\+${PFAM_HEADERS}\|\s\+${PFAM_DATA}" "${PFAM_MD5}" | md5sum -c - && exit 3
echo "Got 'em."

echo "Decompress the gzipped files."
echo "This will take a while."
if test "$(command -v ${GZIP})"; then
  ${GZIP} -vd "${PFAM_HEADERS}.gz"
  ${GZIP} -vd "${PFAM_DATA}.gz"
else
  echo "${GZIP} not found. Exiting."
  exit 3
fi
echo "Done!"

#
# Create the database
#

# The SQL dump on the server is for MySQL (MariaDB), so it needs to be
# converted to a format compatible with SQLite3.

echo "Converting the MySQL dump to SQLite3."
./mysql2sqlite "${PFAM_HEADERS}" | sqlite3 "${PFAM_DB}"

echo "Importing data."
sqlite3 -batch "${PFAM_DB}" << EOF
.separator "\t"
.import pfamseq.txt pfamseq
EOF
echo "Done!"
