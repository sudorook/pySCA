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

PFAMURL="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files"
PFAMHEADERS="pfamseq.sql"
PFAMDATA="pfamseq.txt"
PFAMMD5="md5_checksums"
PFAMDB="pfamseq.db"

GZIP=pigz  # replace this value with whatever GZIP compression tool you use


#
# Download and extract the data
#

echo "Downloading the Pfam database files and generate a SQLite3 database."
echo "Requires ~90 GB of free storage and could take several hours."

echo "Downloading the Pfam annotated sequence data."

wget -Nc "${PFAMURL}/${PFAMHEADERS}.gz"
wget -Nc "${PFAMURL}/${PFAMDATA}.gz"
wget -Nc "${PFAMURL}/${PFAMMD5}"
! grep "\s\+${PFAMHEADERS}\|\s\+${PFAMDATA}" "${PFAMMD5}" | md5sum -c - && exit 3
echo "Got 'em."

echo "Decompress the gzipped files."
echo "This will take a while."
if test "$(command -v ${GZIP})"; then
  ${GZIP} -vd "${PFAMHEADERS}.gz"
  ${GZIP} -vd "${PFAMDATA}.gz"
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
./mysql2sqlite "${PFAMHEADERS}" | sqlite3 "${PFAMDB}"

echo "Importing data."
sqlite3 -batch "${PFAMDB}" << EOF
.separator "\t"
.import pfamseq.txt pfamseq
EOF
echo "Done!"
