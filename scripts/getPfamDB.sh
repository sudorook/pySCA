#! /bin/bash
set -eu

#
# The Pfam annotation script is much, much faster when using a database instead
# of iterating over a 20 GB text file line by line. This script is intended to
# download the text file and convert it into a SQLite3 database.
#
# I recommend running this overnight when you aren't using your computer.
# SQLite3 has to create key-value pairs for over 40 million sequences, and it
# is VERY, VERY slow.
#
# Dependencies: wget, sqlite3, awk, and gzip or pigz
#

#
# Globals
#

pfamurl="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files"
pfamheaders="pfamseq.sql"
pfamdata="pfamseq.txt"
pfamdb="pfamseq.db"

gzip=gzip  # replace this value with whatever gzip compression tool you use


#
# Download and extract the data
#

echo "Downloading the Pfam database files and generate a SQLite3 database."
echo "Requires ~90 GB of free storage and could take several hours."

echo "Downloading the Pfam annotated sequence data."

wget -Nc "${pfamurl}/${pfamheaders}.gz"
wget -Nc "${pfamurl}/${pfamdata}.gz"
echo "Got 'em."

echo "Decompress the gzipped files."
echo "This will take a while."
if test "$(command -v ${gzip})"; then
  ${gzip} -vd "${pfamheaders}.gz"
  ${gzip} -vd "${pfamdata}.gz"
else
  echo "${gzip} not found. Exiting."
  exit 3
fi
echo "Done!"


#
# Create the database
#

# The SQL dump on the server is for MySQL (MariaDB), so it needs to be
# converted to a format compatible with SQLite3.

echo "Converting the MySQL dump to SQLite3."
./mysql2sqlite "${pfamheaders}" | sqlite3 "${pfamdb}"

echo "Importing data."
sqlite3 -batch "${pfamdb}" << "EOF"
.separator "\t"
.import pfamseq.txt pfamseq
EOF
echo "Done!"
