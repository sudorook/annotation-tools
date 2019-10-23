#! /bin/bash
set -eu

#
# Script for downloading all the NCBI files needed for annotation of raw
# sequences.
#
# Requires wget and a gzip decompression tool. I strongly recommend using pigz
# instead of gzip. It's a parallel gzip implementation and saves a lot of time
# when decompressing large files.
#

#
# Globals
#
database_dir=`pwd`
gzip=pigz

# Check that pigz is installed before doing anything.
if ! test $(command -v ${gzip}); then
  echo "${gzip} not installed."
  exit 3
fi

# Simple wrapper function that will check the MD5 sum for a file against the
# sum downloaded from the server. Probably redundant with `wget -nc`.
function check_md5sum {
  local file="${1}"

  # Return failure if md5 file not found.
  if ! [ -f "${file}.md5" ]; then
    return 1
  fi
  
  # Check the md5sum against the md5 file.
  if `cat "${file}.md5" | diff - <(md5sum "${file}") >/dev/null`; then
    return 0
  else
    return 1
  fi
}


#
# Download: nr (non-redundant) proteins
#
nr_url="ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA"
nr_archive="nr.gz"
wget -nc ${nr_url}/${nr_archive}
wget -nc ${nr_url}/${nr_archive}.md5

# Check data integrity before decompressing.
if `check_md4sum "${nr_archive}"`; then
  echo "MD5 success. Extracting ${nr_archive}."
  ${gzip} -d ${nr_archive}
  rm ${nr_archive}.md5
  chmod -w,o-r nr
else
  echo "MD5 mismatch. Exiting."
  exit 1
fi


#
# Download: Accession to TaxID
#
accession2taxid_url="ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid"
dead_prot_archive="dead_prot.accession2taxid.gz"
prot_archive="prot.accession2taxid.gz"
wget -nc ${accession2taxid_url}/${dead_prot_archive}
wget -nc ${accession2taxid_url}/${dead_prot_archive}.md5
wget -nc ${accession2taxid_url}/${prot_archive}
wget -nc ${accession2taxid_url}/${prot_archive}.md5

# Check data integrity for prot database before decompressing.
if `check_md5sum "${dead_prot_archive}"`; then
  echo "MD5 success. Extracting ${dead_prot_archive}."
  ${gzip} -d ${dead_prot_archive}
  rm ${dead_prot_archive}.md5
  chmod -w,o-r dead_prot.accession2taxid
else
  echo "MD5 mismatch. Exiting."
  exit 3
fi

# Check data integrity for prot database before decompressing.
if `check_md5sum "${prot_archive}"`; then
  echo "MD5 success. Extracting ${prot_archive}."
  ${gzip} -d ${prot_archive}
  rm ${prot_archive}.md5
  chmod -w,o-r prot.accession2taxid
else
  echo "MD5 mismatch. Exiting."
  exit 3
fi

#
# Download: Taxonomy
#
taxonomy_url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump"
taxonomy_archive="new_taxdump.tar.gz"
wget -nc ${taxonomy_url}/${taxonomy_archive}
wget -nc ${taxonomy_url}/${taxonomy_archive}.md5

# Check data integrity before decompressing.
if `check_md5sum "${taxonomy_archive}"`; then
  echo "MD5 success. Extracting ${taxonomy_archive}."
  ${gzip} -dc ${taxonomy_archive} | tar xf - fullnamelineage.dmp
  rm ${taxonomy_archive}
  rm ${taxonomy_archive}.md5
  chmod -w,o-r fullnamelineage.dmp
else
  echo "MD5 mismatch. Exiting."
  exit 3
fi
