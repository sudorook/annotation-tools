#! /bin/bash
set -eu

#
# Script for downloading all the NCBI files needed for annotation of raw
# sequences.
#
# Requires wget.
#

#
# Globals
#
database_dir=raw_data

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

mkdir -p ${database_dir}
cd ${database_dir}


#
# Download: nr (non-redundant) proteins
#
nr_url="ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA"
nr_archive="nr.gz"
wget -nc ${nr_url}/${nr_archive}
wget -nc ${nr_url}/${nr_archive}.md5

# Check data integrity before decompressing.
if `check_md5sum "${nr_archive}"`; then
  echo "MD5 match."
  rm ${nr_archive}.md5
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
  echo "MD5 match."
  rm ${dead_prot_archive}.md5
else
  echo "MD5 mismatch. Exiting."
  exit 3
fi

# Check data integrity for prot database before decompressing.
if `check_md5sum "${prot_archive}"`; then
  echo "MD5 match."
  rm ${prot_archive}.md5
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
  echo "MD5 match."
  rm ${taxonomy_archive}.md5
else
  echo "MD5 mismatch. Exiting."
  exit 3
fi


#
# Download: Pfam
#
pfam_url="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files"
# pfam_sql="pfamseq.sql.gz"
pfam_archive="pfamseq.txt.gz"

# wget -nc "${pfam_url}/${pfam_sql}"
wget -nc "${pfam_url}/${pfam_archive}"
