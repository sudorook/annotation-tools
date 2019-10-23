#! /bin/bash
set -eu

#
# Script for extracting already downloaded NCBI files needed for annotation of
# raw sequences.
#
# Requires a gzip decompression tool. I strongly recommend using pigz instead
# of gzip. It's a parallel gzip implementation and saves a lot of time when
# decompressing large files.
#

# Globals
database_dir=raw_data
gzip=pigz

# Check that pigz is installed before doing anything.
if ! test $(command -v ${gzip}); then
  echo "${gzip} not installed."
  exit 3
fi

mkdir -p ${database_dir}
cd ${database_dir}

# Extract: nr (non-redundant) proteins
nr_archive="nr.gz"

${gzip} -d ${nr_archive}
chmod -w,o-r nr


# Extract: Accession to TaxID
dead_prot_archive="dead_prot.accession2taxid.gz"
prot_archive="prot.accession2taxid.gz"

${gzip} -d ${dead_prot_archive}
chmod -w,o-r dead_prot.accession2taxid

${gzip} -d ${prot_archive}
chmod -w,o-r prot.accession2taxid

# Extract: Taxonomy
taxonomy_archive="new_taxdump.tar.gz"

${gzip} -dc ${taxonomy_archive} | tar xf - fullnamelineage.dmp
rm ${taxonomy_archive}
chmod -w,o-r fullnamelineage.dmp
