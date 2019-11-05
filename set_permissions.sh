#! /bin/bash
set -eu

#
# Short script for disabling write permsissions for a set of files. Useful for
# preventing accidental overwriting large files that are time consuming to
# download or generate.
#

data_dirs=(raw_data data)

files=(nr dead_prot.accession2taxid prot.accession2taxid fullnamelineage.dmp
       nr.db taxonomy.db dead_prot_accession2taxid.db prot_accession2taxid.db
       pfamseq.db)

for dir in ${data_dirs[@]}; do
  if [ -d "${dir}" ]; then
    cd "${dir}"
    for file in ${files[@]}; do
      if [ -f "${file}" ]; then
        chmod -w,o-r "${file}"
      fi
    done
    cd ../
  fi
done
