# Annotation Tools

Scripts for downloading annotation data and generating databases used for
speedy annotation of raw sequence data.

Don't run these scripts unless you have enough storage capacity (~700 GB, to be
safe).

1. `download_data.sh` - download the pfamseq, nr, accessiont2id, and taxonomy
   databases. Checks the MD5 sums.
2. `extract_data.sh` - extract the downloaded data.
3. `import_sqlite.py` - generate sqlite3 databases from the extracted data.
   This is very slow and memory-hungry.
4. `set_permissions.sh` - set the permissions for a set of files to be
   non-writable. Prevents accidental writing of changes or deletions of files.
