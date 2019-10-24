#! /usr/bin/env python3
""" Import text files into databases. """

import sqlite3
from Bio import SeqIO


def import_accession2taxid(file_in, db_out, table_name):
    """ Import the dead_prot.accession2taxid file. """
    fields = ["[accession]", "[accession.version]", "[taxid]", "[gi]"]
    fields = ",".join(fields)

    conn = sqlite3.connect(db_out)
    cur = conn.cursor()

    # initialize the database
    cmd = "PRAGMA synchronous = OFF;\n"
    cmd += "PRAGMA journal_mode = MEMORY;\n"
    cmd += "BEGIN TRANSACTION;\n"
    cmd += "DROP TABLE IF EXISTS `%s`;\n" % (table_name)
    cmd += "CREATE TABLE `%s` (\n" % (table_name)
    cmd += "  `accession` varchar(16) NOT NULL\n"
    cmd += ",  `accession.version` varchar(16) NOT NULL\n"
    cmd += ",  `taxid` varchar(16) NOT NULL\n"
    cmd += ",  `gi` varchar(16) NOT NULL\n"
    cmd += ",  PRIMARY KEY (`accession`)\n"
    cmd += ",  UNIQUE (`accession`)\n"
    cmd += ",  UNIQUE (`accession.version`)\n"
    cmd += ",  UNIQUE (`gi`)\n"
    cmd += ");\n"
    cmd += "END TRANSACTION;"
    cur.executescript(cmd)

    with open(file_in, "r") as handle:
        next(handle)  # skip header in first line
        for line in handle.readlines():
            row = tuple([field.strip() for field in line.split("\t")])
            cmd = "INSERT INTO %s(" % (table_name) + fields + ")\n"
            cmd += "VALUES(?,?,?,?);"
            cur.execute(cmd, row)
    conn.commit()
    conn.close()


def import_taxonomy(file_in, db_out, table_name):
    """ Import the taxid lineage file. """
    fields = ",".join(["[taxid]", "[species]", "[lineage]"])

    conn = sqlite3.connect(db_out)
    cur = conn.cursor()

    # initialize the database
    cmd = "PRAGMA synchronous = OFF;\n"
    cmd += "PRAGMA journal_mode = MEMORY;\n"
    cmd += "BEGIN TRANSACTION;\n"
    cmd += "DROP TABLE IF EXISTS `%s`;\n" % (table_name)
    cmd += "CREATE TABLE `%s` (\n" % (table_name)
    cmd += "  `taxid` varchar(16) NOT NULL\n"
    cmd += ",  `species` text NOT NULL\n"
    cmd += ",  `lineage` text NOT NULL\n"
    cmd += ",  PRIMARY KEY (`taxid`)\n"
    cmd += ",  UNIQUE (`taxid`)\n"
    cmd += ");\n"
    cmd += "END TRANSACTION;"
    cur.executescript(cmd)

    with open(file_in, "r") as handle:
        for line in handle.readlines():
            row = tuple([field.strip() for field in line.split("|")][:-1])
            cmd = "INSERT INTO %s(" % (table_name) + fields + ")\n"
            cmd += "VALUES(?,?,?);"
            cur.execute(cmd, row)
    conn.commit()
    conn.close()


def import_nr(file_in, db_out, table_name, fts=False):
    """ Import the nr (non-redundant) protein sequences. """
    fields = ",".join(["[accession.version]", "[sequence]"])

    conn = sqlite3.connect(db_out)
    cur = conn.cursor()

    # initialize the database
    cmd = "PRAGMA synchronous = OFF;\n"
    cmd += "PRAGMA journal_mode = MEMORY;\n"
    if fts:
        cmd += "BEGIN TRANSACTION;\n"
        cmd += "DROP TABLE IF EXISTS `%s`;\n" % (table_name)
        cmd += "CREATE VIRTUAL TABLE `%s` using fts4(\n" % (table_name)
        cmd += "  `accession.version` varchar(16) NOT NULL\n"
        cmd += ",  `sequence` text NOT NULL\n"
        cmd += ");\n"
        cmd += "END TRANSACTION;"
    else:
        cmd += "BEGIN TRANSACTION;\n"
        cmd += "DROP TABLE IF EXISTS `%s`;\n" % (table_name)
        cmd += "CREATE TABLE `%s` (\n" % (table_name)
        cmd += "  `accession.version` varchar(16) NOT NULL\n"
        cmd += ",  `sequence` text NOT NULL\n"
        cmd += ",  PRIMARY KEY (`accession.version`)\n"
        cmd += ",  UNIQUE (`accession.version`)\n"
        cmd += ");\n"
        cmd += "END TRANSACTION;"
    cur.executescript(cmd)

    with open(file_in, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            row = tuple([record.id, str(record.seq)])
            cmd = "INSERT INTO %s(" % (table_name) + fields + ")\n"
            cmd += "VALUES(?,?);"
            cur.execute(cmd, row)
    conn.commit()
    conn.close()


def main():
    """ main """
    dead_prot_file = "raw_data/dead_prot.accession2taxid"
    dead_prot_db = "dead_prot_accession2taxid.db"
    dead_prot_name = "dead_prot"
    prot_file = "raw_data/prot.accession2taxid"
    prot_db = "prot_accession2taxid.db"
    prot_name = "prot"

    print("Importing dead accession to taxid databases.")
    import_accession2taxid(dead_prot_file, dead_prot_db, dead_prot_name)

    print("Importing accession to taxid databases.")
    import_accession2taxid(prot_file, prot_db, prot_name)

    taxonomy_file = "raw_data/fullnamelineage.dmp"
    taxonomy_db = "taxonomy.db"
    taxonomy_name = "taxonomy"

    print("Importing taxonomy fullnamelineage database.")
    import_taxonomy(taxonomy_file, taxonomy_db, taxonomy_name)

    nr_prot_file = "raw_data/nr"
    nr_prot_db = "raw_data/nr.db"
    nr_prot_name = "nr"
    nr_prot_fts_db = "raw_data/nr_fts.db"
    nr_prot_fts_name = "nr_fts"

    print("Importing non-redundant (nr) protein database.")
    import_nr(nr_prot_file, nr_prot_db, nr_prot_name, fts=False)

    print("Importing non-redundant (nr) database with full-text search.")
    import_nr(nr_prot_file, nr_prot_fts_db, nr_prot_fts_name, fts=True)

    print("Done!")


main()
