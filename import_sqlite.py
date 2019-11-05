#! /usr/bin/env python3
""" Import text files into databases. """

import os
import sqlite3
from Bio import SeqIO


def truncate(sequence):
    """ Do nothing. Just a placeholder. """
    string = str(sequence)
    return string.split()[0]


def expand(sequence):
    """ Permute the string """
    string = str(sequence)
    string_length = len(string)
    return " ".join([string[idx:string_length]
                     for idx in range(0, string_length)])


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


def import_nr(file_in, db_out, table_name):
    """ Import the nr (non-redundant) protein sequences. """
    fields = ",".join(["[accession.version]", "[sequence]"])
    table = table_name
    table_fts = table_name + "_fts"

    conn = sqlite3.connect(db_out)
    conn.create_function("compress", 1, truncate)
    conn.create_function("uncompress", 1, expand)
    cur = conn.cursor()

    # Initialize the database
    cmd = "PRAGMA synchronous = OFF;\n"
    cmd += "PRAGMA journal_mode = MEMORY;\n"
    cmd += "BEGIN TRANSACTION;\n"
    cmd += "DROP TABLE IF EXISTS `%s`;\n" % (table)
    cmd += "CREATE TABLE `%s` (\n" % (table)
    cmd += "  `accession.version` varchar(16) NOT NULL\n"
    cmd += ",  `sequence` text NOT NULL\n"
    cmd += ",  PRIMARY KEY (`accession.version`)\n"
    cmd += ",  UNIQUE (`accession.version`)\n"
    cmd += ");\n"
    cmd += "DROP TABLE IF EXISTS `%s`;\n" % (table_fts)
    cmd += "CREATE VIRTUAL TABLE `%s` using fts4(\n" % (table_fts)
    cmd += "  `accession.version` varchar(16) NOT NULL\n"
    cmd += ",  `sequence` text NOT NULL\n"
    cmd += ",  compress=compress\n"
    cmd += ",  uncompress=uncompress\n"
    cmd += ");\n"
    cmd += "END TRANSACTION;"
    cur.executescript(cmd)

    # Insert the data
    with open(file_in, "r") as handle:
        block_size = 10000
        cur.execute("BEGIN TRANSACTION;")
        for i, record in enumerate(SeqIO.parse(handle, "fasta")):
            # Add entry to regular database
            row = tuple([record.id, str(record.seq)])
            cmd = "INSERT INTO %s(" % (table) + fields + ")\n"
            cmd += "VALUES(?,?);"
            cur.execute(cmd, row)

            # Add entry to FTS database
            row = tuple([record.id,
                         expand(str(record.seq))])
            cmd = "INSERT INTO %s(" % (table_fts) + fields + ")\n"
            cmd += "VALUES(?,?);"
            cur.execute(cmd, row)

            # Commit changes to databases after block_size records
            if not i % block_size:
                cur.execute("COMMIT;")
                cur.execute("BEGIN TRANSACTION;")
        cur.execute("END TRANSACTION;")
    conn.commit()
    conn.close()


def import_pfam(file_in, db_out, table_name):
    """ Import the Pfam protein sequences. """
    fields = ",".join(
        [
            "[pfamseq_acc]",
            "[pfamseq_id]",
            "[seq_version]",
            "[crc64]",
            "[md5]",
            "[description]",
            "[evidence]",
            "[length]",
            "[species]",
            "[taxonomy]",
            "[is_fragment]",
            "[sequence]",
            "[updated]",
            "[created]",
            "[ncbi_taxid]",
            "[auto_architecture]",
            "[treefam_acc]",
            "[swissprot]",
        ]
    )
    fields_fts = ",".join(
        [
            "[pfamseq_acc]",
            "[sequence]",
        ]
    )
    table = table_name
    table_fts = table_name + "_fts"

    conn = sqlite3.connect(db_out)
    conn.create_function("compress", 1, truncate)
    conn.create_function("uncompress", 1, expand)
    cur = conn.cursor()

    # Initialize the database
    cmd = "PRAGMA synchronous = OFF;\n"
    cmd += "PRAGMA journal_mode = MEMORY;\n"
    cmd += "BEGIN TRANSACTION;\n"
    cmd += "DROP TABLE IF EXISTS `%s`;\n" % (table)
    cmd += "DROP TABLE IF EXISTS `%s`;\n" % (table_fts)
    cmd += "CREATE VIRTUAL TABLE `%s` using fts4(\n" % (table_fts)
    cmd += "  `pfamseq_acc` varchar(10) NOT NULL\n"
    cmd += ",  `sequence` blob NOT NULL\n"
    cmd += ",  compress=compress\n"
    cmd += ",  uncompress=uncompress\n"
    cmd += ");\n"
    cmd += "CREATE TABLE `%s` (\n" % (table)
    cmd += "  `pfamseq_acc` varchar(10) NOT NULL\n"
    cmd += ",  `pfamseq_id` varchar(16) NOT NULL\n"
    cmd += ",  `seq_version` integer NOT NULL\n"
    cmd += ",  `crc64` varchar(16) NOT NULL\n"
    cmd += ",  `md5` varchar(32) NOT NULL\n"
    cmd += ",  `description` text NOT NULL\n"
    cmd += ",  `evidence` integer NOT NULL\n"
    cmd += ",  `length` integer NOT NULL DEFAULT '0'\n"
    cmd += ",  `species` text NOT NULL\n"
    cmd += ",  `taxonomy` mediumtext\n"
    cmd += ",  `is_fragment` integer DEFAULT NULL\n"
    cmd += ",  `sequence` blob NOT NULL\n"
    cmd += ",  `updated` timestamp NOT NULL DEFAULT current_timestamp\n"
    cmd += ",  `created` datetime DEFAULT NULL\n"
    cmd += ",  `ncbi_taxid` integer  NOT NULL DEFAULT '0'\n"
    cmd += ",  `auto_architecture` integer  DEFAULT NULL\n"
    cmd += ",  `treefam_acc` varchar(8) DEFAULT NULL\n"
    cmd += ",  `swissprot` integer DEFAULT '0'\n"
    cmd += ",  PRIMARY KEY (`pfamseq_acc`)\n"
    cmd += ",  UNIQUE (`pfamseq_acc`)\n"
    cmd += ",  CONSTRAINT `FK_pfamseq_1` FOREIGN KEY (`ncbi_taxid`) REFERENCES `ncbi_taxonomy` (`ncbi_taxid`) ON DELETE CASCADE ON UPDATE NO ACTION\n"
    cmd += ",  CONSTRAINT `fk_pfamseq_evidence1` FOREIGN KEY (`evidence`) REFERENCES `evidence` (`evidence`) ON DELETE NO ACTION ON UPDATE NO ACTION\n"
    cmd += ");\n"
    cmd += 'CREATE INDEX "idx_pfamseq_ncbi_taxid" ON "%s" (`ncbi_taxid`);\n' % (table)
    cmd += 'CREATE INDEX "idx_pfamseq_crc64" ON "%s" (`crc64`);\n' % (table)
    cmd += 'CREATE INDEX "idx_pfamseq_pfamseq_id" ON "%s" (`pfamseq_id`);\n' % (table)
    cmd += 'CREATE INDEX "idx_pfamseq_pfamseq_architecture_idx" ON "%s" (`auto_architecture`);\n' % (table)
    cmd += 'CREATE INDEX "idx_pfamseq_pfamseq_acc_version" ON "%s" (`pfamseq_acc`,`seq_version`);\n' % (table)
    cmd += 'CREATE INDEX "idx_pfamseq_md5" ON "%s" (`md5`);\n' % (table)
    cmd += 'CREATE INDEX "idx_pfamseq_pfamseq_tax_idx" ON "%s" (`taxonomy`);\n' % (table)
    cmd += 'CREATE INDEX "idx_pfamseq_fk_pfamseq_evidence1_idx" ON "%s" (`evidence`);\n' % (table)
    cmd += 'CREATE INDEX "idx_pfamseq_fragment_idx" ON "%s" (`is_fragment`);\n' % (table)
    cmd += 'CREATE INDEX "idx_pfamseq_swissprot_idx" ON "%s" (`swissprot`);\n' % (table)
    cmd += "END TRANSACTION;"
    cur.executescript(cmd)

    with open(file_in, "r") as handle:
        block_size = 10000
        cur.execute("BEGIN TRANSACTION;")
        for i, line in enumerate(handle.readlines()):
            # Add record to regular table.
            row = tuple([field.strip() for field in line.split("\t")])
            cmd = "INSERT INTO %s(" % (table) + fields + ")\n"
            cmd += "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);"
            cur.execute(cmd, row)

            # Add record to FTS table.
            row_data = [field.strip() for field in line.split("\t")]
            row = tuple([row_data[0], expand(row_data[11])])
            cmd = "INSERT INTO %s(" % (table_fts) + fields_fts + ")\n"
            cmd += "VALUES(?,?);"
            cur.execute(cmd, row)

            # Commit changes to databases after block_size records
            if not i % block_size:
                cur.execute("COMMIT;")
                cur.execute("BEGIN TRANSACTION;")
    conn.commit()
    conn.close()


def main():
    """ main """

    if not os.path.exists("data"):
        os.mkdir("data")

    # Accession to TaxID
    dead_prot_file = "raw_data/dead_prot.accession2taxid"
    dead_prot_db = "data/dead_prot_accession2taxid.db"
    dead_prot_name = "dead_prot"
    prot_file = "raw_data/prot.accession2taxid"
    prot_db = "data/prot_accession2taxid.db"
    prot_name = "prot"

    print("Importing dead accession to taxid databases.")
    import_accession2taxid(dead_prot_file, dead_prot_db, dead_prot_name)

    print("Importing accession to taxid databases.")
    import_accession2taxid(prot_file, prot_db, prot_name)

    # Taxonomy
    taxonomy_file = "raw_data/fullnamelineage.dmp"
    taxonomy_db = "data/taxonomy.db"
    taxonomy_name = "taxonomy"

    print("Importing taxonomy fullnamelineage database.")
    import_taxonomy(taxonomy_file, taxonomy_db, taxonomy_name)

    # Non-redundant (nr) proteins
    nr_prot_file = "raw_data/nr"
    nr_prot_db = "data/nr.db"
    nr_prot_name = "nr"

    print("Importing non-redundant (nr) database with full-text search.")
    import_nr(nr_prot_file, nr_prot_db, nr_prot_name)

    # Pfam
    pfam_file = "raw_data/pfamseq.txt"
    pfam_db = "data/pfamseq.db"
    pfam_name = "pfamseq"

    print("Importing Pfam sequences with full-text search.")
    import_pfam(pfam_file, pfam_db, pfam_name)

    print("Done!")


main()
