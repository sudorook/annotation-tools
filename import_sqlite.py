#! /usr/bin/env python3
""" Import text files into databases. """

import os
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

def import_pfam(file_in, db_out, table_name, fts=False):
    """ Import the Pfam protein sequences. """
    fields = ",".join(["[pfamseq_acc]", "[pfamseq_id]", "[seq_version]",
                       "[crc64]", "[md5]", "[description]", "[evidence]",
                       "[length]", "[species]", "[taxonomy]", "[is_fragment]",
                       "[sequence]", "[updated]", "[created]", "[ncbi_taxid]",
                       "[auto_architecture]", "[treefam_acc]", "[swissprot]"])

    conn = sqlite3.connect(db_out)
    cur = conn.cursor()

    # initialize the database
    cmd = ""
    if fts:
        cmd += "PRAGMA synchronous = OFF;\n"
        cmd += "PRAGMA journal_mode = MEMORY;\n"
        cmd += "BEGIN TRANSACTION;\n"
        cmd += "DROP TABLE IF EXISTS `%s`;\n" % (table_name)
        cmd += "CREATE VIRTUAL TABLE `%s` using fts4(\n" % (table_name)
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
        cmd += ",  --  PRIMARY KEY (`pfamseq_acc`)\n"
        cmd += ",  --  UNIQUE (`pfamseq_acc`)\n"
        cmd += ",  --  KEY `ncbi_taxid` (`ncbi_taxid`)\n"
        cmd += ",  --  KEY `crc64` (`crc64`)\n"
        cmd += ",  --  KEY `pfamseq_id` (`pfamseq_id`)\n"
        cmd += ",  --  KEY `pfamseq_architecture_idx` (`auto_architecture`)\n"
        cmd += ",  --  KEY `pfamseq_acc_version` (`pfamseq_acc`,`seq_version`)\n"
        cmd += ",  --  KEY `md5` (`md5`)\n"
        cmd += ",  --  KEY `pfamseq_tax_idx` (`taxonomy`)\n"
        cmd += ",  --  KEY `fk_pfamseq_evidence1_idx` (`evidence`)\n"
        cmd += ",  --  KEY `fragment_idx` (`is_fragment`)\n"
        cmd += ",  --  KEY `swissprot_idx` (`swissprot`)\n"
        cmd += ",  --  CONSTRAINT `FK_pfamseq_1` FOREIGN KEY (`ncbi_taxid`) REFERENCES `ncbi_taxonomy` (`ncbi_taxid`) ON DELETE CASCADE ON UPDATE NO ACTION\n"
        cmd += ",  --  CONSTRAINT `fk_pfamseq_evidence1` FOREIGN KEY (`evidence`) REFERENCES `evidence` (`evidence`) ON DELETE NO ACTION ON UPDATE NO ACTION\n"
        cmd += ");\n"
        cmd += "END TRANSACTION;"
    else:
        cmd += "PRAGMA synchronous = OFF;\n"
        cmd += "PRAGMA journal_mode = MEMORY;\n"
        cmd += "BEGIN TRANSACTION;\n"
        cmd += "DROP TABLE IF EXISTS `%s`;\n" % (table_name)
        cmd += "CREATE TABLE `%s` (\n" % (table_name)
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
        cmd += "CREATE INDEX \"idx_pfamseq_ncbi_taxid\" ON \"%s\" (`ncbi_taxid`);\n"% (table_name)
        cmd += "CREATE INDEX \"idx_pfamseq_crc64\" ON \"%s\" (`crc64`);\n"% (table_name)
        cmd += "CREATE INDEX \"idx_pfamseq_pfamseq_id\" ON \"%s\" (`pfamseq_id`);\n"% (table_name)
        cmd += "CREATE INDEX \"idx_pfamseq_pfamseq_architecture_idx\" ON \"%s\" (`auto_architecture`);\n"% (table_name)
        cmd += "CREATE INDEX \"idx_pfamseq_pfamseq_acc_version\" ON \"%s\" (`pfamseq_acc`,`seq_version`);\n"% (table_name)
        cmd += "CREATE INDEX \"idx_pfamseq_md5\" ON \"%s\" (`md5`);\n"% (table_name)
        cmd += "CREATE INDEX \"idx_pfamseq_pfamseq_tax_idx\" ON \"%s\" (`taxonomy`);\n"% (table_name)
        cmd += "CREATE INDEX \"idx_pfamseq_fk_pfamseq_evidence1_idx\" ON \"%s\" (`evidence`);\n"% (table_name)
        cmd += "CREATE INDEX \"idx_pfamseq_fragment_idx\" ON \"%s\" (`is_fragment`);\n"% (table_name)
        cmd += "CREATE INDEX \"idx_pfamseq_swissprot_idx\" ON \"%s\" (`swissprot`);\n"% (table_name)
        cmd += "END TRANSACTION;"
    print(cmd)
    cur.executescript(cmd)

    with open(file_in, "r") as handle:
        for line in handle.readlines():
            row = tuple([field.strip() for field in line.split("\t")])
            print(row)
            cmd = "INSERT INTO %s(" % (table_name) + fields + ")\n"
            cmd += "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);"
            cur.execute(cmd, row)
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
    nr_prot_fts_db = "data/nr_fts4.db"
    nr_prot_fts_name = "nr_fts"

    print("Importing non-redundant (nr) protein database.")
    import_nr(nr_prot_file, nr_prot_db, nr_prot_name, fts=False)

    print("Importing non-redundant (nr) database with full-text search.")
    import_nr(nr_prot_file, nr_prot_fts_db, nr_prot_fts_name, fts=True)

    # Pfam
    pfam_file = "raw_data/pfamseq.txt"
    pfam_db = "data/pfamseq.db"
    pfam_name = "pfamseq"
    pfam_fts_db = "data/pfamseq_fts4.db"
    pfam_fts_name = "pfamseq_fts"

    print("Importing Pfam sequences.")
    import_pfam(pfam_file, pfam_db, pfam_name, fts=False)

    print("Importing Pfam sequences with full-text search.")
    import_pfam(pfam_file, pfam_fts_db, pfam_fts_name, fts=True)

    print("Done!")


main()
