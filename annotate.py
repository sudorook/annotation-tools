#! /usr/bin/env python3
""" Annotate alignments and sequences using NCBI and Pfam databases. """

import sys
import os
import argparse
import time
import sqlite3
from Bio import SeqIO
from Bio import AlignIO


# Globals

PFAMSEQ_DB = "data/pfamseq.db"
PFAMSEQ_FTS = "data/pfamseq_fts4.db"
DEADPROT_DB = "data/dead_prot_accession2taxid.db"
PROT_DB = "data/prot_accession2taxid.db"
NR_DB = "data/nr.db"
NR_FTS = "data/nr_fts4.db"
TAXONOMY_DB = "data/taxonomy.db"


# Functions


def parse_options():
    """ Make the command line option parser. """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", dest="input", required=True, help="Input file"
    )
    parser.add_argument(
        "-a",
        "--alignment",
        dest="alignment",
        action="store_true",
        default=False,
        help="Flag for whether the input is a set of sequences or an alignment",
    )
    parser.add_argument(
        "-f",
        "--format",
        dest="format",
        default="fasta",
        help="File format of the input file",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help="Output file for annotated sequences",
    )
    parser.add_argument(
        "-t", "--type", dest="type", required=True, help="'pfam' or 'ncbi'"
    )
    return parser.parse_args()


def check_databases():
    """ Check that the required databases exist. """
    flag = True
    databases = [
        PFAMSEQ_DB,
        PFAMSEQ_FTS,
        DEADPROT_DB,
        PROT_DB,
        NR_DB,
        NR_FTS,
        TAXONOMY_DB,
    ]
    for database in databases:
        if not os.path.exists(database):
            print(database + ": not found")
            flag = False
            break
    return flag


def load_alignment(alignment_file, alignment_format="fasta"):
    """ Load the alignment and return it. """
    alignment = AlignIO.read(alignment_file, alignment_format)
    return alignment


def load_sequences(sequence_file, sequence_format="fasta"):
    """ Load the sequence and return it. """
    sequences = SeqIO.read(sequence_file, sequence_format)
    return sequences


def save_alignment(alignment, alignment_file, alignment_format):
    """ Save an alignment to disk. """
    with open(alignment_file, "w") as handle:
        AlignIO.write(alignment, handle, alignment_format)


def save_sequences(sequence, sequence_file, sequence_format):
    """ Save sequences to disk. """
    with open(sequence_file, "w") as handle:
        SeqIO.write(sequence, handle, sequence_format)


def get_sequence_ncbi_id(sequence):
    """ Get sequence accession number """
    with sqlite3.connect(NR_FTS) as conn:
        cur = conn.cursor()
        cur.execute("SELECT `accession.version` FROM nr_fts WHERE sequence MATCH '%s'" % sequence)
        res = cur.fetchall()
        if res:
            if len(res) == 1:
                accession = res[0][0]
            else:
                print("Multiple records found.")
                accession = ",".join([acc[0] for acc in res])
        else:
            accession = "unknown"
    return accession


def get_accession_number(identifier, delimiter="|"):
    """ Parse the accession number from the header """
    fields = identifier.split(delimiter)
    try:
        idx = fields.index("ref") + 1
        return fields[idx]
    except ValueError as error:
        print("Error finding accession: %s" % error)
        return False


def get_gi_number(identifier, delimiter="|"):
    """ Parse the gi number from the header """
    fields = identifier.split(delimiter)
    try:
        idx = fields.index("gi") + 1
        return fields[idx]
    except ValueError as error:
        print("Error finding gi: %s" % error)
        return False


def get_accession_taxid(accession_number):
    """ Get the TaxID using the accession number. """
    with sqlite3.connect(PROT_DB) as conn:
        cur = conn.cursor()
        cur.execute(
            "SELECT taxid FROM prot WHERE `accession.version` = '%s'"
            % accession_number
        )
        res = cur.fetchall()
        if res:
            if len(res) == 1:
                taxid = res[0][0]
            else:
                print("Multiple records found.")
                taxid = ",".join([acc[0] for acc in res])
        else:
            taxid = "unknown"
    if taxid == "unknown":
        print("accession %s not found. Trying dead_prot." % accession_number)
        with sqlite3.connect(DEADPROT_DB) as conn:
            cur = conn.cursor()
            cur.execute(
                "SELECT taxid FROM dead_prot WHERE `accession.version` = '%s'"
                % accession_number
            )
            res = cur.fetchall()
            if res:
                if len(res) == 1:
                    taxid = res[0][0]
                else:
                    print("Multiple records found.")
                    taxid = ",".join([acc[0] for acc in res])
            else:
                taxid = "unknown"
    return taxid


def get_gi_taxid(gi_number):
    """ Get the TaxID using the GI number. """
    with sqlite3.connect(PROT_DB) as conn:
        cur = conn.cursor()
        cur.execute("SELECT taxid FROM prot WHERE gi = '%s'" % gi_number)
        res = cur.fetchall()
        if res:
            if len(res) == 1:
                taxid = res[0][0]
            else:
                print("Multiple records found.")
                taxid = ",".join([i[0] for i in res])
        else:
            taxid = False
    if not taxid:
        print("GI %s not found. Trying dead_prot." % gi_number)
        with sqlite3.connect(DEADPROT_DB) as conn:
            cur = conn.cursor()
            cur.execute(
                "SELECT taxid FROM dead_prot WHERE gi = '%s'" % gi_number
            )
            res = cur.fetchall()
            if res:
                if len(res) == 1:
                    taxid = res[0][0]
                else:
                    print("Multiple records found.")
                    taxid = ",".join([i[0] for i in res])
            else:
                taxid = False
    return taxid


def get_species_lineage(taxid):
    """ Get the species and lineage fields from the taxonomy database. """
    with sqlite3.connect(TAXONOMY_DB) as conn:
        cur = conn.cursor()
        cur.execute(
            "SELECT species,lineage FROM taxonomy WHERE taxid = '%s'" % taxid
        )
        res = cur.fetchall()
        if res:
            if len(res) == 1:
                species = res[0][0]
                lineage = res[0][1].replace(";", ",").strip()
        else:
            species = "unknown"
            lineage = "unknown"
    return species, lineage


def annotate_ncbi(alignment):
    """ Loop over sequences and fill in annotations using nr. """
    for record in alignment:
        print(record.id)
        if not record.id:
            record.id = "ref|" + get_sequence_ncbi_id(
                str(record.seq).replace("-", "")
            )

        # Use an ID (accession or GI) to get the taxID
        seq_id = get_accession_number(record.id)
        if seq_id:
            taxid = get_accession_taxid(seq_id)
        else:
            seq_id = get_gi_number(record.id)
            if seq_id:
                taxid = get_gi_taxid(seq_id)
            else:
                record.id = "ref|" + get_sequence_ncbi_id(
                    str(record.seq).replace("-", "")
                )
                taxid = get_accession_taxid(seq_id)

        # Use the TaxID to get the species and lineage information
        if taxid:
            species, lineage = get_species_lineage(taxid)
        record.id = record.id + "|" + species + "|" + lineage
        record.description = ""
        print(record.id)

    return alignment


def annotate_pfam(alignment):
    """ Loop over sequences and fill in annotations using Pfam. """
    print(alignment)
    return alignment


def main():
    """ Main """
    options = parse_options()
    if not check_databases():
        sys.exit("Cannot find required database files. Exiting.")

    if options.alignment:
        input_data = load_alignment(options.input, options.format)
    else:
        input_data = load_sequences(options.input, options.format)

    start = time.time()
    if options.type == "pfam":
        output_data = annotate_pfam(input_data)
    elif options.type == "ncbi":
        output_data = annotate_ncbi(input_data)
    end = time.time()

    print(output_data)
    print("Annotation time: %f sec" % (end - start))

    if options.format == "fasta":
        output_format = "fasta-2line"
    else:
        output_format = options.format

    if options.alignment:
        save_alignment(output_data, options.output, output_format)
    else:
        save_sequences(output_data, options.output, output_format)


main()
