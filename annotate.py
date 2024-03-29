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

DELIMITER = "|"

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
        "-d",
        "--delimiter",
        dest="delimiter",
        default="|",
        help="Delimiter for fields in FASTA headers",
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
        cur.execute(
            "SELECT `accession.version` FROM nr_fts WHERE sequence MATCH ?",
            (sequence,),
        )
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


def get_accession_number(identifier):
    """ Parse the accession number from the header """
    global DELIMITER
    fields = identifier.split(DELIMITER)
    try:
        idx = fields.index("ref") + 1
        return fields[idx]
    except ValueError as error:
        print("Error finding accession: %s" % error)
        return False


def get_gi_number(identifier):
    """ Parse the gi number from the header """
    global DELIMITER
    fields = identifier.split(DELIMITER)
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
            "SELECT taxid FROM prot WHERE `accession.version` = ?",
            (accession_number,),
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
        print("Accession %s not found. Trying dead_prot." % accession_number)
        with sqlite3.connect(DEADPROT_DB) as conn:
            cur = conn.cursor()
            cur.execute(
                "SELECT taxid FROM dead_prot WHERE `accession.version` = ?",
                (accession_number,),
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
        cur.execute("SELECT taxid FROM prot WHERE gi = ?", (gi_number,))
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
                "SELECT taxid FROM dead_prot WHERE gi = ?", (gi_number,)
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
            "SELECT species,lineage FROM taxonomy WHERE taxid = ?", (taxid,)
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
    return alignment


def get_annotations_pfamid(pfam_id):
    """ Get annotations by using pfam ID as the query. """
    with sqlite3.connect(PFAMSEQ_DB) as conn:
        cur = conn.cursor()
        cur.execute(
            "SELECT pfamseq_id,description,species,taxonomy FROM pfamseq WHERE pfamseq_id = ?",
            (pfam_id,),
        )
        res = cur.fetchall()
        if res:
            if len(res) == 1:
                pfam_id = res[0][0]
                pfam_name = res[0][1]
                pfam_species = res[0][2]
                pfam_lineage = ",".join(
                    [field.strip() for field in res[0][3].split(";")]
                )
                header = "|".join(
                    [pfam_id, pfam_name, pfam_species, pfam_lineage]
                )
            else:
                print("Multiple records found.")
                header = ",".join(["|".join(record) for record in res])
        else:
            header = False
    return header


def get_annotations_no_pfamid(sequence):
    """ Get annotations by querying the sequence. """
    with sqlite3.connect(PFAMSEQ_FTS) as conn:
        cur = conn.cursor()
        cur.execute(
            "SELECT pfamseq_id,description,species,taxonomy FROM pfamseq_fts WHERE sequence MATCH ?",
            (sequence,),
        )
        res = cur.fetchall()
        if res:
            if len(res) == 1:
                pfam_id = res[0][0]
                pfam_name = res[0][1]
                pfam_species = res[0][2]
                pfam_lineage = ",".join(
                    [field.strip() for field in res[0][3].split(";")]
                )
                header = "|".join(
                    [pfam_id, pfam_name, pfam_species, pfam_lineage]
                )
                header = "|".join(res[0])
            else:
                print("Multiple records found.")
                header = ",".join(["|".join(record) for record in res])
        else:
            header = False
    return header


def annotate_pfam(alignment):
    """ Loop over sequences and fill in annotations using Pfam. """
    for record in alignment:
        if not record.id:
            header = get_annotations_no_pfamid(
                str(record.seq).replace("-", "")
            )
            record.id = header
        else:
            pfam_id = record.id.split("|")[0]
            header = get_annotations_pfamid(pfam_id)
            if not header:
                print("Pfam ID %s not found. Querying sequence." % pfam_id)
                header = get_annotations_no_pfamid(
                    str(record.seq).replace("-", "")
                )
            if header:
                record.id = header
                record.name = ""
                record.description = ""
    return alignment


def main():
    """ Main """
    options = parse_options()
    if not check_databases():
        sys.exit("Cannot find required database files. Exiting.")

    global DELIMITER
    DELIMITER = options.delimiter

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
