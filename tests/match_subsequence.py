#! /usr/bin/env python3
""" just testing """

import sqlite3


def truncate(sequence):
    """ Do nothing. Just a placeholder. """
    string = str(sequence)
    return string.split()[0]


def expand(sequence):
    """ Permute the string """
    string = str(sequence)
    string_length = len(string)
    return " ".join(
        [string[idx:string_length] for idx in range(0, string_length)]
    )


def connect_db(database):
    """ return db connection """
    conn = sqlite3.connect(database)
    conn.create_function("compress", 1, truncate)
    conn.create_function("uncompress", 1, expand)
    return conn


def disconnect_db(conn):
    """ disconnect from nr db """
    conn.close()


def search_db(conn, table, sequence):
    """ search nr db """
    cur = conn.cursor()
    cmd = "SELECT * FROM `%s` WHERE `sequence` = '%s';" % (table, sequence)
    cur.execute(cmd)
    res = cur.fetchone()
    if res:
        print("Match:\t" + res[0])
    else:
        print("Match:\tNONE")


def search_db_fts(conn, table, sequence):
    """ search nr db """
    cur = conn.cursor()
    cmd = "SELECT * FROM `%s` WHERE `sequence` MATCH '%s*';" % (
        table,
        sequence,
    )
    cur.execute(cmd)
    res = cur.fetchone()
    if res:
        print("Match:\t" + truncate(res[0]))
    else:
        print("Match:\tNONE")


def main():
    """ stuff """

    #
    # Search the nr database
    #
    database = "../data/nr.db"
    sequences = [
        "MSYLGQTDDRAIRQLACKCSVPREQTVHVAKGGPEDCRGAGGSAVGSSRMGPRWVLRMVTRHFILVQRSRSPVSLLPRRRFGRPTPPHPAPSPTLLASKQHIIISDRQKRPHRCALPKNFSPSKISSLRQAENPKMVP",
        "MSYLGQTDDRAIRQLACKCSVPREQTVHVAKGGPEDCRGAGGSAVGSSRMGPRWVLRMVTRHFILVQRSRSPVSLLPRRRFGRPTPPHPAPSPTLLASKQHIIISDRQKRPHRCALPKNFSPS",
        "CSVPREQTVHVAKGGPEDCRGAGGSAVGSSRMGPRWVLRMVTRHFILVQRSRSPVSLLPRRRFGRPTPPHPAPSPTLLASKQHIIISDRQKRPHRCALPKNFSPS",
        "CSVPREQTVHVAKGGPEDCRGAGRHFILVQRSRSPVSLLPRRRFGRPTPPHPAPSPTLLASKQHIIISDRQKRPHRCALPKNFSPS",
    ]

    conn = connect_db(database)

    table = "nr_fts"
    for sequence in sequences:
        print("Query:\t" + sequence)
        search_db_fts(conn, table, sequence)

    table = "nr"
    conn = connect_db(database)
    for sequence in sequences:
        print("Query:\t" + sequence)
        search_db(conn, table, sequence)
    print()
    disconnect_db(conn)

    #
    # Search the nr database
    #
    database = "../data/pfamseq.db"
    sequences = [
        "MLLGGGNALGQALIRLGAEEDIGFLAPRPPQDGWDAASLTQLLDDTRPDALINLAYYFDWFQAESVAEARLAAQERAVERLAELCQHHNIVLVQPSSYRVFDGSRATAYSEKDEPVPLGLRGQALWRIEQSVRATCPQHVMLRFGWLLDDSVDGTLGRFLARAELPEELLMADDRRGNPTPVDDAARVIISVLKQLDCAAPLWGTYHYAGHEATTPLALGQAILAEARNLHPLAIEAPTAQAHAARPDAAEEPQHAVLACKKILHTFGIKPRAWRAALPGLLDRFYRHG",
        "MLLGGGNALGQALIRLGAEEDIGFLAPRPPQDGWDAASLTQLLDDTRPDALINLAYYFDWFQAESVAEARLAAQERAVERLAELCQHHNIVLVQPSSYRVFDGSRATAYSEKDEPVPLGLRGQALWRIEQSVRATCPQHVMLRFGWLLDDSVDGTLGRFLARAELPEELLMADDRRGNPTPVDDAARVIISVLKQLDCAAPLWGTYHYAGHEATTPLALGQAILAEARNLHPLAIEAPTAQAHAARPDAAEEPQHAVLACKKIL",
        "EEDIGFLAPRPPQDGWDAASLTQLLDDTRPDALINLAYYFDWFQAESVAEARLAAQERAVERLAELCQHHNIVLVQPSSYRVFDGSRATAYSEKDEPVPLGLRGQALWRIEQSVRATCPQHVMLRFGWLLDDSVDGTLGRFLARAELPEELLMADDRRGNPTPVDDAARVIISVLKQLDCAAPLWGTYHYAGHEATTPLALGQAILAEARNLHPLAIEAPTAQAHAARPDAAEEPQHAVLACKKIL",
        "EEDIGFLAPRPPQDGINLAYYFDWFQAESVAEARLAAQERAVERLAELCQHHNIVLVQPSSYRVFDGSRATAYSEKDEPVPLGLRGQALWRIEQSVRATCPQHVMLRFGWLLDDSVDGTLGRFLARAELPEELLMADDRRGNPTPVDDAARVIISVLKQLDCAAPLWGTYHYAGHEATTPLALGQAILAEARNLHPLAIEAPTAQAHAARPDAAEEPQHAVLACKKIL",
    ]
    conn = connect_db(database)
    table = "pfamseq_fts"
    for sequence in sequences:
        print("Query:\t" + sequence)
        search_db_fts(conn, table, sequence)

    table = "pfamseq"
    conn = connect_db(database)
    for sequence in sequences:
        print("Query:\t" + sequence)
        search_db(conn, table, sequence)
    print()
    disconnect_db(conn)


main()
