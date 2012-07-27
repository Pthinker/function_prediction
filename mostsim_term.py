#!/usr/bin/env python

import sqlite3
import sys
import thread

def readIC(fpath):
    term_ic = {}
    ifile = open(fpath, "r")
    for line in ifile:
        line = line.strip()
        arr = line.split(",")
        term_ic[arr[0]] = float(arr[1])
    ifile.close()
    return term_ic

def create_conn():
    try:
        conn = sqlite3.connect(db_fpath)
    except:
        sys.exit(1)
    return conn

def get_most_sim_terms(terms, i):
    print "Starting thread %d" % i
    ofile = open(folder + "10simTerms%d.csv"%i, "w")
    conn = create_conn()
    cur = conn.cursor()
    for term in terms:
        t = (term,)
        cur.execute("select term2, sim from similarity where term1=? order by sim desc limit 10", t)
        for row in cur.fetchall():
            term2 = row[0]
            sim = row[1]
            ofile.write(term + "," + term2 + "," + str(sim) + "\n")
    conn.close()
    ofile.close()
    print "Thread %d finished!" % i

if __name__ == "__main__":
    folder = "yeast_data/"
    #folder = "../pfalciparum_data/"
    THREAD_NUM = 12
    ic_fpath = folder + "ic.csv"
    db_fpath = folder + "data.db"

    term_ic = readIC(ic_fpath)
    terms = term_ic.keys()

    totalLength = len(terms)
    subLength = totalLength/THREAD_NUM;
    for i in range(1, THREAD_NUM+1):
        if i==THREAD_NUM:
            sublist = terms[(i-1)*subLength:]
        else:
            sublist = terms[(i-1)*subLength : i*subLength]
        thread.start_new_thread( get_most_sim_terms, (sublist,i) )

    while 1:
        pass
