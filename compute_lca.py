#!/usr/bin/env python

import threading
from dag import *
from multiprocessing import Process

class LCAThread(threading.Thread):
    def __init__(self, i, sublist):
        self.id = i
        self.sublist = sublist
        threading.Thread.__init__(self)

    def run(self):
        print "Thread %d starts running, handle %d terms" % (self.id, len(self.sublist))
        ofpath = folder + "lca%d.tsv" % self.id
        ofile = open(ofpath, "w")
        for term1 in self.sublist:
            for term2 in terms:
                arr = dag.lca(term1, term2)
                ofile.write(term1 + "\t" + term2 + "\t" + ",".join(arr) + "\n")
        ofile.close()
        print "Thread %d finished" % self.id

def get_lca(i, sublist):
    print "Thread %d starts running, handle %d terms" % (i, len(sublist))
    ofpath = folder + "lca%d.tsv" % i
    ofile = open(ofpath, "w")
    for term1 in sublist:
        for term2 in terms:
            if term1==term2:
                ofile.write(term1 + "\t" + term2 + "\t" + term1 + "\n")
            else:
                arr = dag.lca(term1, term2)
                ofile.write(term1 + "\t" + term2 + "\t" + ",".join(arr) + "\n")
    ofile.close()
    print "Thread %d finished" % i

if __name__ == "__main__":
    THREAD_NUM = 12
    dag = DAG(config.go_fpath)

    terms = []
    ifile = open("terms", "r")
    for line in ifile:
        line = line.strip()
        terms.append(line)
    ifile.close()

    totalLength = len(terms)
    subLength = totalLength/THREAD_NUM+1;

    processes = []
    for i in range(1, THREAD_NUM+1):
        if i==THREAD_NUM:
            sublist = sub_terms[(i-1)*subLength:]
        else:
            sublist = sub_terms[(i-1)*subLength : i*subLength]
        p = Process(target=get_lca, args=(i, sublist))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()
