#!/usr/bin/env python

import re

filename1 = "Freq.txt"
filename2 = "mapChild.txt"

def get_terms(fpath):
    terms = {}
    ifile = open(fpath, "r")
    for line in ifile:
        line = line.strip()
        terms[line] = 1
    ifile.close()
    return terms

def read_mapChild(fpath):
    nodes = {}
    ifile = open(fpath, "r")
    for line in ifile:
        line = line.strip()
        arr = line.split("\t")
        node = arr[0].strip()
        nodes[node]=1
    ifile.close()
    return nodes

def read_geneassociation(fpath):
    terms = {}
    count = 0
    ifile = open(fpath, "r")
    for line in ifile:
        line = line.strip()
        m = re.search('\[(.+)\]', line)
        arr = m.group(1).split(",")
        for term in arr:
            count+=1
            terms[term.strip()] = 1
    ifile.close()
    print count
    return terms

if __name__ == "__main__":
    term_fpath = "terms"
    mapchild_fpath = "../mapChild.txt"
    geneasso_fpath = "../yeast_data/gene_annotation.txt"

    terms = get_terms(term_fpath)
    #nodes = read_mapChild(mapchild_fpath)
    gene_asso_terms = read_geneassociation(geneasso_fpath)

    print len(gene_asso_terms)

    for term in gene_asso_terms:
        if not term in terms:
            print term

