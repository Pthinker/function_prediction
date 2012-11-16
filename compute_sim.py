#!/usr/bin/env python

import sys
import os
import sqlite3
import thread
from multiprocessing import Process

import utils
from dag import *


def wang_sim(term1, term2):
    if term1==term2:
        return 1.0

    node1 = dag.id_node[term1]
    node2 = dag.id_node[term2]

    sv_dict1 = {term1:1}
    nlist = [node1]
    while len(nlist)>0:
        new_nlist = []
        for node in nlist:
            for p in node.parents:
                value = sv_dict1[node.id]*node._parents[p.id]
                if p.id in sv_dict1:
                    if value>sv_dict1[p.id]:
                        sv_dict1[p.id] = value
                else:
                    sv_dict1[p.id] = value
                new_nlist.append(p)
        nlist = new_nlist

    sv_dict2 = {term2:1}
    nlist = [node2]
    while len(nlist)>0:
        new_nlist = []
        for node in nlist:
            for p in node.parents:
                value = sv_dict2[node.id]*node._parents[p.id]
                if p.id in sv_dict2:
                    if value>sv_dict2[p.id]:
                        sv_dict2[p.id] = value
                else:
                    sv_dict2[p.id] = value
                new_nlist.append(p)
        nlist = new_nlist

    set1 = set(sv_dict1.keys())
    set2 = set(sv_dict2.keys())
    common_terms = set1.intersection(set2)
    ssum = 0.0
    for term in common_terms:
        ssum += sv_dict1[term] + sv_dict2[term]

    sv1 = 0.0
    for key in sv_dict1:
        sv1 += sv_dict1[key]
    sv2 = 0.0
    for key in sv_dict2:
        sv2 += sv_dict2[key]

    return ssum/(sv1+sv2)


def modified_wang_sim(term1, term2):
    #if term1 == term2:
    #    return 1.0

    node1 = dag.id_node[term1]
    node2 = dag.id_node[term2]

    sv_dict1 = {term1 : 1}
    nlist = [node1]
    while len(nlist) > 0:
        new_nlist = []
        for node in nlist:
            for p in node.parents:
                value = sv_dict1[node.id] * node._parents[p.id]
                if p.id in sv_dict1:
                    if value > sv_dict1[p.id]:
                        sv_dict1[p.id] = value
                else:
                    sv_dict1[p.id] = value
                new_nlist.append(p)
        nlist = new_nlist

    sv_dict2 = {term2 : 1}
    nlist = [node2]
    while len(nlist) > 0:
        new_nlist = []
        for node in nlist:
            for p in node.parents:
                value = sv_dict2[node.id] * node._parents[p.id]
                if p.id in sv_dict2:
                    if value > sv_dict2[p.id]:
                        sv_dict2[p.id] = value
                else:
                    sv_dict2[p.id] = value
                new_nlist.append(p)
        nlist = new_nlist

    set1 = set(sv_dict1.keys())
    set2 = set(sv_dict2.keys())
    common_terms = set1.intersection(set2)
    ssum = 0.0
    for term in common_terms:
        ssum += sv_dict1[term] + sv_dict2[term]

    sv1 = 0.0
    for key in sv_dict1:
        sv1 += sv_dict1[key]

    sv2 = 0.0
    for key in sv_dict2:
        sv2 += sv_dict2[key]

    return ssum/(sv1 + sv2 + 0.1)


def lin_sim(term1, term2, lca_list):
    if term1 == term2:
        return (1.0, term1)
    max_sim = -1.0
    max_lca = None
    for lca in lca_list:
        sim = (2 * term_ic[lca]) / (term_ic[term1] + term_ic[term2])
        if sim > max_sim:
            max_sim = sim
            max_lca = lca
    return (max_sim, max_lca)


def resnik_sim(term1, term2, lca_list):
    if term1 == term2:
        return (term_ic[term1], term1)
    max_sim = -1.0
    max_lca = None
    for lca in lca_list:
        sim = term_ic[lca]
        if sim > max_sim:
            max_sim = sim
            max_lca = lca
    return (max_sim, max_lca)


def get_term_sim(term_sublist, i):
    print "Starting thread %d\n" % i
    sim_fpath = config.folder + "resnik_sim%d.csv" % i
    ofile = open(sim_fpath, "w")
    conn = utils.create_connection(config.db_fpath)
    for term in term_sublist:
        cur = conn.cursor()
        t = (term,)
        for record in cur.execute('select node2,lca from lca where node1=?', t):
            term2 = record[0]
            lca_str = record[1]
            arr = lca_str.split(",")
            (max_sim, max_lca) = resnik_sim(term, term2, arr)
            ofile.write( term + "," + term2 + "," + str(max_sim) + "," + max_lca + "\n" )
    ofile.close()
    conn.close()
    print "get_term_sim %d finished!" % i


def get_filtered_term_sim(term_sublist, terms, i):
    print "Starting thread %d\n" % i
    sim_fpath = config.folder + "filtered_lin_sim%d.csv" % i
    ofile = open(sim_fpath, "w")
    conn = utils.create_connection(config.db_fpath)
    for term in term_sublist:
        for term2 in terms:
            cur = conn.cursor()
            t = (term, term2)
            for record in cur.execute('select lca from lca where node1=? and node2=?', t):
                lca_str = record[0]
                arr = lca_str.split(",")
                #(max_sim, max_lca) = resnik_sim(term, term2, arr)
                (max_sim, max_lca) = lin_sim(term, term2, arr)
                ofile.write( term + "," + term2 + "," + str(max_sim) + "," + max_lca + "\n" )
            '''
            sim = resnik_sim(term, term2)
            ofile.write( term + "," + term2 + "," + str(sim) + "\n" )
            '''
    ofile.close()
    conn.close()
    print "get_filtered_term_sim%d finished!" % i


def compute_wang_sim(subterms, terms, i):
    print "Starting process %d\n" % i
    sim_fpath = config.folder + "modified_wang_sim%d.csv" % i
    ofile = open(sim_fpath, "w")
    for term1 in subterms:
        for term2 in terms:
            sim = modified_wang_sim(term1, term2)
            ofile.write(term1 + "," + term2 + "," + str(sim) + "\n")
    ofile.close()
    print "process %d finished!" % i
    

# Compute all the pairs of terms similarity
def compute_sim():
    terms = term_ic.keys()

    total_length = len(terms)
    sub_length = total_length/(THREAD_NUM + 1)

    processes = []
    for i in range(1, THREAD_NUM+1):
        if i == THREAD_NUM:
            sublist = terms[(i-1)*sub_length : ]
        else:
            sublist = terms[(i-1)*sub_length : i*sub_length]

        p = Process(target=compute_wang_sim, args=(sublist, terms, i))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()


if __name__ == "__main__":
    THREAD_NUM = 12

    dag = DAG(config.go_fpath)
    gene_annotation = utils.get_annotation(config.annotation_fpath, config.filtered_annotation_fpath, dag.get_root())
    term_ic = utils.calculate_ic(gene_annotation, dag, config.ic_fpath)

    compute_sim()

