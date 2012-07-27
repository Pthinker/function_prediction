#!/usr/bin/env python

import re
import os
import math
import sqlite3
from operator import *
import heapq
from operator import itemgetter
import networkx as nx
import matplotlib.pyplot as plt

from dag import *
import config


# Create sqlite database connection
def create_connection(db_fpath):
    try:
        conn = sqlite3.connect(db_fpath)
    except:
        print "Cannot create database connection"
    return conn


# Return a dict with all the gene names in the network
def get_network_genes(fpath):
    network_genes = {}
    ifile = open(fpath, "r")
    for line in ifile:
        line = line.strip()
        arr = line.split("\t")
        network_genes[arr[0]] = 1
        network_genes[arr[1]] = 1
    ifile.close()
    return network_genes


# Create networkX network from network file
def create_network(fpath):
    print "Reading network from %s" % fpath
    graph = nx.Graph()
    ifile = open(fpath, "r")

    for line in ifile:
        line = line.strip()
        arr = line.split("\t")
        weight = float(arr[2])
        graph.add_edge(arr[0], arr[1], weight=weight)

    '''
    flag = False
    for line in ifile:
        if flag:
            line = line.strip()
            arr = line.split("\t")
            graph.add_edge(arr[2].strip(), arr[3].strip())
        else:
            if line.startswith("INTERACTOR_A"):
                flag = True
    '''
            
    ifile.close()

    return graph


# Visualize a network node and its connected nodes
def visualize_node(graph, node):
    if type(graph) is not nx.Graph:
        graph = create_network(graph)
    hub_ego = nx.ego_graph(graph, node)
    pos = nx.spring_layout(hub_ego)
    nx.draw(hub_ego, pos, node_color='b', node_size=900)
    nx.draw_networkx_nodes(hub_ego, pos, nodelist=[node], node_size=600, node_color='r')
    plt.show()


# Get term annotation for proteins from fpath
# Reture a dict with key as protein and value a list of terms
def get_annotation(fpath, filtered_fpath, root_term):
    annotation = {}
    if filtered_fpath and os.path.exists(filtered_fpath): # filterd_fpath exists, then read directly
        ifile = open(filtered_fpath, "r")
        for line in ifile:
            line = line.strip()
            arr = line.split("\t")
            gene = arr[0].strip()
            annotation[gene] = []
            anno_str = arr[1]
            arr = anno_str.split(",")
            for term in arr:
                annotation[gene].append(term.strip())
        ifile.close()
        return annotation

    ifile = open(fpath, "r")
    for line in ifile:
        if not line.startswith("!"):
            line = line.strip()
            arr = line.split("\t")
            if arr[8].strip()=='P':
                gene = arr[2].strip()
                term = arr[4].strip()
                if term == root_term: # Not consider root term
                    continue
                else:
                    if gene in annotation:
                        annotation[gene].add(term)
                    else:
                        annotation[gene] = set([term])
    ifile.close()

    if filtered_fpath:
        ofile = open(filtered_fpath, "w")
        for gene in annotation:
            ofile.write(gene + "\t" + ",".join(annotation[gene]) + "\n")
        ofile.close()

    return annotation


# Calculate GO term IC according to gene annotation
def calculate_ic(gene_annotation, dag, ic_fpath):
    if os.path.exists(ic_fpath):
        term_ic = {}
        ifile = open(ic_fpath, "r")
        for line in ifile:
            line = line.strip()
            arr = line.split(",")
            term_ic[arr[0]] = float(arr[1])
        ifile.close()
        return term_ic
    else:
        term_ic = {}
        term_count = {}
        for node in dag.nodes:
            term = node.id
            count = 0 
            for gene in gene_annotation:
                for t in gene_annotation[gene]:
                    if t in dag.offspring[term]:
                        count += 1
            term_count[term] = count
        root_count = term_count[dag.get_root().id]
        if config.filtered:
            filtered_term_count = {}
            for term in term_count:
                if term_count[term]>0:
                    filtered_term_count[term] = term_count[term]
            term_count = filtered_term_count

        for term in term_count:
            freq = term_count[term]/float(root_count)
            if freq<=0.0:
                ic = 0.001
            else:
                ic = -math.log(freq)
            term_ic[term] = ic

        ofile = open(ic_fpath, "w")
        for term in term_ic:
            ofile.write(term + "," + str(term_ic[term]) + "\n")
        ofile.close()
        return term_ic


# Get similarity of two terms from database
def get_db_sim(term1, term2, db_fpath):
    conn = create_connection(db_fpath)
    cur = conn.cursor()
    t = (term1, term2)
    cur.execute('select sim from similarity where term1=? and term2=?', t)
    rec = cur.fetchone()
    sim = rec[0]
    conn.close()
    return sim


# Get all similarities from a file
def read_sim(fpath):
    print "Reading from %s" % fpath
    sim_terms = {}
    ifile= open(fpath, "r")
    for line in ifile:
        line = line.strip()
        arr = line.split(",")
        if arr[0] in sim_terms:
            sim_terms[arr[0]][arr[1]] = float(arr[2])
        else:
            sim_terms[arr[0]] = {}
            sim_terms[arr[0]][arr[1]] = float(arr[2])
    ifile.close()
    return sim_terms


# Generate similarity cache
def generate_sim_cache(fpath, annotation_fpath, sim_fpath, db_fpath):
    # If cache exists, then just read
    if os.path.exists(fpath):
        sim_cache = read_sim(fpath)
        return sim_cache
    else: # if not exists, then create
        gene_annotation = get_annotation(config.annotation_fpath, config.filtered_annotation_fpath, dag.get_root().id)
        sim_terms = read_sim(sim_fpath)
        cache_terms = {}
        for gene in gene_annotation:
            for term in gene_annotation[gene]:
                cache_terms[term] = 1

        terms = cache_terms.keys()
        term_count = len(cache_terms)
        while 1:
            for term in terms:
                for t in sim_terms[term]:
                    cache_terms[t] = 1
            if len(cache_terms) == term_count:
                break
            else:
                term_count = len(cache_terms)
        print "Number of cached terms:%d" % len(cache_terms)
        ofile = open(fpath, "w")
        for term1 in cache_terms:
            for term2 in cache_terms:
                sim = get_db_sim(term1, term2, db_fpath)
                ofile.write(term1 + "," + term2 + "," + str(sim) + "\n")
        ofile.close()
        sim_cache = read_sim(fpath)
        return sim_cache


def filtered_most_sim():
    if os.path.exists(config.mostsim_fpath):
        return read_sim(config.mostsim_fpath)
    else:
        term_sim = read_sim(config.simcache_fpath)
        ofile = open(config.mostsim_fpath, "w")
        for term in term_sim:
            top_terms = heapq.nlargest(config.mostsim_num, term_sim[term].iteritems(), itemgetter(1))
            for rec in top_terms:
                ofile.write(term + "," + rec[0] + "," + str(rec[1]) + "\n")
        ofile.close()
        return read_sim(config.mostsim_fpath)

# Return a list of terms used in annotation
def get_annotated_terms(fpath):
    terms = set()
    fh = open(fpath, "r")
    for line in fh:
        line = line.strip()
        arr = line.split("\t")[1].split(",")
        for term in arr:
            terms.add(term)
    fh.close()
    return list(terms)

def main():
    dag = DAG(config.go_fpath)
    print "Root term: %s" % dag.get_root().id

    gene_annotation = get_annotation(config.annotation_fpath, config.filtered_annotation_fpath, dag.get_root().id)
    print "Annotated gene: %d" % len(gene_annotation)

    terms = get_annotated_terms(config.filtered_annotation_fpath)
    print len(terms)

    #term_ic = calculate_ic(gene_annotation, dag, config.ic_fpath)
    #print "Number of term ic: %d" % len(term_ic)

    # generate_sim_cache(config.sim_cache_fpath, config.annotation_fpath, config.sim_fpath, config.db_fpath)

if __name__ == "__main__":
    main()




