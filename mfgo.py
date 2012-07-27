#!/usr/bin/env python

import os
import math
import random
import logging
import networkx as nx

import utils
import config
import dag


def adjacency_matrix(network, fpath):
    print "Reading adjacency matrix..."
    matrix = dict()
    if os.path.exists(fpath):
        fh = open(fpath, "r")
        for line in fh:
            line = line.strip()
            arr = line.split("\t")
            if not arr[0] in matrix:
                matrix[arr[0]] = dict()
            try:
                matrix[arr[0]][arr[1]] = float(arr[2])
            except:
                print arr
        fh.close()
    else:
        for node1 in network.nodes():
            matrix[node1] = dict()
            for node2 in network.nodes():
                if node2 != node1:
                    node1_neighbours = set(network.neighbors(node1))
                    node2_neighbours = set(network.neighbors(node2))
                    if node2 in node1_neighbours:
                        matrix[node1][node2] = 1.0
                    else:
                        intersection_len = len(node1_neighbours & node2_neighbours)
                        union_len = len(node1_neighbours | node2_neighbours)
                        matrix[node1][node2] = float(intersection_len) / union_len
        fh = open(fpath, "w")
        for node1 in matrix:
            for node2 in matrix[node1]:
                fh.write("%s\t%s\t%f\n" % (node1, node2, matrix[node1][node2]))
        fh.close()
    print "Finished reading adjacency matrix."
    return matrix


# Get a list of terms used in annotation
def get_terms(annotated_genes):
    terms = set()
    for gene in annotated_genes:
        for term in annotated_genes[gene]:
            terms.add(term)
    return list(terms)
    

# Normalize an array to make sum(d^2)=1
def normalize(arr):
    sum = 0.0
    for i in arr:
        sum += i*i
    sum = math.sqrt(sum)
    for i in range(len(arr)):
        arr[i] /= sum
        

def distance(arr1, arr2):
    dis = 0.0
    for i in range(len(arr1)):
        dis += (arr1[i] - arr2[i]) ** 2
    return dis


def MFGO(network, annotated_genes, adjacency):
    terms = get_terms(annotated_genes)
    print "Number of terms in configuration: %d" % len(terms)

    term_index = dict()
    for index, term in enumerate(terms):
        term_index[term] = index

    # generate term configuration of each gene
    unannotated_genes = list()
    gene_configures = dict()
    for gene in network.nodes():
        if gene in annotated_genes:
            config = [0.0] * len(terms)
            for term in annotated_genes[gene]:
                config[term_index[term]] = 1.0
            gene_configures[gene] = config
        else:
            unannotated_genes.append(gene)
            config = [0.0] * len(terms)
            for i in range(len(terms)):
                config[i] = random.random() + 1.0
            normalize(config)
            gene_configures[gene] = config
    print "Generated term configurations of genes"
    
    threshold = 0.01
    d = 1
    new_config = [0.0] * len(terms)
    gene2update = None
    while d > threshold:
        for gene in unannotated_genes:
            for i in range(len(terms)):
                sum = 0.0
                for other_gene in adjacency[gene]:
                    sum += adjacency[gene][other_gene] * gene_configures[other_gene][i]
                new_config[i] = sum
            normalize(new_config)
            logging.info(new_config)
            d = distance(new_config, gene_configures[gene])
            if d > threshold:
                gene2update = gene
                break
        if d > threshold:
            for i in range(len(gene_configures[gene2update])):
                gene_configures[gene2update][i] = new_config[i]
            logging.info(gene2update)
            logging.info(d)

    logging.info(gene_configures)
    logging.info("--------------------------------------")


# Distribute annotated genes in the network to different cross validataion group
# Return a dict with key as group id and value as a list of genes in that group
def distribute_cross_validation(network_annotated_gene):
    group_genes = {}
    index = 0
    for gene in network_annotated_gene:
        group = index%config.CV
        if not group in group_genes:
            group_genes[group] = [gene]
        else:
            group_genes[group].append(gene)
        index += 1
    return group_genes


def cross_validation(network, network_annotated_genes, gene_annotation, adjacency):
    group_genes = distribute_cross_validation(network_annotated_genes)
    for i in range(config.CV):
        print "Cross validation %d......" % i
        removed_genes = group_genes[i] # Genes to be removed
        annotated_gene_cv = {} # Remaining annotated genes
        for gene in network_annotated_genes:
            if not gene in removed_genes:
                annotated_gene_cv[gene] = network_annotated_genes[gene]
        print "%d genes annotated and %d genes to predict as test" % ( len(annotated_gene_cv), len(removed_genes) )

        MFGO(network, annotated_gene_cv, adjacency)


def main():
    logging.basicConfig(filename='mfgo.log', level=logging.DEBUG)

    ontology = dag.DAG(config.go_fpath)

    gene_annotation = utils.get_annotation(config.annotation_fpath, config.filtered_annotation_fpath, ontology.get_root().id)

    network = utils.create_network(config.network_fpath)
    print "Number of nodes in network: %d" % network.number_of_nodes()
    # Remove individual nodes by getting the largest indepedent connected component
    network = nx.connected_component_subgraphs(network)[0]
    print "Number of nodes in network after removing individual genes: %d" % network.number_of_nodes()
    #print "Number of edges in network after removing individual genes:%d" % network.number_of_edges()

    # Annotated genes in the network
    network_annotated_genes = {}
    for node in network.nodes():
        if node in gene_annotation:
            network_annotated_genes[node] = gene_annotation[node]
    print "Number of annotated genes in network:%d" % len(network_annotated_genes)

    # adjacency matrix of any pair of nodes base on their neighbour sharing property
    adjacency = adjacency_matrix(network, config.mfgo_adj_fpath)

    cross_validation(network, network_annotated_genes, gene_annotation, adjacency)


if __name__ == "__main__":
    main()
