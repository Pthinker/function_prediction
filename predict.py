#!/usr/bin/env python

from operator import *
import networkx as nx
import heapq

import method
import utils
import config
import dag


def evaluate(predicted_genes, removed_genes, network_annotated_gene):
    recall_avg = 0.0
    precision_avg = 0.0

    recall_size = len(removed_genes)
    precision_size = len(removed_genes)

    for gene in removed_genes:
        # Calculate recall
        recall_base = 0.0
        recall_numerator = 0.0
        for oterm in network_annotated_gene[gene]:
            recall_base += wang_sim_cache[oterm][oterm]

            max_sim = -1.0
            for pterm in predicted_genes[gene]:
                sim = wang_sim_cache[oterm][pterm]
                if sim > max_sim:
                    max_sim = sim
            recall_numerator += max_sim
        recall_avg += recall_numerator / float(recall_base)

        # calculate precision
        precision_base = 0.0
        precision_numerator = 0.0
        for pterm in predicted_genes[gene]:
            precision_base += wang_sim_cache[pterm][pterm]

            max_sim = -1.0
            for oterm in network_annotated_gene[gene]:
                sim = wang_sim_cache[pterm][oterm]
                if sim > max_sim:
                    max_sim = sim
            precision_numerator += max_sim
        precision_avg += precision_numerator/precision_base

    recall_avg /= recall_size
    precision_avg /= precision_size

    return recall_avg, precision_avg


def distribute_cross_validation(network_annotated_gene):
    """distrite all annotated gene to different cross validation group
    """
    group_genes = {}
    index = 0
    for gene in network_annotated_gene:
        group = index % config.CV
        if not group in group_genes:
            group_genes[group] = [gene]
        else:
            group_genes[group].append(gene)
        index += 1
    return group_genes


def cross_validation(network, network_annotated_gene, gene_annotation, go_num): 
    # Put annotated gene into different cross validation groups
    group_genes = distribute_cross_validation(network_annotated_gene)

    recall_avg = 0.0
    precision_avg = 0.0
    for i in range(0, config.CV):
        print "Cross validation %d......" % i
        annotated_gene_cv = {}
        removed_genes = group_genes[i]
        for gene in network_annotated_gene:
            if not gene in removed_genes:
                annotated_gene_cv[gene] = network_annotated_gene[gene]
    
        if config.METHOD == "mv":
            predicted_genes = method.iterate_mv(network, annotated_gene_cv, go_num)
        elif config.METHOD == "weighted_mv":
            predicted_genes = method.iterate_weighted_mv(network, annotated_gene_cv, go_num)

        (recall, precision) = evaluate(predicted_genes, removed_genes, network_annotated_gene)
        print "recall:%f, precision:%f" % (recall, precision)

        recall_avg += recall
        precision_avg += precision

    recall_avg /= config.CV
    precision_avg /= config.CV

    return (recall_avg, precision_avg)


def remove_unannotated_genes(network, gene_annotation):
    """Remove unannotated gene from network
    """
    for node in network.nodes():
        if not node in gene_annotation:
            network.remove_node(node)
    print "Number of nodes in network after removing unannotated genes:%d" % network.number_of_nodes()
    

def main():
    ontology = dag.DAG(config.go_fpath)

    #gene_annotation = utils.get_annotation(config.annotation_fpath, config.filtered_annotation_fpath, ontology.get_root())
    gene_annotation = utils.get_annotation(config.annotation_fpath, config.filtered_annotation_fpath, ontology.get_root())
    print "Number of annotated genes:%d" % len(gene_annotation)

    network = utils.create_network(config.network_fpath)
    print "Number of nodes in network:%d" % network.number_of_nodes()

    '''
    #remove_unannotated_genes(network, gene_annotation)

    # Remove individual nodes by get the largest indepedent connected component
    network = nx.connected_component_subgraphs(network)[0]
    print "Number of nodes in network after removing individual genes:%d" % network.number_of_nodes()
    print "Number of edges in network after removing individual genes:%d" % network.number_of_edges()

    network_annotated_gene = {} # cross validation set of annotated genes in the network
    for gene in gene_annotation:
        if gene in network.nodes():
            network_annotated_gene[gene] = gene_annotation[gene]
    print "Number of annotated genes in network:%d" % len(network_annotated_gene)

    fh = open(config.folder + "result/" + config.METHOD + ".csv", "w")
    for go_num in range(1, 11):
        (recall, precision) = cross_validation(network, network_annotated_gene, gene_annotation, go_num)
        fh.write("%f,%f\n" % (recall, precision))
    fh.close()
    '''

if __name__ == "__main__":
    #wang_sim_cache = utils.read_sim(config.folder + "filtered_wang_sim.csv")

    main()

