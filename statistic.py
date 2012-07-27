#!/usr/bin/env python

import networkx as nx

from dag import *
import config
import utils


# Compute the percentage of genes in the network that has at least one GO term
# in any one of its neighbour
def compute_term_in_neighbour_ratio():
    total_gene = len(network.nodes())
    count = 0
    for gene in network.nodes():
        terms = gene_annotation[gene]
        for neighbour in network.neighbors(gene):
            nterms = gene_annotation[neighbour]
            if len( set(terms).intersection(set(nterms)) )>0:
                count += 1
                break
    print "Percentage of genes in the network that has at least one GO term in any one of its neighbour: %f" % (float(count)/total_gene)


# Compute the average number of terms of genes in the network
def compute_avg_term_num():
    total_gene = len(network.nodes())
    num = 0
    for gene in network.nodes():
        num += len( gene_annotation[gene] )
    print "Average term num: %f" % (float(num)/total_gene)


# Compute the similarity of 2 gene based on their annotations using the maximum
# term sim
def compute_gene_sim_max(terms1, terms2):
    # Use max similarity between two terms as similarity of two genes
    max_sim = -1.0
    for t1 in terms1:
        for t2 in terms2:
            sim = sim_cache[t1][t2]
            if sim > max_sim:
                max_sim = sim
    return max_sim


# Compute the similarity of 2 gene based on their annotations using the total
# term similarities
def compute_gene_sim_total(terms1, terms2):
    total_sim = 0.0
    for t1 in terms1:
        for t2 in terms2:
            sim = sim_cache[t1][t2]
            total_sim += sim
    return total_sim


def compute_avg_sim():
    print "Gene, Neighbor sim avg, Non-neighbor sim avg"
    for gene in network.nodes():
        terms = gene_annotation[gene]
        neighbors = network.neighbors(gene)
        # Compute avg sim with its neighbors
        neighbor_num = len(neighbors)
        sim_avg = 0.0
        for neighbor in neighbors:
            nterms = gene_annotation[neighbor]
            sim = compute_gene_sim_total(nterms, terms)
            #sim = compute_gene_sim_max(nterms, terms)
            sim_avg += sim
        sim_avg /= neighbor_num

        # Compute the avg sim with non-neighbor gene
        count = 0
        non_sim_avg = 0.0
        for node in network.nodes():
            if not node in neighbors:
                count += 1
                nterms = gene_annotation[node]
                sim = compute_gene_sim_total(nterms, terms)
                #sim = compute_gene_sim_max(nterms, terms)
                non_sim_avg += sim
        non_sim_avg /= count
        print "%s, %f, %f" % (gene, sim_avg, non_sim_avg)


if __name__ == "__main__":
    dag = DAG(config.go_fpath)

    gene_annotation = utils.get_annotation(config.annotation_fpath, config.filtered_annotation_fpath, dag.get_root().id)

    term_ic = utils.calculate_ic(gene_annotation, dag, config.ic_fpath)

    network = utils.create_network(config.network_fpath)
    # Remove unannotated gene from network
    for node in network.nodes():
        if not node in gene_annotation:
            network.remove_node(node)
    # Remove individual nodes by get the largest indepedent connected component
    network = nx.connected_component_subgraphs(network)[0]

    sim_cache = utils.read_sim(config.simcache_fpath)

    #compute_term_in_neighbour_ratio()
    #compute_avg_term_num()

    compute_avg_sim()


