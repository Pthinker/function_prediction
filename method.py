#!/usr/bin/env python

import heapq
from operator import *

import utils
import config


def iterate_mv(network, annotated_genes, go_num):
    unannotated_genes = list()
    for gene in network.nodes():
        if not gene in annotated_genes:
            unannotated_genes.append(gene)

    predicted_genes = {}

    while len(unannotated_genes) > 0:
        for gene in network.nodes():
            if not gene in annotated_genes:
                neighbor_terms = {}
                for neighbour in network.neighbors(gene):
                    if neighbour in annotated_genes:
                        for term in annotated_genes[neighbour]:
                            neighbor_terms[term] = neighbor_terms.get(term, 0) + 1
                    elif neighbour in predicted_genes:
                        for term in predicted_genes[neighbour]:
                            neighbor_terms[term] = neighbor_terms.get(term, 0) + 1

                if len(neighbor_terms) > 0:
                    # Select top go_num terms as predicted GO terms for the gene
                    top_terms = heapq.nlargest(go_num, neighbor_terms.iteritems(), itemgetter(1))
                    if gene in predicted_genes:
                        del predicted_genes[gene]
                    predicted_genes[gene] = []
                    for rec in top_terms:
                        predicted_genes[gene].append(rec[0])
                    if gene in unannotated_genes:
                        unannotated_genes.remove(gene)
                else:
                    continue

    return predicted_genes


def iterate_weighted_mv(network, annotated_genes, go_num):
    predicted_genes = {}
    iter = 0
    last_sum = -1.0
    ITERATION = 20

    sim_cache = utils.read_sim("pfalciparum_data/modified_wang_sim.csv")

    while iter < ITERATION:
        total_sum = 0.0
        for gene in network.nodes():
            if not gene in annotated_genes:
                candidate_terms = get_candidate_terms(network, gene, annotated_genes, predicted_genes)
                
                cterm_sim_sum = {}
                for cterm in candidate_terms:
                    sim_sum = 0.0
                    # For each neighbour of gene
                    for neighbour in network.neighbors(gene):
                        if neighbour in annotated_genes:
                            max_sim = -1.0
                            for nterm in annotated_genes[neighbour]:
                                new_sim = sim_cache[cterm][nterm]
                                if new_sim > max_sim:
                                    max_sim = new_sim
                            if gene in predicted_genes:
                                weight = compute_gene_sim(predicted_genes[gene], annotated_genes[neighbour], sim_cache)
                            else:
                                weight = compute_gene_sim([cterm], annotated_genes[neighbour], sim_cache)
                            sim_sum += 1.0 * weight * max_sim
                        elif neighbour in predicted_genes:
                            max_sim = -1.0
                            for nterm in predicted_genes[neighbour]:
                                new_sim = sim_cache[cterm][nterm]
                                if new_sim > max_sim:
                                    max_sim = new_sim
                            if gene in predicted_genes:
                                weight = compute_gene_sim(predicted_genes[gene], predicted_genes[neighbour], sim_cache)
                            else:
                                weight = compute_gene_sim([cterm], predicted_genes[neighbour], sim_cache)
                            sim_sum += 1.0 * weight * max_sim
                        cterm_sim_sum[cterm] = sim_sum

                if len(candidate_terms) > 0:
                    # Select top go_num terms as predicted GO terms for the gene
                    top_terms = heapq.nlargest(go_num, cterm_sim_sum.iteritems(), itemgetter(1))
                    if gene in predicted_genes:
                        del predicted_genes[gene]
                    predicted_genes[gene] = []
                    for rec in top_terms:
                        predicted_genes[gene].append(rec[0])

        total_sum = compute_total_sim(network, annotated_genes, predicted_genes, sim_cache)
        diff = int(total_sum) - int(last_sum)
        if diff==0:
            break
        else:
            last_sum = total_sum
        iter += 1

    return predicted_genes


def get_candidate_terms(network, gene, annotated_genes, predicted_genes):
    # Add neighbor terms to candidate list
    candidate_terms = []
    term_dict = {}
    for neighbour in network.neighbors(gene):
        if neighbour in annotated_genes:
            for term in annotated_genes[neighbour]:
                if not term in term_dict:
                    term_dict[term] = 1
                    candidate_terms.append(term)
        elif neighbour in predicted_genes:
            for term in predicted_genes[neighbour]:
                if not term in term_dict:
                    term_dict[term] = 1
                    candidate_terms.append(term)

    return candidate_terms


def compute_gene_sim(gene1_terms, gene2_terms, sim_cache):
    """Compute the similarity between two genes. Each gene is represented using a list of terms
    """
    # Use max sim between two terms as the sim between two genes
    max_sim = -1.0
    for term1 in gene1_terms:
        for term2 in gene2_terms:
            sim = sim_cache[term1][term2]
            if sim > max_sim:
                max_sim = sim
    return max_sim


def compute_total_sim(network, annotated_genes, predicted_genes, sim_cache):
    sum = 0.0
    for gene in network.nodes():
        if gene in annotated_genes:
            gene_terms = annotated_genes[gene]
        elif gene in predicted_genes:
            gene_terms = predicted_genes[gene]
        
        for neighbor in network.neighbors(gene):
            if neighbor in annotated_genes:
                sim = compute_gene_sim(gene_terms, annotated_genes[neighbor], sim_cache)
            elif neighbor in predicted_genes:
                sim = compute_gene_sim(gene_terms, predicted_genes[neighbor], sim_cache)
            sum += sim
    return sum


def main():
    pass

if __name__ == "__main__":
    main()
