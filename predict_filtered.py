#!/usr/bin/env python

from operator import *
import heapq
import networkx as nx
import logging
from random import choice

import config
import utils
from dag import *


def evaluate(predicted_genes, removed_genes):
    #print "Removed %d genes" % len(removed_genes)
    #print "Predicted %d genes" % len(predicted_genes)

    recall_avg = 0.0;
    precision_avg = 0.0;

    recall_size = len(removed_genes)
    precision_size = len(removed_genes)

    for gene in predicted_genes:
        orig_terms = network_annotated_gene[gene]
        predicted_terms = predicted_genes[gene]

        #print "Predicted(%s):%s" % (gene, predicted_terms)
        #print "Observed(%s):%s" % (gene, orig_terms)

        recall_base = 0.0
        recall_numerator = 0.0

        # Calculate recall
        for oterm in network_annotated_gene[gene]:
            recall_base += wang_sim_cache[oterm][oterm]
        for oterm in network_annotated_gene[gene]:
            max_sim = -1.0
            for pterm in predicted_genes[gene]:
                sim = wang_sim_cache[oterm][pterm]
                if sim > max_sim:
                    max_sim = sim
            recall_numerator += max_sim
        #print "recall:%f" % (recall_numerator/recall_base)
        recall_avg += recall_numerator/float(recall_base)

        # calculate precision
        precision_base = 0.0
        precision_numerator = 0.0
        for pterm in predicted_genes[gene]:
            precision_base += wang_sim_cache[pterm][pterm]
        for pterm in predicted_genes[gene]:
            max_sim = -1.0
            for oterm in network_annotated_gene[gene]:
                sim = wang_sim_cache[pterm][oterm]
                if sim>max_sim:
                    max_sim = sim
            precision_numerator += max_sim
        #print "precision:%f" % (precision_numerator/precision_base)
        precision_avg += precision_numerator/precision_base

    recall_avg /= recall_size
    precision_avg /= precision_size
    return recall_avg, precision_avg


def iterate_weighted_predict(network, annotated_genes, sim_cache, sim_terms):
    predicted_genes = {}
    itr = 0
    last_sum = -1.0

    while itr<ITERATION:
        print "Iteration:%d" % itr
        total_sum = 0.0
        for gene in network.nodes():
            # Predict unannotated gene
            if not gene in annotated_genes:
                candidate_terms = get_candidate_terms(network, gene, annotated_genes, predicted_genes)
                
                cterm_sim_sum = {}
                for cterm in candidate_terms:
                    sim_sum = 0.0
                    # For each neighbour of gene
                    for neighbour in network.neighbors(gene):
                        weight = network.edge[gene][neighbour]['weight']
                        if neighbour in annotated_genes:
                            max_sim = -1.0
                            for nterm in annotated_genes[neighbour]:
                                new_sim = sim_cache[cterm][nterm]
                                if new_sim>max_sim:
                                    max_sim = new_sim
                            sim_sum += (max_sim*weight)
                        elif neighbour in predicted_genes:
                            max_sim = -1.0
                            for nterm in predicted_genes[neighbour]:
                                new_sim = sim_cache[cterm][nterm]
                                if new_sim>max_sim:
                                    max_sim = new_sim
                            sim_sum += (max_sim*weight)
                            
                    cterm_sim_sum[cterm] = sim_sum

                #print heapq.nlargest(len(cterm_sim_sum), cterm_sim_sum.iteritems(), itemgetter(1))

                # Select top GONUMBER terms as predicted GO terms for the gene
                top_terms = heapq.nlargest(GONUMBER, cterm_sim_sum.iteritems(), itemgetter(1))
                if gene in predicted_genes:
                    del predicted_genes[gene]
                predicted_genes[gene] = []
                for rec in top_terms:
                    predicted_genes[gene].append(rec[0])

        total_sum = compute_total_sim(network, annotated_genes, predicted_genes)
        print "total sum:%f" % total_sum
        diff = int(total_sum) - int(last_sum)
        if diff==0:
            break
        else:
            last_sum = total_sum
        itr += 1

    for gene in predicted_genes.keys():
        if len(predicted_genes[gene])==0:
            del predicted_genes[gene]

    return predicted_genes


def iterate_predict(network, annotated_genes, sim_cache, sim_terms, removed_genes):
    predicted_genes = {}
    itr = 0
    last_sum = -1.0
    
    result_list = list()

    while itr<ITERATION:
        print "Iteration:%d" % itr
        total_sum = 0.0
        for gene in network.nodes():
            # Predict unannotated gene
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
                                if new_sim>max_sim:
                                    max_sim = new_sim
                            if gene in predicted_genes:
                                weight = compute_gene_sim(predicted_genes[gene], annotated_genes[neighbour])
                            else:
                                weight = compute_gene_sim([cterm], annotated_genes[neighbour])
                            sim_sum += 1.0 * weight * max_sim
                        elif neighbour in predicted_genes:
                            max_sim = -1.0
                            for nterm in predicted_genes[neighbour]:
                                new_sim = sim_cache[cterm][nterm]
                                if new_sim>max_sim:
                                    max_sim = new_sim
                            if gene in predicted_genes:
                                weight = compute_gene_sim(predicted_genes[gene], predicted_genes[neighbour])
                            else:
                                weight = compute_gene_sim([cterm], predicted_genes[neighbour])
                            sim_sum += 1.0 * weight * max_sim
                    cterm_sim_sum[cterm] = sim_sum

                # Select top GONUMBER terms as predicted GO terms for the gene
                top_terms = heapq.nlargest(GONUMBER, cterm_sim_sum.iteritems(), itemgetter(1))
                if gene in predicted_genes:
                    del predicted_genes[gene]
                predicted_genes[gene] = []
                for rec in top_terms:
                    predicted_genes[gene].append(rec[0])
        '''
        total_sum = compute_total_sim(network, annotated_genes, predicted_genes)
        print "total sum:%f" % total_sum
        diff = int(total_sum) - int(last_sum)
        if diff==0:
            break
        else:
            last_sum = total_sum
        '''
        itr += 1

        for gene in predicted_genes.keys():
            if len(predicted_genes[gene])==0:
                del predicted_genes[gene]
        (recall, precision) = evaluate(predicted_genes, removed_genes)
        print "recall:%f, precision:%f" % (recall, precision)
        result_list.append(str(recall) + "," + str(precision))

    '''
    for gene in predicted_genes.keys():
        if len(predicted_genes[gene])==0:
            del predicted_genes[gene]
    return predicted_genes
    '''
    return result_list


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
    
    '''
    level2neighors = {}
    for neighbor in network.neighbors(gene):
        for n in network.neighbors(neighbor):
            if n != gene:
                if n in annotated_genes:
                    for term in annotated_genes[n]:
                        if not term in term_dict:
                            term_dict[term] = 1
                            candidate_terms.append(term)
                elif n in predicted_genes:
                    for term in predicted_genes[n]:
                        if not term in term_dict:
                            term_dict[term] = 1
                            candidate_terms.append(term)
    '''
    # Add similar terms to the candidate list
    cterms = term_dict.keys()
    for term in cterms:
        for t in sim_terms[term]:
            if not t in term_dict:
                term_dict[t] = 1
                candidate_terms.append(t)

    return candidate_terms 


def iterate_majority_voting(network, annotated_genes, sim_cache, sim_terms):
    predicted_genes = {}
    itr = 0
    last_sum = -1.0
    while itr<ITERATION:
        print "Iteration:%d" % itr
        total_sum = 0.0
        for gene in network.nodes():
            # Predict unannotated gene
            if not gene in annotated_genes:
                # Count neighbor terms
                neighbor_terms = {}
                for neighbour in network.neighbors(gene):
                    if neighbour in annotated_genes:
                        for term in annotated_genes[neighbour]:
                            neighbor_terms[term] = neighbor_terms.get(term, 0) + 1
                    elif neighbour in predicted_genes:
                        for term in predicted_genes[neighbour]:
                            neighbor_terms[term] = neighbor_terms.get(term, 0) + 1
                    
                # Select top GONUMBER terms as predicted GO terms for the gene
                top_terms = heapq.nlargest(GONUMBER, neighbor_terms.iteritems(), itemgetter(1))
                if gene in predicted_genes:
                    del predicted_genes[gene]
                predicted_genes[gene] = []
                for rec in top_terms:
                    predicted_genes[gene].append(rec[0])

        total_sum = compute_total_sim(network, annotated_genes, predicted_genes)
        diff = int(total_sum) - int(last_sum)
        if diff==0:
            break
        else:
            last_sum = total_sum
        itr += 1

    for gene in predicted_genes.keys():
        if len(predicted_genes[gene])==0:
            del predicted_genes[gene]

    return predicted_genes


# Compute the similarity between two genes. Each gene is represented using a
# list of terms
def compute_gene_sim(gene1_terms, gene2_terms):
    # Use max sim between two terms as the sim between two genes
    max_sim = -1.0
    for term1 in gene1_terms:
        for term2 in gene2_terms:
            sim = sim_cache[term1][term2]
            if sim > max_sim:
                max_sim = sim
    return max_sim
    '''

    # Use reciprocal max average
    row_avg = 0.0
    for term1 in gene1_terms:
        row_max = -1.0
        for term2 in gene2_terms:
            sim = sim_cache[term1][term2]
            if sim>row_max:
                row_max = sim
        row_avg += row_max

    col_avg = 0.0
    for term2 in gene2_terms:
        col_max = -1.0
        for term1 in gene1_terms:
            sim = sim_cache[term2][term1]
            if sim>col_max:
                col_max = sim
        col_avg += col_max

    return (row_avg+col_avg)/2
    '''

def compute_total_sim(network, annotated_genes, predicted_genes):
    sum = 0.0
    for gene in network.nodes():
        if gene in annotated_genes:
            gene_terms = annotated_genes[gene]
        elif gene in predicted_genes:
            gene_terms = predicted_genes[gene]

        for neighbor in network.neighbors(gene):
            if neighbor in annotated_genes:
                sim = compute_gene_sim(gene_terms, annotated_genes[neighbor])
            elif neighbor in predicted_genes:
                sim = compute_gene_sim(gene_terms, predicted_genes[neighbor])
            sum += sim
    return sum


# distrite all annotated gene to different cross validation group
def distribute_cross_validation(network_annotated_gene):
    group_genes = {}
    index = 0
    for gene in network_annotated_gene:
        group = index%CV
        if not group in group_genes:
            group_genes[group] = [gene]
        else:
            group_genes[group].append(gene)
        index += 1
    return group_genes


def remove_predict():
    for num in range(200, 3601, 200):
        removed_genes = []
        gene_list = network_annotated_gene.keys()
        for i in range(num):
            gene = choice(gene_list)
            removed_genes.append(gene)
            gene_list.remove(gene)

        annotated_gene_cv = {}
        for gene in gene_list:
            annotated_gene_cv[gene] = network_annotated_gene[gene]

        predicted_genes = iterate_predict(network, annotated_gene_cv, sim_cache, sim_terms)
        #predicted_genes = iterate_weighted_predict(network, annotated_gene_cv, sim_cache, sim_terms)
        #predicted_genes = iterate_majority_voting(network, annotated_gene_cv, sim_cache, sim_terms)
        
        '''
        (recall, precision) = evaluate(predicted_genes, removed_genes)
        print "removednum:%d, recall:%f, precision:%f" % (num,recall, precision)
        '''

def cross_validation():
    group_genes = distribute_cross_validation(network_annotated_gene)
    if ITERATION>1:
        iterate_str = "iterated_"
    else:
        iterate_str = ""

    result_folder = config.folder + "results/graphlet95/%s%s%s/" % (config.filtered, iterate_str, config.sim_metric)
    #result_folder = config.folder + "results/%s%smv_%s/" % (config.filtered, iterate_str, config.sim_metric)
    #result_folder = config.folder + "results/%s%sweighted_%s/" % (config.filtered, iterate_str, config.sim_metric)
    
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)

    result_fpath = result_folder + "%d.csv" % (GONUMBER)

    ofile = open(result_fpath, "w")
    for i in range(CV):
        print "Cross validation %d......" % i
        annotated_gene_cv = {}
        removed_genes = group_genes[i]
        for gene in network_annotated_gene:
            if not gene in removed_genes:
                annotated_gene_cv[gene] = network_annotated_gene[gene]
        print "%d genes annotated and %d genes to predict" % ( len(annotated_gene_cv), len(removed_genes) )

        result_list = iterate_predict(network, annotated_gene_cv, sim_cache, sim_terms, removed_genes)
        #predicted_genes = iterate_majority_voting(network, annotated_gene_cv, sim_cache, sim_terms)
        #predicted_genes = iterate_weighted_predict(network, annotated_gene_cv, sim_cache, sim_terms)
        
        #print_predictions(predicted_genes, removed_genes)

        #(recall, precision) = evaluate(predicted_genes, removed_genes)
        #print "recall:%f, precision:%f" % (recall, precision)

        #ofile.write(str(recall) + "," + str(precision) + "\n")
        for result in result_list:
            ofile.write(result + "\n")
        ofile.write("\n")
    ofile.close()


def predict_one(gene):
    annotated_gene = {}
    for g in network_annotated_gene:
        if g != gene:
            annotated_gene[g] = network_annotated_gene[g]
    predicted_genes = iterate_predict(network, annotated_gene, sim_cache, sim_terms, term_ic)
    return predicted_genes


def print_predictions(predicted_genes, removed_genes):
    for gene in predicted_genes:
        logging.info( "Predicted %s : %s" % (gene, str(predicted_genes[gene])) )
        logging.info( "Original %s : %s" % (gene, str(network_annotated_gene[gene])) )
        logging.info( "Neighbors:" )
        for neighbor in network.neighbors(gene):
            if neighbor in removed_genes:
                logging.info( "%s : []" % (neighbor) )
            else:
                logging.info( "%s : %s" % (neighbor, str(network_annotated_gene[neighbor])) )
        logging.info("----------------------------------------------------")


if __name__ == "__main__":
    CV = 10
    ITERATION = 10
    GONUMBER = 1

    logging.basicConfig(filename='predict.log', level=logging.DEBUG)

    dag = DAG(config.go_fpath)

    gene_annotation = utils.get_annotation(config.annotation_fpath, config.filtered_annotation_fpath, dag.get_root().id)
    print "Number of annotated genes:%d" % len(gene_annotation)

    term_ic = utils.calculate_ic(gene_annotation, dag, config.ic_fpath)

    network = utils.create_network(config.network_fpath)
    print "Number of nodes in network:%d" % network.number_of_nodes()

    # Remove unannotated gene from network
    for node in network.nodes():
        if not node in gene_annotation:
            network.remove_node(node)
    print "Number of nodes in network after removing unannotated genes:%d" % network.number_of_nodes()

    # Remove individual nodes by get the largest indepedent connected component
    network = nx.connected_component_subgraphs(network)[0]
    print "Number of nodes in network after removing individual genes:%d" % network.number_of_nodes()
    print "Number of edges in network after removing individual genes:%d" % network.number_of_edges()

    network_annotated_gene = {} # cross validation set of annotated genes
    for gene in gene_annotation:
        if gene in network.nodes():
            network_annotated_gene[gene] = gene_annotation[gene]
    print "Number of annotated genes in network:%d" % len(network_annotated_gene)

    sim_terms = utils.filtered_most_sim()
    sim_cache = utils.read_sim(config.simcache_fpath)
    wang_sim_cache = utils.read_sim(config.folder + "filtered_wang_sim.csv")

    #remove_predict()

    for GONUMBER in range(1, 11):
        cross_validation()
