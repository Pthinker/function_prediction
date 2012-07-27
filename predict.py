#!/usr/bin/env python

from operator import *
import heapq
import networkx as nx
import utils
from dag import *


def get_similarity(term1, term2):
    try:
        return sim_cache[term1][term2]
    except:
        sim = utils.get_db_sim(term1, term2, config.db_fpath)
        return sim

def evaluate(predicted_genes, network_annotated_gene, removed_genes, term_ic):
    recall_avg = 0.0;
    precision_avg = 0.0;

    recall_size = len(removed_genes)
    precision_size = len(removed_genes)

    for gene in removed_genes:
        orig_terms = network_annotated_gene[gene]
        predicted_terms = predicted_genes[gene]

        #print "Predicted(%s):%s" % (gene, predicted_terms)
        #print "Observed(%s):%s" % (gene, orig_terms)

        recall_base = 0.0
        recall_numerator = 0.0

        if len(predicted_genes[gene])>0:
            # Calculate recall
            for oterm in network_annotated_gene[gene]:
                recall_base += term_ic[oterm]
                max_sim = -1.0
                for pterm in predicted_genes[gene]:
                    sim = get_similarity(oterm, pterm)
                    if sim > max_sim:
                        max_sim = sim
                recall_numerator += max_sim
            recall_avg += recall_numerator/float(recall_base)

            # calculate precision
            precision_base = 0.0
            precision_numerator = 0.0
            for pterm in predicted_genes[gene]:
                precision_base += term_ic[pterm]
                max_sim = -1.0
                for oterm in network_annotated_gene[gene]:
                    sim = get_similarity(oterm, pterm)
                    if sim>max_sim:
                        max_sim = sim
                precision_numerator += max_sim
            precision_avg += precision_numerator/precision_base
        else:
            recall_size-=1
            precision_size-=1

    #recall_avg /= len(removed_genes)
    #precision_avg /= len(removed_genes)
    recall_avg /= recall_size
    precision_avg /= precision_size

    return recall_avg, precision_avg


def iterate_predict(network, annotated_genes, sim_cache, sim_terms, term_ic):
    predicted_genes = {}
    itr = 0
    last_sum = -1.0
    while itr<ITERATION:
        print "Iteration:%d" % itr
        total_sum = 0.0
        for gene in network.nodes():
            # Predict unannotated gene
            if not gene in annotated_genes:
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
                # Add similar terms to the candidate list
                cterms = term_dict.keys()
                for term in cterms:
                    for t in sim_terms[term]:
                        if not t in term_dict:
                            term_dict[t] = 1
                            candidate_terms.append(t)
                '''

                cterm_sim_sum = {}
                for cterm in candidate_terms:
                    sim_sum = 0.0
                    weight = 0
                    # For each neighbour of gene
                    for neighbour in network.neighbors(gene):
                        if neighbour in annotated_genes:
                            max_sim = -1.0
                            for nterm in annotated_genes[neighbour]:
                                if nterm == cterm:
                                    weight+=1
                                new_sim = get_similarity(cterm, nterm) + LAMDA * term_ic[cterm]
                                if new_sim>max_sim:
                                    max_sim = new_sim
                            sim_sum += max_sim
                        elif neighbour in predicted_genes:
                            max_sim = -1.0
                            for nterm in predicted_genes[neighbour]:
                                if nterm==cterm:
                                    weight+=1
                                new_sim = get_similarity(cterm, nterm) + LAMDA * term_ic[cterm]
                                if new_sim>max_sim:
                                    max_sim = new_sim
                            sim_sum += max_sim
                            
                    cterm_sim_sum[cterm] = sim_sum * weight

                #print heapq.nlargest(len(cterm_sim_sum), cterm_sim_sum.iteritems(), itemgetter(1))

                # Select top GONUMBER terms as predicted GO terms for the gene
                top_terms = heapq.nlargest(GONUMBER, cterm_sim_sum.iteritems(), itemgetter(1))
                if gene in predicted_genes:
                    del predicted_genes[gene]
                predicted_genes[gene] = []
                for rec in top_terms:
                    predicted_genes[gene].append(rec[0])
                    total_sum += rec[1]

        print "total sum:%f" % total_sum
        diff = int(total_sum) - int(last_sum)
        if diff==0:
            break
        else:
            last_sum = total_sum
        itr += 1

    return predicted_genes


def iterate_majority_voting(network, annotated_genes, sim_cache, sim_terms, term_ic):
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
                            if not term in neighbor_terms:
                                neighbor_terms[term] = 1
                            else:
                                neighbor_terms[term] += 1
                    elif neighbour in predicted_genes:
                        for term in predicted_genes[neighbour]:
                            if not term in neighbor_terms:
                                neighbor_terms[term] = 1
                            else:
                                neighbor_terms[term] += 1
                    
                # Select top GONUMBER terms as predicted GO terms for the gene
                top_terms = heapq.nlargest(GONUMBER, neighbor_terms.iteritems(), itemgetter(1))
                if gene in predicted_genes:
                    del predicted_genes[gene]
                predicted_genes[gene] = []
                for rec in top_terms:
                    predicted_genes[gene].append(rec[0])

                for rec in top_terms:
                    cterm = rec[0]
                    sim_sum = 0.0
                    for neighbour in network.neighbors(gene):
                        if neighbour in annotated_genes:
                            max_sim = -1.0
                            for nterm in annotated_genes[neighbour]:
                                new_sim = get_similarity(cterm, nterm) + LAMDA * term_ic[cterm]
                                if new_sim>max_sim:
                                    max_sim = new_sim
                            sim_sum += max_sim
                        elif neighbour in predicted_genes:
                            max_sim = -1.0
                            for nterm in predicted_genes[neighbour]:
                                new_sim = get_similarity(cterm, nterm) + LAMDA * term_ic[cterm]
                                if new_sim>max_sim:
                                    max_sim = new_sim
                            sim_sum += max_sim
                    total_sum += sim_sum

        print "total sum:%f" % total_sum
        diff = int(total_sum) - int(last_sum)
        if diff==0:
            break
        else:
            last_sum = total_sum
        itr += 1

    return predicted_genes


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


def cross_validation():
    group_genes = distribute_cross_validation(network_annotated_gene)

    #result_fpath = folder + "cv_%d_%s" % (GONUMBER, str(LAMDA))
    result_fpath = folder + "cv_%d_resnik" % (GONUMBER)
    ofile = open(result_fpath, "w")
    for i in range(CV):
        print "Cross validation %d......" % i
        annotated_gene_cv = {}
        removed_genes = group_genes[i]
        for gene in network_annotated_gene:
            if not gene in removed_genes:
                annotated_gene_cv[gene] = network_annotated_gene[gene]

        predicted_genes = iterate_predict(network, annotated_gene_cv, sim_cache, sim_terms, term_ic)
        #predicted_genes = iterate_majority_voting(network, annotated_gene_cv, sim_cache, sim_terms, term_ic)

        (recall, precision) = evaluate(predicted_genes, network_annotated_gene, removed_genes, term_ic)
        print "recall:%f, precision:%f" % (recall, precision)
        ofile.write(str(recall) + "," + str(precision) + "\n")
    ofile.close()


def predict_one(gene):
    annotated_gene = {}
    for g in network_annotated_gene:
        if g != gene:
            annotated_gene[g] = network_annotated_gene[g]
    predicted_genes = iterate_predict(network, annotated_gene, sim_cache, sim_terms, term_ic)
    return predicted_genes


if __name__ == "__main__":
    CV = 10
    ITERATION = 1
    LAMDA = 0.0;
    GONUMBER = 1
    result_fpath = config.folder + "cv_%d_%s" % (GONUMBER, str(LAMDA))

    dag = DAG(config.go_fpath)

    term_ic = utils.calculate_ic(gene_annotation, dag, config.ic_fpath)

    gene_annotation = get_annotation(config.annotation_fpath, config.filtered_annotation_fpath, dag.get_root())
    print "Number of annotated genes:%d" % len(gene_annotation)

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
    
    network_annotated_gene = {} # cross validation set of annotated genes
    for gene in gene_annotation:
        if gene in network.nodes():
            network_annotated_gene[gene] = gene_annotation[gene]
    print "Number of annotated genes in network:%d" % len(network_annotated_gene)

    sim_terms = utils.read_sim(config.sim_fpath)
    sim_cache = utils.generate_sim_cache(config.sim_cache_fpath, annotation_fpath, config.sim_fpath, config.db_fpath)

    for GONUMBER in range(1, 61):
        cross_validation()

    '''
    validate_gene = "YDL119C"
    predicted_genes = predict_one(validate_gene)
    (recall, precision) = evaluate(predicted_genes, network_annotated_gene, [validate_gene], term_ic)
    print "recall:%f, precision:%f" % (recall, precision)
    '''
