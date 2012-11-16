#!/usr/bin/env python


DEBUG = True

CV = 10

# Gene Ontology related
#go_fpath = "../GO_data/gene_ontology_ext.obo"
go_fpath = '../GO_data/goslim_generic.obo'
NAMESPACE = "biological_process"
typedef_tag = "[Typedef]"
term_tag = "[Term]"

offspring_fpath = "data/bp_offspring.txt"


#METHOD = "mv"
METHOD = "weighted_mv"


# Intermedia data 
folder = "yeast_data/"
#folder = "pfalciparum_data/"

# Network
if folder.startswith("yeast"):
    network_fpath = "../Yeast_data/yeastnet2.gene.txt"
elif folder.startswith("pfalciparum"):
    network_fpath = "../P.Falciparum_data/above-3-links.3D7.filtered"

# Gene Annotation
if folder.startswith("yeast"):
    #annotation_fpath = "../Yeast_data/gene_association.sgd"
    annotation_fpath = "../Yeast_data/go_slim_mapping.tab"
elif folder.startswith("pfalciparum"):
    annotation_fpath = "../P.Falciparum_data/gene_association.GeneDB_Pfalciparum"


# Whether to remove terms that do not exist in annotation
filtered = "filtered_"

# which similarity metric to use
sim_metric = "wang"

filtered_annotation_fpath = folder + "gene_annotation.txt"

ic_fpath = folder + filtered + "ic.csv"

mostsim_num = 10
mostsim_fpath = folder + filtered + sim_metric + "_" + "%dmostsim.csv" % mostsim_num

if filtered:
    simcache_fpath = folder + filtered + sim_metric + "_" + "sim.csv"
else:
    simcache_fpath = folder + sim_metric + "_" + "simcache.csv"

db_fpath = folder + "data.db"

mfgo_adj_fpath = folder + "mfgo_adj.csv"

