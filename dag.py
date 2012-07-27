#!/usr/bin/env python

import sys
import os

import config


def after_colon(line):
    # macro for getting anything after the :
    return line.split(":", 1)[1].strip()

def read_until(handle, start):
    # read each line until it has a certain start, and then puts the start tag back
    while 1:
        pos = handle.tell()
        line = handle.readline()
        if not line:
            break
        if line.startswith(start):
            handle.seek(pos)
            return
    raise EOFError, "%s tag cannot be found" % start


class OBOReader:
    """
    parse obo file, usually the most updated can be downloaded from
    http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo

    >>> reader = OBOReader()
    >>> for rec in reader:
            print rec
    """

    def __init__(self, obo_file="../../GO_data/gene_ontology_ext.obo"):
        try:
            self._handle = file(obo_file)
        except:
            print >> sys.stderr
            sys.exit(1)

    def __iter__(self):
        line = self._handle.readline()
        if not line.startswith(config.term_tag):
            read_until(self._handle, config.term_tag)
        while 1:
            yield self.next()

    def next(self):
        lines = []
        line = self._handle.readline()
        if not line or line.startswith(config.typedef_tag):
            raise StopIteration
        
        # read until the next tag and save everything in between
        while 1:
            pos = self._handle.tell() # save current postion for roll-back
            line = self._handle.readline()
            if line.startswith(config.typedef_tag) or line.startswith(config.term_tag):
                self._handle.seek(pos) # roll-back
                break
            lines.append(line)

        rec = TermNode()
        for line in lines:
            if line.startswith("id:"):
                rec.id = after_colon(line)
            if line.startswith("alt_id:"):
                rec.alt_ids.append(after_colon(line))
            elif line.startswith("name:"):
                rec.name = after_colon(line)
            elif line.startswith("namespace:"):
                rec.namespace = after_colon(line)
            elif line.startswith("is_a:"):
                #rec._parents.append(after_colon(line).split()[0])
                rec._parents[after_colon(line).split()[0]] = 0.8
            elif line.startswith("relationship: part_of"):
                #rec._parents.append(after_colon(line).split()[1])
                rec._parents[after_colon(line).split()[1]] = 0.6
            elif line.startswith("is_obsolete:") and after_colon(line)=="true":
                rec.is_obsolete = True
        return rec


class TermNode:
    def __init__(self):
        self.id = ""             # GO:xxxxxx
        self.name = ""           # description
        self.namespace = ""      # BP, CC, MF
        self._parents = {}       # parents term dict: key is parent term id, value is weight, if is_a then weight=0.8, else if part_of, weight=0.6 
        self.parents  = []       # parent records
        self.children = []       # children records
        self.level = -1          # distance from root node
        self.is_obsolete = False # is_obsolete
        self.alt_ids = []        # alternative identifiers

    def __str__(self):
        obsolete = "obsolete" if self.is_obsolete else ""
        return "%s\tlevel-%02d\t%s [%s] %s" % (self.id, self.level, self.name, self.namespace, obsolete)

    def __repr__(self):
        return self.id

    def has_parent(self, term):
        """ 
        Decides if term is ancester of this node
        """
        for p in self.parents:
            if p.id==term or p.has_parent(term):
                return True
        return False

    def has_child(self, term):
        """ 
        Decides if term is descendant of this node
        """
        for p in self.children:
            if p.id==term or p.has_child(term):
                return True
        return False

    def get_ancestors(self):
        """ 
        Get all ancestors of the node
        """
        all_parents = set()
        for p in self.parents:
            all_parents.add(p.id)
            all_parents |= p.get_ancestors()
        return all_parents

    def get_descendants(self):
        """
        Get all descendants of the node
        """
        all_children = set()
        for p in self.children:
            all_children.add(p.id)
            all_children |= p.get_descendants()
        return all_children

    def get_all_parent_edges(self):
        all_parent_edges = set()
        for p in self.parents:
            all_parent_edges.add((self.id, p.id))
            all_parent_edges |= p.get_all_parent_edges()
        return all_parent_edges

    def get_all_child_edges(self):
        all_child_edges = set()
        for p in self.children:
            all_child_edges.add((p.id, self.id))
            all_child_edges |= p.get_all_child_edges()
        return all_child_edges

    def is_root(self):
        return self.level==0


class DAG:
    def __init__(self, fpath):
        self.nodes = []
        self.id_node = {}
        self.load_file(fpath)
        self.offspring = {}
        self.populate_offspring()

    def load_file(self, fpath):
        print >>sys.stderr, "load obo file %s" % fpath
        obo_reader = OBOReader(fpath)
        for rec in obo_reader:
            if (rec.namespace==config.NAMESPACE) and (not rec.is_obsolete):
                self.nodes.append(rec)
                self.id_node[rec.id] = rec
        self.populate_terms()
        print >>sys.stderr, "%d terms in NAMESPACE %s" % (len(self.id_node), config.NAMESPACE)

    def populate_terms(self):
        def depth(rec):
            if rec.level < 0:
                if not rec.parents:
                    rec.level = 0
                else:
                    rec.level = min(depth(rec) for rec in rec.parents) + 1
            return rec.level

        for id in self.id_node.keys():
            rec = self.id_node[id]
            rec.parents = [self.id_node[p] for p in rec._parents]

        # populate children and levels
        for id in self.id_node:
            rec = self.id_node[id]
            for p in rec.parents:
                p.children.append(rec)
            if rec.level < 0:
                depth(rec)

    # Generate a dict of offsprings including the node itself for each node
    def populate_offspring(self):
        fpath = config.offspring_fpath
        if os.path.exists(fpath):
            ifile = open(fpath, "r")
            for line in ifile:
                line = line.strip()
                arr = line.split("\t")
                term = arr[0]
                self.offspring[term] = {}
                arr = arr[1].split(",")
                for t in arr:
                    self.offspring[term][t] = 1
            ifile.close()
        else:
            ofile = open(fpath, "w")
            for id in self.id_node:
                self.offspring[id] = {}
                rec = self.id_node[id]
                children = list(rec.get_descendants())
                children.append(id)
                for c in children:
                    self.offspring[id][c] = 1
                ofile.write(id + "\t" + ",".join(list(children)) + "\n")
            ofile.close()

    # draw AMIGO style network, lineage containing one query record
    def draw_lineage(self, recs, nodecolor="mediumseagreen", edgecolor="lightslateblue", dpi=96, verbose=False, lineage_img="GO_lineage.png"):
        pass

    # find the Lowest Common Ancestor of node1 and node2
    def lca(self, term1, term2):
        arr = []
        if term1==term2:
            arr.append(self.id_node[term1])
            return arr

        node1 = self.id_node[term1]
        node2 = self.id_node[term2]
        for node in self.nodes:
            child_dict = self.offspring[node.id]
            if (term1 in child_dict) and (term2 in child_dict):
                arr.append(node.id)
        return arr

    def get_root(self):
        for node in self.nodes:
            if node.is_root():
                return node


def main():
    dag = DAG(config.go_fpath)

if __name__ == "__main__":
    main()
    
