# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 11:08 2020
@author: Corey Williams
"""
## Import packages
import pandas as pd
import leidenalg
import sklearn.neighbors
import igraph as ig
import csv
import os

## Input variables
INPUT_FOLDER = os.path.dirname(os.path.realpath(__file__))
EXPRESSION_FILENAME = "/Concat_Transformed.csv"
PANEL_FILENAME = "/panel.csv"
CLUSTERS_FILENAME = "/clusters.csv"
KNN_K = 100
N_ITERATIONS = -1
SEED = None

## Get knn graph
#load exprs file
cell_exprs_in = pd.read_csv(INPUT_FOLDER+EXPRESSION_FILENAME)
#take only clustering vars
panel = pd.read_csv(INPUT_FOLDER+PANEL_FILENAME)
clustering_antigens = panel['Antigen'][panel['Clustering']==1]
clustering_metals = panel['Metal'][panel['Clustering']==1]
clustering_vars = clustering_antigens+"_"+clustering_metals
#remove problematic characters
clustering_vars = clustering_vars.str.replace("-","_")
clustering_vars = clustering_vars.str.replace(".","_")
#subset for clustering variables
cell_exprs = cell_exprs_in[clustering_vars.values]
#get knn graph as sparse matrix
knn_out = sklearn.neighbors.kneighbors_graph(cell_exprs, KNN_K,
                                             mode = "distance")

# ## Do Leiden clustering
#convert graph matrix to igraph object
sources, targets = knn_out.nonzero()
edgelist = list(zip(sources.tolist(), targets.tolist()))
weights = 1/knn_out[sources,targets]
leiden_graph = ig.Graph(edgelist,edge_attrs={'weight': weights})
#Do Leiden clustering
leiden_out = leidenalg.find_partition(leiden_graph,
                                      leidenalg.ModularityVertexPartition,
                                      n_iterations=N_ITERATIONS,seed=SEED)
    
## Output csv file
with open((INPUT_FOLDER + CLUSTERS_FILENAME),'w',newline='') as csvfile:
    cluster_writer = csv.writer(csvfile, delimiter = ',')
    cluster_writer.writerow(leiden_out.membership)