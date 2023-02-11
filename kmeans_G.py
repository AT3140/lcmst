#run src_Gh_kmeans.c => pauses => run ipynb => continue src
# python -u "d:\lcmst\kmeans_G.py"
from sklearn.cluster import KMeans
import numpy as np
import math
import matplotlib.pyplot as plt

CURR_GRAPH="curr_graph.txt"
CURR_INDS="curr_inds.txt"
KMeansClusterCenters="km_centres.txt"
KMeansLabels="k_labels.txt"
MAX=250

G=np.loadtxt(CURR_GRAPH,dtype='float')
I=np.loadtxt(CURR_INDS,dtype='float')

count_in=math.floor(MAX*0.1)
#kmeans=KMeans(n_clusters=count_in,n_init=MAX*10).fit(G)
kmeans=KMeans(n_clusters=count_in,init=I,n_init=1).fit(G)

np.savetxt(KMeansLabels,kmeans.labels_,fmt='%d')
np.savetxt(KMeansClusterCenters,kmeans.cluster_centers_,fmt='%.7f')