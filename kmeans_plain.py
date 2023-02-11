# run this file => run src_kmeans_plain.c
from sklearn.cluster import KMeans
import numpy as np
import math
import matplotlib.pyplot as plt

CURR_GRAPH="./inst/lcmste250.15"
MAX=250
KMeansClusterCenters="km_centres_plain.txt"
KMeansLabels="k_labels_plain.txt"

G=np.loadtxt(CURR_GRAPH,dtype='float')

count_in=math.floor(MAX*0.1)
kmeans=KMeans(n_clusters=count_in,n_init=1).fit(G)
#kmeans=KMeans(n_clusters=count_in,init=I,n_init=1).fit(G)

# to be used as root[] in src.c
np.savetxt(KMeansLabels,kmeans.labels_,fmt='%d')
# to be used as ax[], ay[] in src.c
np.savetxt(KMeansClusterCenters,kmeans.cluster_centers_,fmt='%.7f')