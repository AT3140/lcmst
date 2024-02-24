#run src_Gh_kmeans.c => pauses => run ipynb => continue src
# python -u "d:\lcmst\kmeans_G.py"
from sklearn.cluster import KMeans
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree

CURR_GRAPH="curr_graph.txt"
CURR_INDS="curr_inds.txt"
KMeansClusterCenters="km_centres.txt"
KMeansLabels="k_labels.txt"
MAX=250

G=np.loadtxt(CURR_GRAPH,dtype='float')
I=np.loadtxt(CURR_INDS,dtype='float')
#n_iter=math.floor(MAX)
count_in=math.floor(MAX*0.1)
n_iter=math.floor(MAX*6)
n_perturb=math.ceil(count_in/5)
# n_perturb=1
print(n_iter)

# def connect_G(G,I):    



#np.savetxt(KMeansLabels,kmeans.labels_,fmt='%d')
#np.savetxt(KMeansClusterCenters,kmeans.cluster_centers_,fmt='%.7f')
def kcc_to_inds(kcc,G):
    # Find nearest neighbors
    kdtree = cKDTree(G)
    _, nn_indices = kdtree.query(kcc)
    # Create internal nodes
    internal_nodes = set()
    for nn_index in nn_indices:
        internal_nodes.add(nn_index)
    return list(internal_nodes)

def buildTree(inds, G):
    num_G=G.shape[0]
    count_in = len(inds)

    # Create edges between internal nodes using Prim's algorithm
    # making changes here take care of indexing of nn_indices 0 based
    Icd = G[inds, :]
    visited = np.zeros(count_in, dtype=bool)
    distances = np.full(num_G, np.inf)
    parent = np.full(num_G, -1, dtype=int)
    
    start_index = 0
    distances[inds[start_index]] = 0
    
    for _ in range(count_in):
        # Find the unvisited node with the smallest distance
        unvisited_indices = np.where(~visited)[0]

        # TODO: obt distance from num_G array but return index enumerated from num_kcc array
        mi=-1
        min_= sys.maxsize
        pos=0
        for i in unvisited_indices :
            if(distances[inds[i]]<min_) :
                min_=distances[inds[i]]
                mi=pos
            pos=pos+1
        u = unvisited_indices[mi]
        visited[u] = True
        
        # Update distances and parent for neighbors of u
        kdtree2=cKDTree(Icd)
        neighbors = kdtree2.query_ball_point(Icd[u], 1)
        
        for v in neighbors:
            if not visited[v]:
                dist = np.linalg.norm(Icd[u] - Icd[v])
                if dist < distances[inds[v]]:
                    distances[inds[v]] = dist
                    #print("u:"+str(u)+" v:"+str(v))
                    parent[inds[v]] = inds[u]
    
    # Create edges between leaves and nearest internal nodes
    leaf_indices = np.where(~np.isin(np.arange(num_G), inds))[0]

    for leaf_index in leaf_indices:
        nn_index = inds[np.argmin(np.linalg.norm(G[leaf_index] - G[inds], axis=1))]
        parent[leaf_index] = nn_index
        distances[leaf_index]=np.linalg.norm(G[leaf_index] - G[nn_index])
    
    # Create adjacency matrix
    T = np.zeros((num_G, num_G))
    for i in range(num_G):
        if parent[i] != -1:
            T[i, parent[i]] = T[parent[i], i] = distances[i]
    
    return T

def computeTreeCost(T):
    """
    This function takes an adjacency matrix T as a NumPy array, and computes
    the cost of the graph represented by T, which is the sum of all the
    values in the upper triangle of the matrix (excluding the diagonal).
    """
    n = len(T)
    cost = 0
    
    for i in range(n):
        for j in range(i+1, n):
            cost += T[i][j]
    return cost

import numpy as np

def perturb(inds, G, nm):
    """
    This function takes in a list of 5 integers `inds` in the range [0, 50),
    a NumPy array `G` of shape 50,2 containing coordinates of 50 points, and
    an integer `nm`. It randomly picks `nm` elements in `inds` and replaces
    them with any other valid integer such that each integer in `inds` is
    always unique. The function then returns a NumPy array of shape 5,2
    containing coordinates extracted from `G` corresponding to integers in
    `inds` in the same order.
    """
    
    # Randomly choose `nm` indices to replace
    indices_to_replace = np.random.choice(count_in, nm, replace=False)
    
    # Replace the chosen indices with new integers
    for i in indices_to_replace:
        new_int = np.random.choice([x for x in range(MAX) if x not in inds])
        inds[i] = new_int
    
    # Extract the corresponding coordinates from G
    coords = G[inds]
    
    return coords

min_w=sys.maxsize
last_improv=-1 #test
for i in range(n_iter):
    kmeans= KMeans(n_clusters=count_in,init=I,n_init=1).fit(G)
    inds=kcc_to_inds(kmeans.cluster_centers_,G)
    T=buildTree(inds,G) #T is Adj Matrix as 2d numpy array 
    w=computeTreeCost(T)
    #print(w)
    if(w<min_w):
        min_w=w
        np.savetxt("T.txt",T,fmt='%f')
        last_improv=i #test
    I=perturb(inds,G,1) #given inds and Graph coord, change 1 of the inds and return inds coord 

print("Tree Size: "+str(MAX))
print("Last Improvement: "+str(last_improv))
print("MST cost: "+str(min_w))