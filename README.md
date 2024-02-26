<h1>Single Objective Optimization</h1>
<h2>Leaf-Constrained Minimum Spanning Tree Problem</h2>
<div>
<h3>Problem Statement</h3>
<p>Given a connected, weighted and undirected graph G, spanning tree is an acyclic subgraph which connects all the vertices. Find a spanning tree on G such that there are k leaves and the total cost is minimum.</p>
</div>

<div>
<h3>Input Format</h3>

<b>Euclidean Instances : </b> lcmste\<test_case_id>

All the nodes are located within a unit square and uniformly distributed. Therefore each coordinate is a floating-point between 0 and 1.

<center><i> x, y &epsilon; (0,1) </i></center><br>

List of space-separated x and y coordinates

> \<x1> \<y1> <br>
> \<x2> \<y2> <br>
> \<x3> \<y3> <br>
> ... <br>
> \<x<sub>n</sub>> \<y<sub>n</sub>> <br>

<b>Randomized Instances : </b> lcmstr\<test_case_id>.txt

Adjacency Matrix Representation of the graph G wherein each entry A[ i ][ j ] represents the cost of the edge from node i to j.

</div>

<h3>Algorithm Used</h3>

<b>Kmeans Augmentation with Perturbation and Search Space Partitioning</b>

Time Complexity : <i>O(n<sup>2</sup>)</i>

<ol>
<li>Iteratively partition the search space such that the nodes are fairly distributed into regions.
<ol type='a'>
<li>The initial search space is divided into two halves vertically. Label the nodes as per locality.</li>
<li>Split the region with the greatest population along the shorter side and relabel the nodes accordingly.</li>
<li>Repeat step 2 till the number of regions(partitons) is equal to the number of internal nodes required.</li>
</ol>
</li>
<li>Obtain the centroids of each region by averaging the coordinates of the respective member nodes.</li> 
<li>Identify nodes nearest to the centroids as identified in previous step as Internal Nodes.</li>
<li>Pass these Internal Nodes as seeds to the Kmeans algorithm in the pipeline.</li>
<li>Connect the Internal Nodes, as obtained in step 4., by Prim's Algorithm.</li>
<li>Connect rest of the nodes to the nearest internal nodes thus completing the Spanning Tree as required.</li>
<li>Introduce perturbations by swapping one of the leaf nodes randomly with an internal node and repeating step 4-7. Record the best Minimum Spanning Tree obtained each iteration.</li>
</ol>

<h3>References</h3>

[1] Julstrom, B.A., 2004. Codings and operators in two genetic algorithms for the leaf-constrained minimum spanning tree problem. International Journal of Applied Mathematics and Computer Science, 14(3), pp.385-396.

[2] Prakash, V.P. and Patvardhan, C., 2020. Novel Heuristics for the Euclidean Leaf-Constrained Minimum Spanning Tree Problem. SN Computer Science, 1, pp.1-19.