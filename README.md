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
<li>The search space is divided into two halves vertically. Labelling the left and the right region as 0 and 1 respectively.</li>
<li>Split the side with the greatest population across the longer side.</li>
<li>Repeat step 2 till the number of regions is equal to the number of internal nodes required.</li>
<li>Obtain the centroids of each region by averaging the coordinates of the respective member nodes.</li> 
<li>Identify nodes nearest to respective centroids identified in step 4 as Internal Nodes.</li>
<li>Pass these Internal Nodes as seeds to the Kmeans algorithm in the pipeline.</li>
<li>Connect the Internal Nodes by Prim's Algorithm.</li>
<li>Connect rest of the nodes to the nearest internal nodes each thus completing the Spanning Tree as required.</li>
<li>Introduce perturbations by swapping one of the leaf nodes in the search space with an internal node and repeating step 6-8. Record the best Minimum Spanning Tree obtained so far with each iteration.</li>
</ol>

<h3> Implementation Details </h3>

<b>Data Structures</b>

<ul>
<li>Visited Arrays</li>
<li>Structures</li>
</ul>