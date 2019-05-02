# AMG-Graph-Hierarchy
Algebraic multigrid graph coarsening

This repository contains several versions of algebraic multigrid coarsening for graphs used for solving such combinatorial optimization problems on graphs as the minimum linear arrangement, 2-sum, bandwidth, k-partitioning, network compression, sparsification, optimal response to attacks, and more.

If you use this code, please cite

Dorit Ron, Ilya Safro, Achi Brandt, "Relaxation-based coarsening and multiscale graph organization", SIAM Multiscale Modeling and Simulations, Vol. 9, No. 1, pp. 407-423, 2011

Jie Chen, Ilya Safro, "Algebraic distance on graphs", SIAM Journal on Scientific Computing, Vol. 33, No. 6, pp. 3468-3490, 2011

The code contains many redundant procedures that were used to optimize different objectives on graphs and definitely requires some cleaning. The dependencies include only LEDA C++ library free edition (see http://www.algorithmic-solutions.com/index.php/products/leda-free-edition). The simplest way to use it is running

> make

> ./buildAMGhierarchy ./[graph_file] ./[parameter_file] [random_seed]

The resulting AMG hierarchhy with coarsened graphs will be stored in files level[X].dat
The input is accepted in several formats that can be found in readprint.cpp but the most straightforward is just a list of weighted edges with a single header row for the numbers of nodes and edges in a graph (extension .edges is mandatory). The node numbering is started with 1.
Examples of parameters are available in coarsening_param1 (aggregates full nodes) and coarsening_param2 (splits fine nodes in at most two parts)

Input graph example: can_1072.edges

1072 5686   <-- number of nodes [space] number of edges

1 2 1.0   <-- there is an edge from node 1 to node 2 with weight 1.0

1 46 1.0   <-- there is an edge from node 1 to node 46 with weight 1.0

...

Output file with a level k in AMG hierarchy:
line 1: number of nodes [space] number of edges
lines i in [2..end]: information about node i-1 and its adjacency
line i structure:
1. Starts with either 'c' (for a seed or coarse node that will be a center of some aggregate at level k+1) or 'f' (a fine node, i.e., a node that is fully or partially included in one or more of its 'c' neighbors)
2. If there is a 'c'-node in line i, then the next integer is the id of its corresponsing aggregate at level k+1 followed by node weight (aka volume). If there is an 'f'-node in line i then next integer is immediately the node weight.
3. Adjacency information that follows 'c' or 'f' node info:
If node i-1 is a 'c' node, then the adjacency list includes pairs of the current level neighbor id and connecting edge weight. If node i-1 is an 'f' node, then the adjacency list includes either 
-- pairs of the current level neighbor id and connecting edge weight if the neighbor is also 'f' node or
-- triples of the current level neighbor id, connecting edge weight, and AMG interpolation weight denoted by 'iw:' if the neighbor is a 'c' node. Note that some interpolation weights could be zero which means that this 'f' node does not participate in the aggregate of this 'c' neighbor.

Example:
261 1493 <-- number of nodes [space] number of edges

f 5.78261 65 10.5501 iw:0 121 9.15459 iw:0.523749 153 0.166667 46 0.26087 18 0.623819 iw:0 87 0.130435 9 8.32437 iw:0.476251 66 0.375 iw:0 77 0.564 516 iw:0 82 0.125 iw:0 130 5.33333 <-- node 1 is an 'f' node with weight 5.78261; it is connected to node 65 with edge weight 10.5501 but interpolation weight is 0; it is connected to node 121 with edge weight 9.15459 with  interpolation weight 0.523749 ...

c 1 8.11429 5 5.69048 125 7.04 170 6.42295 106 2.26642 216 2.21158 220 5.47048 246 7.75 195 4.22619 173 1.5 60 2.14 92 3.23743 109 0.11 113 0.75 20 2 7.32019 215 6.64905 <-- node 2 is a 'c' node whose aggregate's id at the next level is 1 and current level weight is 8.11429; its first neighbor is 5, the edge between them is 5.69048; its second neighbor is 125 with edge weight 7.04

If you have any questions feel free to email at isafro@clemson.edu

Explanation of basic AMG terminology on graphs is available in 
Safro, Ron, Brandt, "Graph Minimum Linear Arrangement by Multilevel Weighted Edge Contractions", Journal of Algorithms, vol. 60/1, pp. 24-41, 2006


