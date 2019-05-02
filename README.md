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


