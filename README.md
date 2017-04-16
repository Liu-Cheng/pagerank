# Clusteing the vertices in a graph for pagerank
This is a verification of a simple graph clustering idea.
Basically, we want to merge some of the vertices in a 
graph as super vertex and try to change the computing 
method accordingly for each super vertex. There are a 
few potential benefits for pagerank computing.  

1) Reduce function call of the page rank computing function.  
2) Potential data reuse within the clustered pagerank computing function.  
3) With compilation support, the clustering may be done transparently to the users.
Eventually, pagerank may be optimized automatically.  

At the same time, we also try to optimize the data layout on top 
of the clustering. To that end, we define a similarity for each vertex pair based on the 
overlapping of their incoming neighbors. Then we reoranize 
the graph layout based on the similarity information. 
However, according to the experiments, clusetring has negligible influence on the 
computing time and the data layout reorganization is beneficial. 

