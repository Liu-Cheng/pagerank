# Clusteing the vertices in a graph for pagerank
This is a verification of a simple graph clustering idea.
Basically, we want to cluster some of the vertices in a 
graph as super vertices and try to change the computing 
method accordingly for the super vertices. Here are a 
few potential benefits of the clustering:  

1) It reduces function call of the page rank's vertex computing function 
assuming a vertex based graph processing framework.  
2) There is potential data reuse within the clustered pagerank computing function when combined.  
3) With compilation support, the clustering may be done transparently to the users.
Eventually, pagerank may be optimized automatically.  

At the same time, we also tried to optimize the data layout on top 
of the clustering. To that end, we define a similarity meric for each vertex pair based on the 
overlapping of their incoming neighbors. Then we reoranize 
the graph layout based on the similarity information. 
According to the experiments, clusetring has negligible influence on the 
computing time while the data layout reorganization is relativelly beneficial. 

