#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <vector>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cstdlib>
#include <sstream>
#include <cmath>

#define ERROR(FMT, ARG...) do {fprintf(stderr,"File=%s, Line=%d: \
        "FMT" \n",__FILE__, __LINE__,##ARG); exit(1);} while(0)

#define CL_SIZE 1
#define CLUSTER_ON 1 
#define FIRST_ID 0 //Assume the vertex ID starts from 0

class Vertex {
    public:
        int idx;
        double pr; 
        int in_deg;
        int out_deg;
        bool converged;

        std::vector<Vertex*> in_vertices;
        std::vector<Vertex*> out_vertices;

        bool updateDeg();

        Vertex(int _idx) {
            idx = _idx;
            converged = false;
        }

        ~Vertex(){
            std::vector<Vertex*>().swap(in_vertices);
            std::vector<Vertex*>().swap(out_vertices);
        }

};


class Edge{
    public:
        int src_idx;
        int dst_idx;
        double weight;

        Edge(int _src_idx, int _dst_idx, double _weight){
            src_idx = _src_idx;
            dst_idx = _dst_idx;
            weight = _weight;
        }
        ~Edge(){};
};


class Graph{
    public:
        std::vector<Vertex*> basic_vertices;
        std::vector<Vertex*> vertices; 
        std::vector<Edge*> edges; 
        int vertex_num;
        int edge_num;
        bool converged;
        int pair;

        Graph(){
            converged = false;
            pair = 0;
        }

        // Load edge description file and create graph
        bool createGraphFromFile(const std::string& file_name);
        bool clustering();
        int calPagerank(double alpha, int max_it, double convergence);
        bool printPR();
        bool dumpDot(const std::string& file_name);

        ~Graph(){
            std::vector<Edge*>::iterator eit;
            for(eit = edges.begin(); eit != edges.end(); eit++){
                delete (*eit);
            }
            std::vector<Edge*>().swap(edges);
            std::vector<Vertex*>::iterator vit;
            for(vit = basic_vertices.begin(); vit != basic_vertices.end(); vit++){
                delete (*vit);
            }
            std::vector<Vertex*>().swap(basic_vertices);
            std::vector<Vertex*>().swap(vertices);
        }

    private:
        void duplicateVertices();
        void loadEdgeData(const std::string& file_name, 
                std::vector<std::vector<double> >& data);
        int getVertexNum();
        int getEdgeNum();
        bool update(const std::vector<double>& edge_info);
        int getMaxVertexNum(const std::vector<std::vector<double> >& data);
        bool insertEmptyVertices(const int& vertex_num);
        bool isConverged();
        bool updatePR(const std::vector<double>& latest_pr);
        double getTwoVerSimilarity(const int& i, const int& j);
        int getNextIdx(const std::vector<int>& vertices_idx, 
                const std::vector<bool>& vertices_allocated);
        
        int oldIdx2NewIdx(int idx);
        int getMostSimilarVertices(
                const int& center_idx, 
                std::vector<bool>& vertices_allocated);

        void getMostSimilarVertices(
                const int& center_idx, 
                int to_be_clustered[CL_SIZE], 
                std::vector<bool>& vertices_allocated);
        double vertexPagerank(Vertex* vptr, const double& convergence, const double& alpha);
        void vertexPagerank(Vertex* vA, Vertex*vB, 
                double& curr_pr_of_A, double& curr_pr_of_B, 
                const double& convergence, 
                const double& alpha);
        void vertexPagerank(Vertex* vA, Vertex*vB, 
                Vertex* vC, Vertex* vD, 
                double& curr_pr_of_A, double& curr_pr_of_B,
                double& curr_pr_of_C, double& curr_pr_of_D, 
                const double& convergence, 
                const double& alpha);
        void vertexPagerank(Vertex* vA, Vertex*vB,
                Vertex* vC, Vertex* vD,
                Vertex* vE, Vertex* vF,
                Vertex* vG, Vertex* vH,
                double& curr_pr_of_A, double& curr_pr_of_B, 
                double& curr_pr_of_C, double& curr_pr_of_D,
                double& curr_pr_of_E, double& curr_pr_of_F,
                double& curr_pr_of_G, double& curr_pr_of_H, 
                const double& convergence,
                const double& alpha);


};

#endif
