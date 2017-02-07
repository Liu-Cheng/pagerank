#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <list>
#include <vector>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cstdlib>
#include <sstream>
#include <cmath>

#define ERROR(FMT, ARG...) do {fprintf(stderr,"File=%s, Line=%d: \
        "FMT" \n",__FILE__, __LINE__,##ARG); exit(1);} while(0)

//#define SIMPLE_TEST

class Vertex{
    public:
        int idx;
        double pr; // pagerank of the vertex
        std::list<Vertex*> in_vertices;
        std::list<Vertex*> out_vertices;
        bool converged;

        Vertex(int _idx){
            idx = _idx;
            converged = false;
        }

        int getInDeg();
        int getOutDeg();
        double calPagerank(const double& alpha, const int& vertex_num, const double& convergence);
};

double Vertex::calPagerank(const double& alpha, const int& vertex_num, const double& convergence){

    double curr_pr = (1 - alpha)/vertex_num;
    std::list<Vertex*>::iterator vit;
    for(vit = in_vertices.begin(); vit != in_vertices.end(); vit++){
        int out_deg = (*vit)->out_vertices.size();
        curr_pr += alpha * (*vit)->pr * 1.0/out_deg;
    }

    double curr_eps = fabs(curr_pr - pr);
    if(curr_eps < convergence){
        converged = true;
    }

    return curr_pr;

}

int Vertex::getInDeg(){
    return (int) in_vertices.size();
}

int Vertex::getOutDeg(){
    return (int) out_vertices.size();
}

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
};

class Graph{
    public:
        std::vector<Vertex*> vertices;
        std::vector<Edge*> edges; 
        int vertex_num;
        bool converged;

        Graph(){
            converged = false;
        }

        // Load edge description file and create graph
        bool createGraphFromFile(const std::string& file_name);
        int getVertexNum();
        int getEdgeNum();
        bool calPagerank(double alpha, int max_it, double convergence);
        bool printPR();
        bool dumpDot(const std::string& file_name);

    private:
        bool update(const std::vector<double>& edge_info);
        int getMaxVertexIdx(const std::vector<std::vector<double> >& data);
        bool insertEmptyVertices(const int& vertex_num);
        bool isConverged();
        bool updatePR(const std::vector<double>& latest_pr);
};

bool Graph::isConverged(){

    std::vector<Vertex*>::iterator vit;
    for(vit = vertices.begin(); vit != vertices.end(); vit++){
        if((*vit)->converged == false){
            return false;
        }
    }

    return true;
}

bool Graph::insertEmptyVertices(const int& vertex_num){
    
    for(int i = 0; i < vertex_num; i++){
        Vertex* vptr = new Vertex(i);
        vptr->pr = 1.0/vertex_num;
        vertices.push_back(vptr);
    }

    return true;
}

bool Graph::update(const std::vector<double>& edge_info){
    
    double weight;
    int src_idx, dst_idx;

    if(edge_info.size() == 2){
        src_idx = (int) edge_info[0];
        dst_idx = (int) edge_info[1];
        weight = 0;
    }
    else if(edge_info.size() == 3){
        src_idx = (int) edge_info[0];
        dst_idx = (int) edge_info[1];
        weight = edge_info[2];
    }
    else{
        ERROR("Unexpected edge infomation!\n");
    }

    //update vertex based on the edge information
    vertices[src_idx]->out_vertices.push_back(vertices[dst_idx]);
    vertices[dst_idx]->in_vertices.push_back(vertices[src_idx]);

    //Create edge and add it to the edge list
    Edge* eptr = new Edge(src_idx, dst_idx, weight);
    edges.push_back(eptr);

    return true;
    
}

int Graph::getVertexNum(){
    return (int) vertices.size();
}

int Graph::getEdgeNum(){
    return (int) edges.size();
}

int Graph::getMaxVertexIdx(const std::vector<std::vector<double> >& data){

    std::vector<std::vector<double> >::const_iterator it1;
    std::vector<double>::const_iterator it2;
    double max_element = data[0][0];
    for(it1 = data.begin(); it1 != data.end(); it1++){

        std::vector<double>::const_iterator last;
        if(it1->size() == 3) 
            last = it1->end() - 1;
        else if(it1->size() == 2)
            last = it1->end();
        else
            ERROR("Unexpected data size!");

        for(it2 = it1->begin(); it2 != last; it2++){
            if(*it2 > max_element)
                max_element = *it2;
        }
    }

    return (int) max_element;

}

bool Graph::createGraphFromFile(const std::string& file_name){

    std::vector<std::vector<double> > data;
    std::ifstream file_handle(file_name.c_str());
    if(!file_handle.is_open()){
        ERROR("Failed to open %s!", file_name.c_str());
    }

    std::string line;
    while(std::getline(file_handle, line)){
        std::istringstream iss(line);
        data.push_back(
                std::vector<double>(std::istream_iterator<double>(iss),
                    std::istream_iterator<double>()));
    }
    file_handle.close();

    vertex_num = getMaxVertexIdx(data) + 1;

    insertEmptyVertices(vertex_num);

    std::vector<std::vector<double> >::iterator it1;
    for(it1 = data.begin(); it1 != data.end(); it1++){
        update(*it1);
    }

#ifdef SIMPLE_TEST
    std::vector<double>::iterator it2;
    for(it1 = data.begin(); it1 != data.end(); it1++){
        for(it2 = it1->begin(); it2 != it1->end(); it2++){
            std::cout << *it2 << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "End of the edge print." << std::endl;
#endif

    return true;
}

bool Graph::updatePR(const std::vector<double>& latest_pr){
    std::vector<Vertex*>::iterator vit;
    int id = 0;
    for(vit = vertices.begin(); vit != vertices.end(); vit++){
        (*vit)->pr = latest_pr[id];
        id++;
    }

    return true;
}

bool Graph::calPagerank(double alpha, int max_it, double convergence){

    int it_num = 0;
    std::vector<Vertex*>::iterator vit;
    while(it_num < max_it && !converged){
        std::vector<double> latest_pr;
        for(vit = vertices.begin(); vit != vertices.end(); vit++){
            double tmp = (*vit)->calPagerank(alpha, vertex_num, convergence);
            latest_pr.push_back(tmp);
        }
        it_num++;
        converged = isConverged();
        updatePR(latest_pr);
    }

    return true;

};

bool Graph::printPR(){

    std::vector<Vertex*>::iterator vit;
    double sum = 0;
    for(vit = vertices.begin(); vit != vertices.end(); vit++){
        sum += (*vit)->pr;
        std::cout << "Pagerank of vertex[" << (*vit)->idx << "] = " << (*vit)->pr << std::endl;
    }
    std::cout << "sum = " << sum << std::endl;

    return true;
}

bool Graph::dumpDot(const std::string& file_name){

    std::ofstream file_handle(file_name.c_str());
    if(!file_handle.is_open()){
        ERROR("Failed to open %s!", file_name.c_str());
    }

    // created directed graph
    file_handle << "digraph G {" << std::endl;
    std::vector<Edge*>::iterator eit;
    for(eit = edges.begin(); eit != edges.end(); eit++){
        file_handle << "  V" << (*eit)->src_idx << " -> " 
            << "V" << (*eit)->dst_idx << ";" << std::endl;
    }
    file_handle << "}" << std::endl;

    file_handle.close();

    return true;
}

#endif
