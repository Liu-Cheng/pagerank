#include <ctime>
#include "graph.h"

int main(int argc, char** argv){

    double alpha = 0.85;
    double convergence = 0.000001;
    int max_it = 100;
    int it_num = 0;
    std::clock_t begin;
    std::clock_t end;
    double elapsed_sec;

    Graph* gptr = new Graph();
    //gptr->createGraphFromFile("../test/folkman.txt");
    //gptr->createGraphFromFile("../test/barabasi-100000.txt");
    //gptr->createGraphFromFile("./rmat.1mv.16me");
    //gptr->createGraphFromFile("./wikitalk.txt");
    gptr->createGraphFromFile("./rmat1.txt");

    begin = clock();
    gptr->clustering();
    end = clock();
    elapsed_sec = double(end - begin)/CLOCKS_PER_SEC;
    std::cout << "Clustering time: " << elapsed_sec << std::endl;

    begin = clock();
    it_num = gptr->calPagerank(alpha, max_it, convergence);
    end = clock();
    elapsed_sec = double(end - begin)/CLOCKS_PER_SEC;

    std::cout << "Pagerank rum time: " << elapsed_sec << std::endl;
    std::cout << "iteration number: " << it_num << std::endl;
    std::cout << "# of vertex pair with identical input: " << gptr->pair << std::endl;
    gptr->printPR();
    //gptr->dumpDot("./graph.dot");

    return 0;

}
