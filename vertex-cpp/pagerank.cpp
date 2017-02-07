#include "graph.h"

int main(int argc, char** argv){

    double alpha = 0.85;
    double convergence = 0.000001;
    int max_it = 5;

    Graph* gptr = new Graph();
    gptr->createGraphFromFile("../test/folkman.txt");
    gptr->calPagerank(alpha, max_it, convergence);
    gptr->printPR();
    gptr->dumpDot("./folkman.dot");

    return 0;

}
