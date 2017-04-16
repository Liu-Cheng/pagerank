#include "graph.h"

/* Vertex member functions */
bool Vertex::updateDeg(){
    in_deg = (int) in_vertices.size();
    out_deg = (int) out_vertices.size();
    return true;
}

/* Graph member functions */
int Graph::getNextIdx(const std::vector<int>& vertices_idx, 
                const std::vector<bool>& vertices_allocated){
    for(int i = 0; i < vertex_num; i++){
        int idx = vertices_idx[i];
        if(vertices_allocated[idx] == false){
            return idx;
        }
    }
    ERROR("Failed to find the unallocated vertex");
}

void Graph::vertexPagerank(Vertex* vA, Vertex*vB,
                Vertex* vC, Vertex* vD,
                Vertex* vE, Vertex* vF,
                Vertex* vG, Vertex* vH,
                double& curr_pr_of_A, double& curr_pr_of_B, 
                double& curr_pr_of_C, double& curr_pr_of_D,
                double& curr_pr_of_E, double& curr_pr_of_F,
                double& curr_pr_of_G, double& curr_pr_of_H, 
                const double& convergence,
                const double& alpha){
    std::vector<Vertex*>::iterator vit;
    for(vit = vA->in_vertices.begin(); vit != vA->in_vertices.end(); vit++){
        curr_pr_of_A += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    for(vit = vB->in_vertices.begin(); vit != vB->in_vertices.end(); vit++){
        curr_pr_of_B += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    for(vit = vC->in_vertices.begin(); vit != vC->in_vertices.end(); vit++){
        curr_pr_of_C += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    for(vit = vD->in_vertices.begin(); vit != vD->in_vertices.end(); vit++){
        curr_pr_of_D += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    for(vit = vE->in_vertices.begin(); vit != vE->in_vertices.end(); vit++){
        curr_pr_of_E += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    for(vit = vF->in_vertices.begin(); vit != vF->in_vertices.end(); vit++){
        curr_pr_of_F += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    for(vit = vG->in_vertices.begin(); vit != vG->in_vertices.end(); vit++){
        curr_pr_of_G += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    for(vit = vH->in_vertices.begin(); vit != vH->in_vertices.end(); vit++){
        curr_pr_of_H += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    double curr_eps_of_A = fabs(curr_pr_of_A - vA->pr);
    double curr_eps_of_B = fabs(curr_pr_of_B - vB->pr);
    double curr_eps_of_C = fabs(curr_pr_of_C - vC->pr);
    double curr_eps_of_D = fabs(curr_pr_of_D - vD->pr);
    double curr_eps_of_E = fabs(curr_pr_of_E - vE->pr);
    double curr_eps_of_F = fabs(curr_pr_of_F - vF->pr);
    double curr_eps_of_G = fabs(curr_pr_of_G - vG->pr);
    double curr_eps_of_H = fabs(curr_pr_of_H - vH->pr);

    if(curr_eps_of_A < convergence){
        vA->converged = true;
    }
    if(curr_eps_of_B < convergence){
        vB->converged = true;
    }
    if(curr_eps_of_C < convergence)
        vC->converged = true;

    if(curr_eps_of_D < convergence)
        vD->converged = true;

    if(curr_eps_of_E < convergence){
        vE->converged = true;
    }
    if(curr_eps_of_F < convergence){
        vF->converged = true;
    }
    if(curr_eps_of_G < convergence)
        vG->converged = true;

    if(curr_eps_of_H < convergence)
        vH->converged = true;

}

void Graph::vertexPagerank(Vertex* vA, Vertex*vB, 
        Vertex* vC, Vertex* vD, 
        double& curr_pr_of_A, double& curr_pr_of_B,
        double& curr_pr_of_C, double& curr_pr_of_D, 
        const double& convergence,
        const double& alpha){
    std::vector<Vertex*>::iterator vit;
    for(vit = vA->in_vertices.begin(); vit != vA->in_vertices.end(); vit++){
        curr_pr_of_A += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    for(vit = vB->in_vertices.begin(); vit != vB->in_vertices.end(); vit++){
        curr_pr_of_B += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    for(vit = vC->in_vertices.begin(); vit != vC->in_vertices.end(); vit++){
        curr_pr_of_C += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    for(vit = vD->in_vertices.begin(); vit != vD->in_vertices.end(); vit++){
        curr_pr_of_D += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    double curr_eps_of_A = fabs(curr_pr_of_A - vA->pr);
    double curr_eps_of_B = fabs(curr_pr_of_B - vB->pr);
    double curr_eps_of_C = fabs(curr_pr_of_C - vC->pr);
    double curr_eps_of_D = fabs(curr_pr_of_D - vD->pr);

    if(curr_eps_of_A < convergence){
        vA->converged = true;
    }
    if(curr_eps_of_B < convergence){
        vB->converged = true;
    }
    if(curr_eps_of_C < convergence)
        vC->converged = true;

    if(curr_eps_of_D < convergence)
        vD->converged = true;

}

void Graph::vertexPagerank(Vertex* vA, Vertex*vB, 
        double& curr_pr_of_A, double& curr_pr_of_B, 
        const double& convergence,
        const double& alpha){
    std::vector<Vertex*>::iterator vit;
    for(vit = vA->in_vertices.begin(); vit != vA->in_vertices.end(); vit++){
        curr_pr_of_A += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    for(vit = vB->in_vertices.begin(); vit != vB->in_vertices.end(); vit++){
        curr_pr_of_B += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    double curr_eps_of_A = fabs(curr_pr_of_A - vA->pr);
    double curr_eps_of_B = fabs(curr_pr_of_B - vB->pr);
    if(curr_eps_of_A < convergence){
        vA->converged = true;
    }
    if(curr_eps_of_B < convergence){
        vB->converged = true;
    }
}

double Graph::getTwoVerSimilarity(const int& i, const int& j){
    int in_deg1 = basic_vertices[i]->in_deg;
    int in_deg2 = basic_vertices[j]->in_deg;
    int avg_in_deg = (in_deg1 + in_deg2 + 1)/2;

    // Add some constrains
    if(in_deg1 == 0 || in_deg2 == 0 || avg_in_deg == 0){
        return 0;
    }

    int common_data_num = 0;
    std::vector<Vertex*>::iterator vit1;
    std::vector<Vertex*>::iterator vit2;
    for(vit1 = basic_vertices[i]->in_vertices.begin(); vit1 != basic_vertices[i]->in_vertices.end(); vit1++){
        for(vit2 = basic_vertices[j]->in_vertices.begin(); vit2 != basic_vertices[j]->in_vertices.end(); vit2++){
            if((*vit1)->idx == (*vit2)->idx){
                common_data_num++;
                break;
            }
        }
    }

    if(avg_in_deg > 0 && common_data_num > 10){
        pair++;
    }

    return (common_data_num * 1.0/avg_in_deg);
}

int Graph::getMostSimilarVertices(
        const int& center_idx, 
        std::vector<bool>& vertices_allocated){

    double max_similarity = -1;
    int idx = -1;
    for(int j = 0; j < vertex_num; j++){
        if(vertices_allocated[j] == true){
            continue;
        }

        double similarity = getTwoVerSimilarity(center_idx, j);
        if(max_similarity < similarity){
            max_similarity = similarity;
            idx = j;
        }
    }

    return idx;
}


void Graph::getMostSimilarVertices(
        const int& center_idx, 
        int to_be_clustered[CL_SIZE], 
        std::vector<bool>& vertices_allocated){
    if(vertices_allocated[center_idx] == true){
        std::cout << center_idx << " has been allocated." << std::endl;
    }
    vertices_allocated[center_idx] = true;
    to_be_clustered[0] = center_idx;

    for(int i = 1; i < CL_SIZE; i++){
        double max_similarity = -1;
        int idx = -1;

        for(int j = 0; j < vertex_num; j++){
            if(vertices_allocated[j] == true){
                continue;
            }

            double similarity = getTwoVerSimilarity(center_idx, j);
            if(max_similarity < similarity){
                max_similarity = similarity;
                idx = j;
            }
        }
        to_be_clustered[i] = idx;
        vertices_allocated[idx] = true;
    }
}

bool Graph::clustering(){

#if(CL_SIZE > 0 && CLUSTER_ON == 1)
  std::vector<int> vertices_in_deg(vertex_num, 0);
  std::vector<int> vertices_idx(vertex_num, 0);
  for(int i = 0; i < vertex_num; i++){
      vertices_idx[i] = i;
      vertices_in_deg[i] = basic_vertices[i]->in_deg;
  }

  // Reorder the vertex_id based on the vertex in degree
  for(int i = 0; i < vertex_num; i++){
      for(int j = i+1; j < vertex_num; j++){
          if(vertices_in_deg[i] < vertices_in_deg[j]){
              double deg_tmp;
              deg_tmp = vertices_in_deg[i];
              vertices_in_deg[i] = vertices_in_deg[j];
              vertices_in_deg[j] = deg_tmp;

              double idx_tmp;
              idx_tmp = vertices_idx[i];
              vertices_idx[i] = vertices_idx[j];
              vertices_idx[j] = idx_tmp;
          }
      }
  }


  // Clustering based on the similarity
  std::vector<bool> vertices_allocated(vertex_num, false);
#if (CL_SIZE > 1)
  int res = vertex_num%CL_SIZE;
  int cl_num = 0;
  int max_cl_num = vertex_num/CL_SIZE;
  while(cl_num < max_cl_num){
      int idx = getNextIdx(vertices_idx, vertices_allocated);
      int to_be_clustered[CL_SIZE];
      getMostSimilarVertices(idx, to_be_clustered, vertices_allocated);
      for(int i = 0; i < CL_SIZE; i++){
          Vertex* vptr = basic_vertices[to_be_clustered[i]];
          vertices.push_back(vptr);
      }
      cl_num++;
  }

  // The vertices that haven't been allocated will be processed in vertex granularity.
  for(int i = 0; i < res; i++){
      int idx = getNextIdx(vertices_idx, vertices_allocated); 
      vertices_allocated[idx] = true;
      vertices.push_back(basic_vertices[idx]);
  }
#else
  int idx = vertices_idx[0];
  for(int i = 0; i < vertex_num; i++){
      vertices_allocated[idx] = true;
      vertices.push_back(basic_vertices[idx]);
      idx = getMostSimilarVertices(idx, vertices_allocated);
  }
#endif

  // Duplicate graph storage based on the new order in vertices.
  duplicateVertices();

#else
  for(int i = 0; i < vertex_num; i++){
      vertices.push_back(basic_vertices[i]);
  }

#endif

    return true;

}

int Graph::oldIdx2NewIdx(int idx){
    for(int i = 0; i < vertex_num; i++){
        if(vertices[i]->idx == idx){
            return i;
        }
    }
    ERROR("Failed to find the vertx %d .", idx);
}

void Graph::duplicateVertices(){

    std::vector<Vertex*> shadow_vertices;
    shadow_vertices.resize(vertex_num);
    for(int i = 0; i < vertex_num; i++){
        Vertex* vptr = new Vertex(vertices[i]->idx);
        vptr->pr = vertices[i]->pr;
        vptr->converged = vertices[i]->converged;
        vptr->in_deg = vertices[i]->in_deg;
        vptr->out_deg = vertices[i]->out_deg;
        shadow_vertices[i] = vptr;
    }

    // Copy neighboring information
    for(int i = 0; i < vertex_num; i++){
        std::vector<Vertex*>::iterator vit;
        for(vit = vertices[i]->in_vertices.begin(); vit != vertices[i]->in_vertices.end(); vit++){
            int new_idx = oldIdx2NewIdx((*vit)->idx);
            shadow_vertices[i]->in_vertices.push_back(shadow_vertices[new_idx]);
        }
    }

    // Update data to vertices.
    for(int i = 0; i < vertex_num; i++){
        vertices[i] = shadow_vertices[i];
    }

}

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
    
    basic_vertices.resize(vertex_num);
    for(int i = 0; i < vertex_num; i++){
        Vertex* vptr = new Vertex(i);
        vptr->pr = 1.0/vertex_num;
        basic_vertices[i] = vptr;
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
    basic_vertices[src_idx]->out_vertices.push_back(basic_vertices[dst_idx]);
    basic_vertices[dst_idx]->in_vertices.push_back(basic_vertices[src_idx]);

    //Create edge and add it to the edge vector
    Edge* eptr = new Edge(src_idx, dst_idx, weight);
    edges.push_back(eptr);

    return true;
    
}

int Graph::getVertexNum(){
    return (int) basic_vertices.size();
}

int Graph::getEdgeNum(){
    return (int) edges.size();
}

int Graph::getMaxVertexNum(const std::vector<std::vector<double> >& data){

    std::vector<std::vector<double> >::const_iterator it1;
    std::vector<double>::const_iterator it2;
    double max_element = data[0][0];
    double min_element = data[0][0];
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
            if(*it2 < min_element)
                min_element = *it2;
        }
    }

    return (int) (max_element - min_element + 1);

}

void Graph::loadEdgeData(const std::string& file_name, 
        std::vector<std::vector<double> >& data){

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

    //Shift vertex ID by one such that vertex ID matches the vector index
#if (FIRST_ID == 1)
    int num = (int) data.size();
    for(int i = 0; i < num; i++){
        for(int j = 0; j < 2; j++){
            data[i][j]--;
        }
    }
#endif

}

bool Graph::createGraphFromFile(const std::string& file_name){

    std::vector<std::vector<double> > data;
    loadEdgeData(file_name, data);
    vertex_num = getMaxVertexNum(data);
    edge_num = (int)data.size();
    std::cout << "Total vertex num: " << vertex_num << std::endl;
    std::cout << "Total edge num: " << edge_num << std::endl;

    insertEmptyVertices(vertex_num);

    // Update edge related infomration
    std::vector<std::vector<double> >::iterator it1;
    for(it1 = data.begin(); it1 != data.end(); it1++){
        update(*it1);
    }

    // Update in_deg and out_deg
    std::vector<Vertex*>::iterator vit;
    for(vit = basic_vertices.begin(); vit != basic_vertices.end(); vit++){
        (*vit)->updateDeg();
    }

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

double Graph::vertexPagerank(Vertex* vptr, const double& convergence, const double& alpha){
    double curr_pr = (1 - alpha)/vertex_num;
    std::vector<Vertex*>::iterator vit;
    for(vit = vptr->in_vertices.begin(); vit != vptr->in_vertices.end(); vit++){
        curr_pr += alpha * (*vit)->pr * 1.0/(*vit)->out_deg;
    }

    double curr_eps = fabs(curr_pr - vptr->pr);
    if(curr_eps < convergence){
        vptr->converged = true;
    }

    return curr_pr;

}

int Graph::calPagerank(double alpha, int max_it, double convergence){

    int it_num = 0;
    // Ordered by the vertices[] and it doesn't match with the vertex idx
    std::vector<double> latest_pr; 
    latest_pr.resize(vertex_num);
    int res = vertex_num%CL_SIZE;

    while(it_num < max_it && !converged){
#if(CL_SIZE == 1)
        for(int i = 0; i < vertex_num; i++){
            Vertex* vptr = vertices[i];
            latest_pr[i] = vertexPagerank(vptr, convergence, alpha);
        }
#elif(CL_SIZE == 2)
        for(int i = 0; i < vertex_num - res; i = i + CL_SIZE){
            Vertex* vA = vertices[i];
            Vertex* vB = vertices[i+1];
            double curr_pr_of_A = (1 - alpha)/vertex_num;
            double curr_pr_of_B = curr_pr_of_A;
            vertexPagerank(vA, vB, curr_pr_of_A, curr_pr_of_B, convergence, alpha);
            latest_pr[i] = curr_pr_of_A;
            latest_pr[i+1] = curr_pr_of_B;
        }
#elif(CL_SIZE == 4)
        for(int i = 0; i < vertex_num - res; i = i + CL_SIZE){
            double curr_pr_of_A = (1 - alpha)/vertex_num;
            double curr_pr_of_B = curr_pr_of_A;
            double curr_pr_of_C = curr_pr_of_A;
            double curr_pr_of_D = curr_pr_of_A;
            Vertex* vA = vertices[i];
            Vertex* vB = vertices[i+1];
            Vertex* vC = vertices[i+2];
            Vertex* vD = vertices[i+3];
            vertexPagerank(vA, vB, vC, vD,
                    curr_pr_of_A, curr_pr_of_C, 
                    curr_pr_of_C, curr_pr_of_D, 
                    convergence, alpha);
            latest_pr[i] = curr_pr_of_A;
            latest_pr[i+1] = curr_pr_of_B;
            latest_pr[i+2] = curr_pr_of_C;
            latest_pr[i+3] = curr_pr_of_D;
        }

#elif(CL_SIZE == 8)
        for(int i = 0; i < vertex_num - res; i = i + CL_SIZE){
            double curr_pr_of_A = (1 - alpha)/vertex_num;
            double curr_pr_of_B = curr_pr_of_A;
            double curr_pr_of_C = curr_pr_of_A;
            double curr_pr_of_D = curr_pr_of_A;
            double curr_pr_of_E = curr_pr_of_A;
            double curr_pr_of_F = curr_pr_of_A;
            double curr_pr_of_G = curr_pr_of_A;
            double curr_pr_of_H = curr_pr_of_A;
            Vertex* vA = vertices[i];
            Vertex* vB = vertices[i+1];
            Vertex* vC = vertices[i+2];
            Vertex* vD = vertices[i+3];
            Vertex* vE = vertices[i+4];
            Vertex* vF = vertices[i+5];
            Vertex* vG = vertices[i+6];
            Vertex* vH = vertices[i+7];
            vertexPagerank(vA, vB, vC, vD, 
                    vE, vF, vG, vH,
                    curr_pr_of_A, curr_pr_of_B,
                    curr_pr_of_C, curr_pr_of_D,
                    curr_pr_of_E, curr_pr_of_F,
                    curr_pr_of_G, curr_pr_of_H,
                    convergence, alpha);
            latest_pr[i] = curr_pr_of_A;
            latest_pr[i+1] = curr_pr_of_B;
            latest_pr[i+2] = curr_pr_of_C;
            latest_pr[i+3] = curr_pr_of_D;
            latest_pr[i+4] = curr_pr_of_E;
            latest_pr[i+5] = curr_pr_of_F;
            latest_pr[i+6] = curr_pr_of_G;
            latest_pr[i+7] = curr_pr_of_H;
        }
#endif
        for(int i = vertex_num - res; i < vertex_num; i++){
            Vertex* vptr = vertices[i];
            latest_pr[i] = vertexPagerank(vptr, convergence, alpha);
        }

        // Update status after all the vertices' pagerank calculation
        it_num++;
        //converged = isConverged();
        updatePR(latest_pr);
    }

    return it_num;

};

bool Graph::printPR(){

    std::vector<double> pr(vertex_num, 0);
    std::vector<Vertex*>::iterator vit;
#if (CL_SIZE > 1 && CLUSTER_ON == 1)    
    //double sum = 0;
    for(vit = vertices.begin(); vit != vertices.end(); vit++){
        pr[(*vit)->idx] = (*vit)->pr;
        //sum += (*vit)->pr;
    }
#else
    for(vit = basic_vertices.begin(); vit != basic_vertices.end(); vit++){
        pr[(*vit)->idx] = (*vit)->pr;
    }
#endif
    //std::cout << "sum = " << sum << std::endl;
    for(int i = 0; i < vertex_num; i++){
        std::cout << pr[i] << std::endl;
    }

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


