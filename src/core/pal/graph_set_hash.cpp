#include "graph_set_hash.h"
#include <float.h>
//bool gplDebugger = true;
double e;
double top;
// nblp: number of vertices
//all_nblp: maxID of vertices (no use in this representaion)
Graph::Graph(int nblp, int all_nblp){
    numV = nblp+1;
    weights.push_back(-1);
    adList = new edgeList[numV];
    table = new lookupTable(nblp);
}
// the solution_prev stores the avaliable vertexIDs in this round, which are in the previous solution.
//l1: labelID of the vertex
// l1 can be cached after it is added in the lookUptable (by addVertex). 
//TODO: combinae addCache and addVertex together
void Graph::addCache(int l1){
  int v1 = table->lookUpVID(l1);
  solution_prev.insert(v1);
};
// debug code
void printWeights(vector<double> weights){
  cout<< "the weights size is " << weights.size()<< endl;
    for (const auto &p : weights){
      cout<< p << " ";
    }
}
// add a new vertex in the graph
// label: labelID
// weight: initial label weight (priority)
void Graph::addVertex(int label, double weight){
  if(gplPrinter){
    printWeights(weights);
  }
  int v1 = table->insert(label);
  weights.push_back(weight);
  assert(weights[v1] == weight);
};
// To give higher weights to labels in previous solution,
// their weights will be added by e.
// Set e as half of the minimum weight 
void Graph:: setE(){
  double min = DBL_MAX;
  for (int i = 1; i < weights.size(); i++){
    if(weights[i]< min) {
      min = weights[i];
    }
  }
  e = min/2;
  assert(e > 0);
}
// TOP is the weight for hard constrains (like hard clause in maxHS)
// Set top as the sum of all weights +1;
// Be careful: set this after the weights are finally adjusted. 
void Graph:: setTOP(){
  int sum = 0;
  for (int i = 1; i < weights.size(); i++){
    sum += weights[i];
  }
  top = sum +1;
}

// give extra weight e to vertex v.
inline void Graph::increaseWeight(int v){
  weights[v]+= e;
}
// increase weights for previous solution
void Graph::adjustWeights(){
  setE();
  for (const auto &p : solution_prev) {
    increaseWeight(p);
  }
  setTOP();
  // for maxHS (top > sum(weights))
};
// add edge between label l1 and label l2
// Be careful: these two vertices are already in the garph (added by addVertex)
void Graph::addEdge(int l1, int l2){
    int v1 = table->lookUpVID(l1);
    int v2 = table->lookUpVID(l2);
    assert(v1 != -1);
    assert(v2 != -1);
    assert(v1< numV);
    assert(v2 < numV);
    adList[v1].insert(v2);
    adList[v2].insert(v1);
};
// Never use yet
// may use after TODO: modification inside pal directly
void Graph::deleteEdge(int l1, int l2){
    int v1 = table->lookUpVID(l1);
    int v2 = table->lookUpVID(l2);
    assert(v1 != -1);
    assert(v2 != -1);
    assert(v1< numV);
    assert(v2 < numV);
    if(adList[v1].erase(v2)>0){
        adList[v2].erase(v1);
    }

};
// check if one vertex for this label in the graph 
//label: labelID
bool Graph::containVertex(int label){
    int u = table->lookUpVID(label);
    return (u> 0 && u<= numV);
}
// check if one edge between two labels exists
//l1: labelID 
//l2: labelID
bool Graph::containEdge_label(int l1, int l2){
    int v1 = table->lookUpVID(l1);
    int v2 = table->lookUpVID(l2);
    assert(v1 != -1);
    assert(v2 != -1);
    assert(v1< numV);
    assert(v2 < numV);
    return containEdge(v1, v2);
};
// check if one edge between two vertices exists
// v1: vertexID 
// v2: vertexID
bool Graph::containEdge(int v1, int v2){
    return (adList[v1].count(v2) != 0);
};
// check vertexID<< numV
// check if edge bidirected (in an undirected graph)
void Graph::debugGraph(){
    table->debug();
    edgeList::iterator it;
    for(int i =0; i < numV;i++){
        it= adList[i].begin();
        while (it != adList[i].end()){
            //undirected graph (u,v)->(v,u)
            assert(containEdge(*it,i));
            //index < numV
            assert(*it <= numV);
            it++;
        }
    }
};
// print adjecency list
void Graph::printGraph(){
  cout<< "&&&&&&&&&&&&&print Graph &&&&&&&&&&&&&&&"<< endl;
    table->print();
    for (int i = 1; i < numV; ++i){  
        cout << endl << i<< ": " ; 
        for (auto itr = adList[i].begin(); itr != adList[i].end(); ++itr){ 
            cout << *itr << " "; 
        }
        cout << endl; 
    } 
}
// output Weighted dimas format for maxHS
void Graph::outputDIMACS(string const &  fileName){
  ofstream outdata;
  // unweighted first, using  2 as top value for hard clause, 1 for soft clause
  int cost = 1;
  outdata.open (fileName.c_str());
  if(!outdata){
    // file couldn't be opened
    cerr << "Error: file could not be opened" << endl;
    exit(1);  
  }
  outdata << "c "<<fileName<<endl;
  int eSize= 0;
  for (int i =1 ; i < numV; i++) {
      eSize+=adList[i].size();
  }
  // unweighted first, using  2 as top value for hard clause
  outdata<<"p wcnf "<< numV-1<<" "<< numV-1+ eSize/2<<" "<< top<< endl;
  for (int i =1 ; i < numV; i++) {
    for(const auto &q : adList[i]){
      if(q < i){
        outdata<<top<<" "<<-i<<" "<< -q<<" "<< 0<< endl;
      }
    }
  }
  for (int i =1 ; i < numV; i++){
    // use weight
    cost = weights[i];
    outdata<<cost<<" "<<i<<" "<< 0<< endl;
  }
  outdata.close();
}
// output the graph in metis format for kamis
// Be careful: kamis needs ordered edgelists.
void Graph::outputMetis(string const & fileName){
  ofstream outdata;
  outdata.open (fileName.c_str());
  if(!outdata){
    // file couldn't be opened
    cerr << "Error: file could not be opened" << endl;
    exit(1);  
  }
  outdata << "% "<<fileName<<endl;
  int eSize= 0;
  for (int i =1 ; i < numV; i++) {
      eSize+=adList[i].size();
  }
  outdata<< numV-1<<" "<< eSize/2<<endl;
  for (int i =1 ; i < numV; i++) {
    for(const auto &q : adList[i]){
        outdata<<q<<" ";
    }
    outdata<<endl;
  }
  outdata.close();
}
// TODO: add weights in the formular
// It dosent work if divided by (weight+1). 
//may beed to check the inplementaion of priority queue"
inline double getPriority_weighted(int degree, double weight){
  return (double)degree/(weight+1.0);
};
inline double getPriority(int degree){
  return (double)degree;
};
void Graph::setPriorityQueue_weighted(pal::PriorityQueue& list, vector<int>& degrees){
  int degree;
  for (int i =1 ; i < numV; i++) {
    degree = adList[i].size();
    degrees.push_back(degree);
    if(degree == 0) continue;
    list.insert(i, getPriority_weighted(degree, weights[i]));
  }
  if(gplPrinter){
    list.print();
  }
}
void Graph::setPriorityQueue(pal::PriorityQueue * list){
  int degree;
  for (int i =1 ; i < numV; i++) {
    degree = adList[i].size();
    if(degree == 0) continue;
    list->insert(i, getPriority(degree));
  }
  if(gplPrinter){
    list->print();
  }
}
// get vertex cover ( priority of vertex u is  deg(u)/w(u)+1 ) 
unordered_set<int> Graph:: getVertexCover_weighted(int nblp, int all_nblp){
  vector<int> degrees; 
  degrees.push_back(-1);
  // false: higher value, higher priority
  pal::PriorityQueue list( nblp, nblp, false);
  int label;
  unordered_set<int> vertexCover;
  unordered_set<int> labelCover;
  setPriorityQueue_weighted(list,degrees);
  if(gplPrinter){
    for(int i = 0; i < degrees.size(); i++){
      cout<< i << " " << degrees[i]<< endl;
    }
  }
  //
  unordered_set<int> covered;
  unsigned int size = list.getSize();
  while ( covered.size() != size){
      label = list.getBest(); 
      vertexCover.insert(label);
      covered.insert(label);
      if(covered.size() == size) break;
      for(const auto &p: adList[label]){
        if(degrees[p] == 1){
          covered.insert(p);
          list.remove(p);
        }
        else{
          degrees[p]--;
          list.changePriority(p, getPriority_weighted(degrees[p], weights[p]));
        }
        if(covered.size() == size) goto check_pointer;
      }
  }
  check_pointer:{
    if(gplDebugger){
      debugVertexCover(vertexCover);
    }
  }
  for (const auto& elem: vertexCover){ 
    labelCover.insert(table->lookUpLID(elem));
  }
    return labelCover;
} 
// debug code 
void Graph::debugDegree(int vertex, vector<int> degrees,   unordered_set<int>& vertexCover){
  int label2;
  int r = 0;
  for(const auto &q : adList[vertex]){
      label2 = q;
      if(vertexCover.count(q)>0){
        r++;
      }
  }
  if(degrees[vertex] + r != adList[vertex].size()){
    cout<< endl<<"***********And the degrees until now***********"<< endl;
    for(int i = 0; i < degrees.size(); i++){
      cout<<i<< " "<<  degrees[i] << endl;
    }
    cout<< "***********And the vertexCover until now***********"<< endl;
    int i = 0;
    for(const auto &p : vertexCover){
      cout<<i<< " "<<  p << endl;
      i++;
    }
    cout<< "vertex: "<< vertex<< endl;
    cout<< "r: "<< r << endl;
    cout<< "adList[vertex].size(): "<< adList[vertex].size()<< endl;
    cout<< "degrees[vertex]: "<< degrees[vertex]<< endl;
  }
  assert(degrees[vertex] + r == adList[vertex].size());
} 
// get unweighted vertex cover ( priority of vertex u is  deg(u)) 
// use this in mis(), just change mis() function in problem.cpp (2814)
unordered_set<int> Graph:: getVertexCover(int nblp, int all_nblp){
  pal::PriorityQueue *list = nullptr;
  // true: sort by growth
  //list = new pal::PriorityQueue( nblp, nblp, true );
  list = new pal::PriorityQueue( nblp, nblp, false);
  int label;
  unordered_set<int> vertexCover;
  unordered_set<int> labelCover;
  setPriorityQueue(list);
  unordered_set<int> covered;
  unsigned int size = list->getSize();
  while ( covered.size() != size){
      label = list->getBest(); 
      vertexCover.insert(label);
      covered.insert(label);
      if(covered.size() == size) break;
      for(const auto &p: adList[label]){
        if(list->decreaseKey_remove(p,0.0)){
          covered.insert(p);
          if(covered.size() == size) goto check_pointer;
        }
      }
  }
  check_pointer:{
    if(gplDebugger){
      debugVertexCover(vertexCover);
    }
  }
  for (const auto& elem: vertexCover){ 
    labelCover.insert(table->lookUpLID(elem));
  }
    return labelCover;
} 
// bebug code
void Graph:: debugCover(int vertex, unordered_set<int>& vertexCover, vector<int>& degrees){
    int i,label2;
    i = vertex;
    if(adList[i].size() == 0 || vertexCover.find(i) != vertexCover.end()) return;
    for(const auto &q : adList[i]){
      label2 = q;
      if(vertexCover.find(label2) == vertexCover.end()){
        cout<< endl<< "label 1 "<< i<< endl;
        cout<< "label 2 "<< label2<< endl;
        cout<< endl<<"***********And the degrees until now***********"<< endl;
        for(int i = 0; i < degrees.size(); i++){
          cout<<i<< " "<<  degrees[i] << endl;
        }
          cout<< "***********And the vertexCover until now***********"<< endl;
          int i = 0;
          for(const auto &p : vertexCover){
              cout<<i<< " "<<  p << endl;
              i++;
          }
      }
        assert(vertexCover.find(label2) != vertexCover.end());
    }
}
// debug code
void Graph::debugVertexCover(unordered_set<int>& vertexCover){
  cout<< "***********&&&&&&&&&&& debugVertexCover***********"<< endl;
  for(const auto &p : vertexCover){
      cout<< p <<" "<<endl;
  }
  int label1,label2;
    for (int i =1; i < numV; i++) {
      label1 = i;
      if(adList[i].size() == 0 || vertexCover.find(label1) != vertexCover.end()) continue;
      for(const auto &q : adList[i]){
        label2 = q;
        if(vertexCover.find(label2) == vertexCover.end()){
          cout<< "label 1 "<< label1<< endl;
          cout<< "label 2 "<< label2<< endl;
        }
        assert(vertexCover.find(label2) != vertexCover.end());
      }
    }
}
// get MIS from kamis file
// In kamis, one line for each node. 1 for taken, 0 for not taken 
void Graph::readKAMIS(vector<int>& KAMIS,string const & fileName){
  vector<int> vertexMIS;
  ifstream indata;
  indata.open(fileName.c_str());
  if(!indata){
  // file couldn't be opened
  cerr << "Error: file could not be opened" << endl;
  exit(1);  
  }
  int line = 1;
  int in;
  while (indata>>in){
    if(in ==1){
      vertexMIS.push_back(line);
    }
    line++;
  }
  if(gplDebugger){
    debugMIS(vertexMIS);
  }
  for(const auto elem: vertexMIS){
        KAMIS.push_back(table->lookUpLID(elem));
  }
}
// read outputfile from  maxHS-like SAT programm
void Graph::readWCNF(vector<int>& KAMIS,string const & fileName){
  vector<int> vertexMIS;
  ifstream indata;
  indata.open(fileName.c_str());
  if(!indata){
  // file couldn't be opened
  cerr << "Error: file could not be opened" << endl;
  exit(1);  
  }
  string line;
  char head;
  while (getline(indata, line)){
    if(line.empty()) continue;
		head =line.at(0);
		if(head == 'v'){
			break;
		}
  }
  int i;
  char* str = strdup(line.c_str());
  const char s[2] = " ";
  char* token = strtok(str, s);
  token = strtok(NULL, s);
    while(token != NULL){
      i = atoi(token);
      if(i>0){
        vertexMIS.push_back(i);
      }
      token = strtok(NULL, s);
    }
  if(gplDebugger){
    debugMIS(vertexMIS);
  }
  for(const auto elem: vertexMIS){
        KAMIS.push_back(table->lookUpLID(elem));
  }
}
// debug code 
// to check if the set is independent
// to check if its complement a vertex cover (maximal)
void Graph::debugMIS(vector<int>& vertexMIS){
  cout<< "***********&&&&&&&&&&&debugMIS***********"<< endl;
  for(const auto & elem : vertexMIS){
    cout<< elem << " ";
  }
  cout<< endl;
  std::vector<int> f(numV);
  std::iota(f.begin(), f.end(), 1);
  vector<int> cover;
  std::set_difference(f.begin(),f.end(),vertexMIS.begin(), vertexMIS.end(), back_inserter(cover));
  std::unordered_set<int> Cover(cover.begin(), cover.end());
  debugVertexCover(Cover);
  for(int i =0; i < vertexMIS.size(); i++){
    for(int j = 0; j < i; j++){
      assert(!containEdge(vertexMIS[i],vertexMIS[j]));
    }
  }
}
// get MIS by kamis
// KAMIS contains the labelID in original problem.
void Graph::getKAMIS(vector<int>& KAMIS){
  outputMetis("metis_set_map.txt");
  system("../redumis metis_set_map.txt --output=b_set_map.txt");
  readKAMIS(KAMIS, "b_set_map.txt");
}
// get MIS by maxHS
// KAMIS contains the labelID in original problem.
void Graph::getMAXHS(vector<int>& KAMIS){
  outputDIMACS("dimacs_set_map.txt");
  system("../maxhs dimacs_set_map.txt >b_set_map.txt");
  readWCNF(KAMIS, "b_set_map.txt");
}
// debug code
/*int main_metis (int argc, char *argv[]) {
  Graph graph(3,3);
  graph.addVertex(100);
  graph.addVertex(200);
  graph.addVertex(300);
  graph.addEdge(100,200);
  graph.addEdge(300,200);
  graph.addEdge(200,100);
  graph.printGraph();
  graph.debugGraph();
  graph.outputMetis("metis_set_map.txt");
  vector<int> KAMIS;
  graph.getKAMIS(KAMIS);
  cout<< "KAMIS"<< endl;
  for(const auto elem: KAMIS){
    cout<< elem<< " ";
  }
  cout<<endl;
  return 0;
}
*/
/*
// debug code
int main(int argc, char *argv[]){
  Graph graph(3,3);
  graph.addVertex(100);
  graph.addVertex(200);
  graph.addVertex(300);
  graph.addEdge(100,200);
  graph.addEdge(300,200);
  graph.addEdge(200,100);
  graph.printGraph();
  graph.debugGraph();
  graph.outputDIMACS("dimacs_set_map.txt");
  vector<int> KAMIS;
  graph.getMAXHS(KAMIS);
  cout<< "KAMIS"<< endl;
  for(const auto elem: KAMIS){
    cout<< elem<< " ";
  }
  cout<<endl;
  return 0;
}
*/