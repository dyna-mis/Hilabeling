#include "graph_set.h"
//bool gplDebugger = true;
Graph::Graph(int nblp, int all_nblp){
    numV = all_nblp+1;
    adList = new edgeList[numV];
}
void Graph::addVertex(int v){

};
void Graph::addEdge(int u, int v){
    assert(u< numV);
    assert(v < numV);
    adList[u].insert(v);
    adList[v].insert(u);
};
void Graph::deleteEdge(int u, int v){
    assert(u< numV);
    assert(v < numV);
    if(adList[u].erase(v)>0){
        adList[v].erase(u);
    }

};
bool Graph::containVertex(int u){
    return u< numV;
}
bool Graph::containEdge(int source, int target){
    return (adList[source].count(target) != 0);
};
bool Graph::containEdge_label(int source, int target){
    return (adList[source].count(target) != 0);
};
void Graph::debugGraph(){
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
void Graph::printGraph(){
    cout<< "&&&&&&&&&&&&&print Graph &&&&&&&&&&&&&&&"<< endl;
    for (int i = 0; i < numV; ++i){  
        cout << endl << i<< ": " ; 
        for (auto itr = adList[i].begin(); itr != adList[i].end(); ++itr){ 
            cout << *itr << " "; 
        }
        cout << endl; 
    } 
}
void Graph::outputDIMACS(string const &  fileName){
  ofstream outdata;
  int top = 2;
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
  outdata<<"p wcnf "<< numV<<" "<< eSize/2<<endl;
  for (int i =1 ; i < numV; i++) {
    for(const auto &q : adList[i]){
        if(i < q){
         outdata<<top<<" "<<-i<<" "<< -q<<" "<< 0<< endl;
        }
    }
  }
   for (int i =1 ; i < numV; i++) {
      outdata<<cost<<" "<<i<<" "<< 0<< endl;
   }
  outdata.close();
}
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
  for (int i =0 ; i < numV; i++) {
      eSize+=adList[i].size();
  }
  outdata<< numV<<" "<< eSize/2<<endl;
  for (int i =0; i < numV; i++) {
    for(const auto &q : adList[i]){
        outdata<<q<<" ";
    }
    outdata<<endl;
  }
  outdata.close();
}
inline double getPriority(int degree){return degree;};

void Graph::setPriorityQueue(pal::PriorityQueue * list){
    int degree;
    for (int i =1 ; i < numV; i++) {
    degree = adList[i].size();
    if(degree == 0) continue;
    list->insert(i, getPriority(degree));
  }
  if(gplDebugger){
    list->print();
  }
}
unordered_set<int> Graph:: getVertexCover(int nblp, int all_nblp){
  pal::PriorityQueue *list = nullptr;
  list = new pal::PriorityQueue( nblp, all_nblp, false );
  int label;
  unordered_set<int> vertexCover;
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
    return vertexCover;
  }
} 
void Graph::debugVertexCover(unordered_set<int>& vertexCover){
  cout<< "***********&&&&&&&&&&&***********"<< endl;
  for(const auto &p : vertexCover){
      cout<< p << endl;
  }
  int label1,label2;
    for (int i =1; i < numV; i++) {
      label1 = i;
      if(adList[i].size() == 0 || vertexCover.find(label1) != vertexCover.end()) continue;
      for(const auto &q : adList[i]){
        label2 = q;
        if(vertexCover.find(label2) == vertexCover.end()){
          cout<< "label 1 :"<< label1 << endl;
          cout<< "label 2 :"<< label2 << endl;

        }
        assert(vertexCover.find(label2) != vertexCover.end());
      }
    }
}
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
        KAMIS.push_back(elem);
  }
}
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
  cout<< "POINT A"<< endl;
  debugVertexCover(Cover);
  for(int i =0; i < vertexMIS.size(); i++){
    for(int j = 0; j < i; j++){
      if(containEdge(vertexMIS[i],vertexMIS[j])){
        cout<< "i: "<< vertexMIS[i]<< endl;
        cout<< "j: "<< vertexMIS[j]<< endl;

      }
      assert(!containEdge(vertexMIS[i],vertexMIS[j]));
    }
  }
}
void Graph::getKAMIS(vector<int>& KAMIS){
  outputMetis("metis_set.txt");
  system("../redumis metis_set.txt --output=b_set.txt");
  readKAMIS(KAMIS, "b_set.txt");
}
void Graph::readWCNF(vector<int>& KAMIS,string const & fileName){
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
       KAMIS.push_back(i);
      }
      token = strtok(NULL, s);
    }
  if(gplDebugger){
    debugMIS(KAMIS);
  }
}
void Graph::getMAXHS(vector<int>& KAMIS){
  outputDIMACS("dimacs_set.txt");
  system("../maxhs dimacs_set.txt >b_set.txt");
  readWCNF(KAMIS, "b_set.txt");
}
/*int main(int argc, char *argv[]) {
  Graph graph(3,3);
  graph.addEdge(1,2);
  graph.addEdge(3,2);
  graph.addEdge(2,1);
  graph.printGraph();
  graph.debugGraph();
  vector<int> KAMIS;
  graph.getKAMIS(KAMIS);
  cout<< "KAMIS"<< endl;
  for(const auto elem: KAMIS){
    cout<< elem<< " ";
  }
  cout<<endl;
  return 0;
}
int main(int argc, char *argv[]){
  Graph graph(3,300);
  graph.addVertex(100);
  graph.addVertex(200);
  graph.addVertex(300);
  graph.addEdge(100,200);
  graph.addEdge(300,200);
  graph.addEdge(200,100);
  graph.printGraph();
  cout<< "POINT A"<<endl;
  graph.debugGraph();
  cout<< "POINT B"<<endl;
  vector<int> KAMIS;
  graph.outputDIMACS("dimacs_set.txt");
  cout<< "POINT C"<<endl;
  system("../../../maxhs dimacs_set.txt >b.txt");
    cout<< "POINT D"<<endl;
  graph.readWCNF(KAMIS, "b.txt");
  cout<< "MAXHS"<< endl;
  for(const auto elem: KAMIS){
    cout<< elem<< " ";
  }
  cout<<endl;
  return 0;
}
*/