#include "graph_map.h"
//bool gplDebugger = true;
Graph::Graph(int nblp, int all_nblp){
  numV = all_nblp;
}
void Graph::addVertex(int u){
    adList[u];
}
void Graph::addEdge(int u, int v){
    adList[u].insert(v);
    adList[v].insert(u);
};
void Graph:: deleteVertex(int u){
    adList.erase(u);
    /*for (const auto &p : adList) {
        adList[p.first].erase(u);
    }
    */

}
void Graph::deleteEdge(int u, int v){
    assert(adList.find(u) != adList.end());
    assert(adList.find(v) != adList.end());

    adList[u].erase(v);
    adList[v].erase(u);
};
bool Graph::containEdge(int source, int target){
    return (adList[source].count(target) != 0);
};
bool Graph::containEdge_label(int source, int target){
    return (adList[source].count(target) != 0);
};
bool Graph::containVertex(int u){
    return adList.find(u) != adList.end();
};
void Graph::debugGraph(){
        for (const auto &p : adList) {
            for(const auto &q : p.second){
                assert(containEdge(q,p.first));
            }
    }
};
void Graph::printGraph(){
    for (const auto &p : adList) {
            std::cout << p.first << " => ";
            for(const auto &q : p.second){
                        std::cout << q << " ";
            }
            cout<<'\n';
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
  int vSize= 0;
  int eSize= 0;
  vSize = adList.size();
  for (const auto &p : adList){
      eSize+=p.second.size();
  }
  outdata<<"p wcnf"<< vSize<<" "<< eSize/2<<endl;
  for (const auto &p : adList) {
    outdata<<endl;
    for(const auto &q : p.second){
        if(p.first < q){
              outdata<<top<<" "<<-p.first<<" "<< -q<<" "<< 0<< endl;
        }
    }
  }
  for (const auto &p : adList){
    outdata<<cost<<" "<<p.first<<" "<< 0<< endl;
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
  map<int,edgeList>::const_iterator it;
  for (const auto &p : adList){
      eSize+=p.second.size();
  }
  outdata<< numV<<" "<< eSize/2<<endl;
  int i = 1;
  for ( it = adList.begin(); it != adList.end(); it++)
  {
    for(; i < it->first; i++){
      outdata<< endl;      
    }
    for(const auto &q : it->second){
        outdata<<q<<" ";
    }
    i+=2;
    outdata<<endl;
  }
  outdata.close();
}
inline double getPriority(int degree){return degree;};

void Graph::setPriorityQueue(pal::PriorityQueue * list){
    int label, degree;
    for (const auto &p : adList) {
    label = p.first;
    degree = p.second.size();
    if(degree == 0) continue;
    list->insert(label, getPriority(degree));
  }
  if(gplDebugger){
    list->print();
  }
}
unordered_set<int> Graph:: getVertexCover(int nblp, int all_nblp){
  pal::PriorityQueue *list = nullptr;
  list = new pal::PriorityQueue( nblp, all_nblp, false );
  int label;
  int degree;
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
    for (const auto &p : adList) {
      label1 = p.first;
      if(p.second.size() == 0 || vertexCover.find(label1) != vertexCover.end()) continue;
      for(const auto &q : p.second){
        label2 = q;
        if(vertexCover.find(label2) == vertexCover.end()){
          cout<< "label 1: "<< label1 << endl;
          cout<< "label 2: "<< label2 << endl;
        }
        assert(vertexCover.find(label2) != vertexCover.end());
      }
    }
};
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
  debugVertexCover(Cover);
  for(int i =0; i < vertexMIS.size(); i++){
    for(int j = 0; j < i; j++){
      if(containEdge(vertexMIS[i],vertexMIS[j])){
        cout<< "i:"<< i<<endl;
        cout<< "j:"<< j<<endl;
      }
      assert(!containEdge(vertexMIS[i],vertexMIS[j]));
    }
  }
}
void Graph::getKAMIS(vector<int>& KAMIS){
  outputMetis("metis_map.txt");
  system("../redumis metis_map.txt --output=b_map.txt");
  readKAMIS(KAMIS, "b_map.txt");
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
  outputDIMACS("dimacs_map.txt");
  system("../maxhs dimacs_map.txt >b_map.txt");
  readWCNF(KAMIS, "b_map.txt");
}
/*
int main_metis (int argc, char *argv[]) {
  Graph graph(3,30);
  graph.addVertex(10);
  graph.addVertex(20);
  graph.addVertex(30);
  graph.addEdge(10,20);
  graph.addEdge(30,20);
  graph.addEdge(20,10);
  graph.printGraph();
  graph.debugGraph();
  graph.outputMetis("a_map.txt");
  vector<int> KAMIS;
  graph.getKAMIS(KAMIS);
  cout<< "KAMIS"<< endl;
  for(const auto elem: KAMIS){
    cout<< elem<< " ";
  }
  cout<<endl;
  return 0;
}*/
/*
int main(int argc, char *argv[]){
  Graph graph(3,300);
  graph.addVertex(100);
  graph.addVertex(200);
  graph.addVertex(300);
  graph.addEdge(100,200);
  graph.addEdge(300,200);
  graph.addEdge(200,100);
  graph.printGraph();
  graph.debugGraph();
  vector<int> KAMIS;
  graph.outputDIMACS("dimacs_map.txt");
  system("../../../maxhs dimacs_map.txt >b_map.txt");
  graph.readWCNF(KAMIS, "b_map.txt");
  cout<< "MAXHS"<< endl;
  for(const auto elem: KAMIS){
    cout<< elem<< " ";
  }
  cout<<endl;
  return 0;
}
*/