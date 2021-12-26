#include "graph_boost.h"
Graph::Graph(void){
  aList = unG();
}
void Graph::printGraph(){
  graph_traits < unG>::vertex_iterator i, end;
  graph_traits < unG>::adjacency_iterator ai, a_end;
  property_map < unG, vertex_index_t >::type index_map = get(vertex_index, aList);
  int countE = 0;
  for (boost::tie(i, end) = vertices(aList); i != end; ++i) {
    std::cout << get(index_map, *i);
    boost::tie(ai, a_end) = adjacent_vertices(*i, aList);
    if (ai == a_end)
      std::cout << " has no children";
    else
      std::cout << " is the parent of ";
    for (; ai != a_end; ++ai) {
      countE++;
      std::cout <<get(index_map, *ai);
      if (boost::next(ai) != a_end){
        std::cout << ", ";
      }
    }
    std::cout << std::endl;
  }
}
void Graph::addEdge(int u, int v){
    add_edge(u,v,aList);
}
/*int main (int argc, char *argv[]) {
  unG aList(2);
  boost::add_edge(1,2,aList);
  graph_traits < unG>::vertex_iterator i, end;
  graph_traits < unG>::adjacency_iterator ai, a_end;
  property_map < unG, vertex_index_t >::type index_map = get(vertex_index, aList);
  int countV = 0; 
  int countE = 0;
  for (boost::tie(i, end) = vertices(aList); i != end; ++i) {
    countV++;
    std::cout << get(index_map, *i);
    boost::tie(ai, a_end) = adjacent_vertices(*i, aList);
    if (ai == a_end)
      std::cout << " has no children";
    else
      std::cout << " is the parent of ";
    for (; ai != a_end; ++ai) {
      std::cout <<get(index_map, *ai);
      if (boost::next(ai) != a_end){
        countE++;
        std::cout << ", ";
      }
    }
    std::cout << std::endl;
  }
  std::cout<< "****************"<<std::endl;
  Graph graph;
  graph.addEdge(0,1);
  //graph.addEdge(3,4);
  graph.printGraph();
  return 0;
}*/