#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
using namespace boost;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property> unG;
class Graph{
    public:
        Graph();
        void printGraph();
        bool checkGraph(){return true;};
        void addEdge(int u, int v);
    private:
         unG aList;

};
//TODO: using labelposition's id may cause a lot of nullptr (try build structure for vertex descriptor)