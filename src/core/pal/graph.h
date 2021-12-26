//self-defined graph representation, building adjacency_list using array of set
//TODO: sparse graph may fit set of set
//TODO: a map of index of labelLabel could save more memory
//TODO: unordered set?
#ifndef GRAPH_H
#define GRAPH_H
#include <assert.h>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <assert.h>
#include<vector>
#include<numeric>
#include<algorithm>
#include "priorityqueue.h"
#include "debugger.h"
#include "string.h"
using namespace std;
typedef unordered_set<int> edgeList;
typedef unordered_map<int,edgeList> adjecencyList; 
class Graph{
    public:
        Graph(int nblp, int all_nblp);
        void debugGraph();
        void printGraph();
        void addVertex(int u);
        void addEdge(int source, int target);
        void deleteVertex(int u);
        void deleteEdge(int source, int target);
        bool containEdge(int source, int target);
        bool containEdge_label(int source, int target);
        bool containVertex(int u);
        void outputDIMACS(string const &  fileName);
        void outputMetis(string const &  fileName);
    //******************vertex cover***********************
        void setPriorityQueue(pal::PriorityQueue * list);
        unordered_set<int> getVertexCover(int nblp, int all_nblp);
        void debugVertexCover(unordered_set<int>& vertexCover);
    //******************maxi. independent set************** 
        void getKAMIS(vector<int>& KAMIS);
        void debugMIS(vector<int>& vertexMIS);
        void readKAMIS(vector<int>& KAMIS,string const & fileName);
        void readWCNF(vector<int>& KAMIS,string const & fileName);
        void getMAXHS(vector<int>& KAMIS);
        adjecencyList adList; 
        int numV;
};
#endif