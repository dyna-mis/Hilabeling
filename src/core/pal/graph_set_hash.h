// Self-defined graph representation.
// the labelIDs in QGIS are discrete in a large range.
// To build a small graph on it, I compare different ways.
// The current one in use first maps the labelID to vertexID in (much) smaller range, 
// Then uses an array of sets at the convenience of ordered adjecency list for weighted DIMACS.
// It takes maybe more time for conversion and lookup between vertexID and labelID.
// Considering the extern algorithms we choose may be sensitive to the size of verices, we use this way.
// TODO: A more efficient representaion is expected

#ifndef GRAPH_SET_HASH
#define GRAPH_SET_HASH
#include <assert.h>
#include <set>
#include <iostream>
#include<fstream>
#include "debugger.h"
#include <unordered_set>
#include "priorityqueue.h"
#include "gpl_datastructure.h"
#include <numeric>
#include <algorithm>
#include <string.h>
#include <stdlib.h>
using namespace std;
typedef set<int> edgeList;
class Graph{
    public:
        Graph(int nblp, int all_nblp);
        void debugGraph();
        void printGraph();
        void addVertex(int u, double wight);
        void addEdge(int source, int target);
        void deleteEdge(int source, int target);
        bool containVertex(int u);
        bool containEdge(int source, int target);
        bool containEdge_label(int l1, int l2);
        void setPriorityQueue_weighted(pal::PriorityQueue& list, vector<int>& degrees);
        void setPriorityQueue(pal::PriorityQueue * list);
        void outputDIMACS(string const &  fileName);
        void outputMetis(string const & fileName);
        unordered_set<int> getVertexCover_weighted(int nblp, int all_nblp);
        unordered_set<int> getVertexCover(int nblp, int all_nblp);
        void debugCover(int vertex, unordered_set<int>& vertexCover,vector<int>& degrees);
        void debugVertexCover(unordered_set<int>& vertexCover);
        void getKAMIS(vector<int>& KAMIS);
        void debugMIS(vector<int>& vertexMIS);
        void readKAMIS(vector<int>& KAMIS,string const & fileName);
        void readWCNF(vector<int>& KAMIS,string const & fileName);
        void getMAXHS(vector<int>& KAMIS);
        void addCache(int l1);
        void adjustWeights();
        inline void increaseWeight(int v);
        void setE();
        void setTOP();
        void debugDegree(int vertex, vector<int> degrees,   unordered_set<int>& vertexCover);
        int numV;
        vector<double> weights;
        edgeList* adList;
        lookupTable* table; 
        // nodes choosen before in solution
        set<int> solution_prev;
};
#endif
