#ifndef GPL_DATASTRUCTURE_H
#define GPL_DATASTRUCTURE_H
#include<iostream>
#include <vector>
#include <unordered_map>
#include<assert.h>
using namespace std;
class lookupTable{
    public:
        lookupTable(int numV);
        int insert(int lID);
        int lookUpLID(int vID);
        int lookUpVID(int lID);
        void debug();
        void print();
    private:
        int numV;
        vector<int> labelID;
        unordered_map<int,int> vertexID;
};
#endif 