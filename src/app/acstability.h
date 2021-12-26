#ifndef ACSTABILITY_H
#define ACSTABILITY_H
#include <stdio.h>
#include <iostream>
#include "fifo_map.hpp"
#include "json.hpp"
#include <unordered_map>
template<class K, class V, class dummy_compare, class A> using json_fifo_map = nlohmann::fifo_map<K, V, nlohmann::fifo_map_compare<K>, A>; 
using json = nlohmann::basic_json<json_fifo_map>;
using namespace std;
namespace test{
    struct Performance{
        string name;
        long time;
        int solutionSize;
        int remainingLabels = 0;
        double solutionWeight;
        //+++++++++++++++++++++++++debug code++++++++++++++++++++++++++++++++++++
        void print(){
            cout<<"name: "<< name << endl; 
            cout<< "time: "<< time << endl; 
            cout<< "solutionSize: "<< solutionSize<< endl; 
            cout<< "remainingLables: "<< remainingLabels << endl; 
            cout<< "solutionWeight: "<< solutionWeight << endl; 
        }
        json convertJSON(){
            json node;
            node["name"] = name;
            node["time"] = time; 
            node["solutionSize"] = solutionSize;
            node["remainingLabels"] = remainingLabels;
            node["solutionWeight"] =  solutionWeight;
            return node;
        }
        void appendJSON(json node, string name){
            node[name] = convertJSON();
        }
        //-----------------------debug code--------------------------------------
    };
}
#endif