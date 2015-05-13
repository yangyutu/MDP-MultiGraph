
#include <iostream>
#include <fstream>
#include <string>
#include "utilities.h"
#include "Solution.h"
// include input and output archivers
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>


//using namespace std;

	// we are going to serialize some components of the Solution to for the optimization
//	void testSerialization(const Solution &sol);

int main(){
    
    int readFlag = 0;
    int xtarget = 5;
    int ytarget = 5;
    int maxIter = 2000;
// some other initialization parameters
        int rodLen=3;
	int phiNMax=8;
	double gamma=0.99;
	int numActuation = 3;
	double probThresh = 0.5;
        double defaultBigCost = 1000;
        
//    std::string maptag="map2bitsimple_extend/map2bitsimple_extend";
        std::string maptag="map2bitsimple/map2bitsimple";
    if(!readFlag){
	Solution sol;
	sol.initialize(rodLen,phiNMax,gamma, numActuation, probThresh, defaultBigCost);
	sol.readBitmap(maptag+".txt");
	sol.readJumpMatrix("controljm");
	sol.constructObtacle();
	sol.constructConfigSet();
	sol.outputConfigSet(maptag+"configdata.txt");

	sol.constructGraph();
        sol.outputGraph(maptag+"graph.txt");

        sol.initialCost();
	sol.optimize(maptag+"opti_log",maxIter,xtarget,ytarget);

	sol.outputSolution(maptag+"policydata.txt");
        sol.outputNodeInfo(86);
        
    

     {
		   std::ofstream ofs(maptag+"copyCoor.ser");
		   boost::archive::text_oarchive oa(ofs);
		   oa & sol;
     }   
    }else{
        Solution sol2;

	sol2.initialize( rodLen, phiNMax, gamma,  numActuation,  probThresh, defaultBigCost); 
        
        
 {
		   std::ifstream ifs(maptag+"copyCoor.ser");
		   boost::archive::text_iarchive ia(ifs);
		   ia & sol2;
 }

	
	std::cout << sol2.connectTo[0].size() << std::endl;
	std::cout << sol2.connectTo[1].size() << std::endl;
	
	sol2.optimize(maptag+"opti_log_dese",xtarget,ytarget,maxIter);
	sol2.outputConfigSet(maptag+"configdata_dese.txt");
	sol2.outputSolution(maptag+"policydata_dese.txt");
	
    }
//	system("pause");
    std::cout << "finish!" << std::endl;
	return(0);
}



