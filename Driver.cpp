
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
    
    int readFlag = 1;
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
// initial parameter for simulation
    int nstep = 500;
    int outputFreq = 50;
    int x0 = 5;
    int y0 = 5;
        
//    std::string maptag="map2bitsimple_extend/map2bitsimple_extend";
        std::string maptag="map2bitsimple/map2bitsimple";
    if(!readFlag){
	Solution sol;
	sol.initialize(rodLen,phiNMax,gamma, numActuation, probThresh, defaultBigCost);
	sol.readBitmap(maptag+".txt");
	sol.readJumpMatrix("controljm");
	sol.constructObtacle();
	sol.constructConfigSet();
//	sol.outputConfigSet(maptag+"configdata.txt");

	sol.constructGraph();
//        sol.outputGraph(maptag+"graph.txt");
      {
		   std::ofstream ofs(maptag+"copyCoor.ser");
		   boost::archive::text_oarchive oa(ofs);
		   oa & sol;
     }        
        

        sol.initialCost();
	sol.optimize(maptag+"opti_log",maxIter,xtarget,ytarget);

	sol.outputPolicy(maptag+"policydata.txt");
        sol.outputNodeInfo(618);
        
        sol.simulate(maptag+"policydata.txt",maptag+"probSolution",x0,y0, nstep, outputFreq);

  
    }else{
        int optimizeFlag = 1;
        int simulateFlag = 1;
        Solution sol2;

	sol2.initialize( rodLen, phiNMax, gamma,  numActuation,  probThresh, defaultBigCost); 
        
        
 {
		   std::ifstream ifs(maptag+"copyCoor.ser");
		   boost::archive::text_iarchive ia(ifs);
		   ia & sol2;
 }

	
	std::cout << sol2.connectTo[0].size() << std::endl;
	std::cout << sol2.connectTo[1].size() << std::endl;
	if(optimizeFlag == 1){
	sol2.optimize(maptag+"opti_log_dese",maxIter,xtarget,ytarget);
	sol2.outputConfigSet(maptag+"configdata_dese.txt");
	sol2.outputPolicy(maptag+"policydata_dese.txt");
        }
	if(simulateFlag == 1){ 
        sol2.simulate(maptag+"policydata_dese.txt",maptag+"probSolution",x0,y0, nstep, outputFreq);
        }
    }
//	system("pause");
    std::cout << "finish!" << std::endl;
	return(0);
}



