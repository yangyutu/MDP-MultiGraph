#include <iostream>
#include <fstream>
#include <string>
#include "utilities.h"
#include "Solution_est.h"
// include input and output archivers
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>


//using namespace std;

// we are going to serialize some components of the Solution to for the optimization
//	void testSerialization(const Solution &sol);

/* this code can only used for estimation of first passage time for diffusion
 cannot be used to find shortest path
 */
int main() {

    int readFlag = 0;
    // for map5bit 68,7; map3bit 2, 25; map6bit 3,3
    int xtarget = 75;
    int ytarget = 75;
    int maxIter = 8000;
    int serializeFlag = 0;
    // some other initialization parameters
    int rodLen = 3;
    int phiNMax = 1;
    double gamma = 0.99;
    int numActuation = 1;
    double probThresh = 0.5;
    double defaultBigCost = 1000;
    // initial parameter for simulation
    int nstep = 100000;
    int outputFreq = 30;
    // for map5bit 4,4, map3bit 3,3, map6bit 45,45
    int x0 = 10;
    int y0 = 60;
    int simulateFlag = 0;
//    int optimizeFlag = 1;
    int calFirstPassageTimeFlag = 1;

    Solution::actionMode FirstPassageTimeOpt = Solution::diffusion;
    std::string maptag = "map7bit/map7bit";
// if construct the graph from stratch
    if (!readFlag) {
        Solution sol;
        sol.initialize(rodLen, phiNMax, gamma, numActuation, probThresh, defaultBigCost);
        sol.readBitmap(maptag + ".txt");
        sol.readJumpMatrix("controljm");
        sol.constructObtacle();
        sol.constructConfigSet();
        //	sol.outputConfigSet(maptag+"configdata.txt");

        sol.constructGraph();
        //        sol.outputGraph(maptag+"graph.txt");
        
        if(serializeFlag){
            std::ofstream ofs(maptag + "copyCoor.ser");
            boost::archive::text_oarchive oa(ofs);
            oa & sol;
        }

//        if (optimizeFlag){
//            sol.initialCost();
//            sol.optimize(maptag + "opti_log", maxIter, xtarget, ytarget);
//            sol.outputPolicy(maptag + "policydata.txt");
//        }
        if (simulateFlag){
            sol.simulate(maptag + "policydata.txt", maptag + "probSolution", x0, y0, nstep, outputFreq);
        }
        if (calFirstPassageTimeFlag){
            sol.calFirstPassageTime(maptag + "policydata.txt", x0, y0, xtarget, ytarget, 2, nstep, FirstPassageTimeOpt);
        }

    } else {
//if directly read the map from the serialized data        


        Solution sol2;

        sol2.initialize(rodLen, phiNMax, gamma, numActuation, probThresh, defaultBigCost);
        {
            std::ifstream ifs(maptag + "copyCoor.ser");
            boost::archive::text_iarchive ia(ifs);
            ia & sol2;
        }


        std::cout << sol2.connectTo[0].size() << std::endl;
        std::cout << sol2.connectTo[1].size() << std::endl;
//        if (optimizeFlag == 1) {
//            sol2.optimize(maptag + "opti_log_dese", maxIter, xtarget, ytarget);
//           sol2.outputConfigSet(maptag + "configdata_dese.txt");
//            sol2.outputPolicy(maptag + "policydata_dese.txt");
//        }

        if (simulateFlag == 1) {
            sol2.simulate(maptag + "policydata_dese.txt", maptag + "probSolution", x0, y0, nstep, outputFreq);
        }
        if (calFirstPassageTimeFlag){
        sol2.calFirstPassageTime(maptag + "policydata_dese.txt", x0, y0, xtarget, ytarget, 2, nstep, FirstPassageTimeOpt);
        }
        
    }
    //	system("pause");
    std::cout << "finish!" << std::endl;
    return (0);
}



