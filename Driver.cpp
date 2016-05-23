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

int main() {

    int readFlag = 0;
    int xtarget = 2;
    int ytarget = 25;
    int maxIter = 8000;
    int serializeFlag = 0;
    // some other initialization parameters
    int rodLen = 3;
    int phiNMax = 8;
    double gamma = 0.99;
    int numActuation = 3;
    double probThresh = 0.5;
    double defaultBigCost = 1000;
    // initial parameter for simulation
    int nstep = 12000;
    int outputFreq = 30;
    int x0 = 3;
    int y0 = 3;
    int simulateFlag = 0;
    int optimizeFlag = 0;
    int calFirstPassageTimeFlag = 1;
    int target_width;
    Solution::actionMode FirstPassageTimeOpt = Solution::slow;
    int FirstPassageTimeOpt_num;
    std::string maptag = "map3bit/map3bit";
    {
    std::string line;
    std::ifstream runfile;
    runfile.open("run.txt");
    getline(runfile, line);
    runfile >> readFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> xtarget;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> ytarget;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> maxIter;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> serializeFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> rodLen;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> phiNMax;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> gamma;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> numActuation;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> probThresh;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> defaultBigCost;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> nstep;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> outputFreq;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> x0;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> y0;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> simulateFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> optimizeFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> calFirstPassageTimeFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >>  FirstPassageTimeOpt_num;
    getline(runfile, line);    
     getline(runfile, line);
    getline(runfile, maptag);
    getline(runfile, line);    
     runfile >> target_width;
    
    FirstPassageTimeOpt = static_cast<Solution::actionMode>( FirstPassageTimeOpt_num);

    }
    
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

        if (optimizeFlag){
            sol.initialCost();
            sol.optimize(maptag + "opti_log", maxIter, xtarget, ytarget);
            sol.outputPolicy(maptag + "policydata.txt");
        }
        if (simulateFlag){
            sol.simulate(maptag + "policydata.txt", maptag + "probSolution", x0, y0, nstep, outputFreq);
        }
        if (calFirstPassageTimeFlag){
            sol.calFirstPassageTime(maptag + "policydata.txt", x0, y0, xtarget, ytarget, target_width, nstep, FirstPassageTimeOpt, outputFreq);
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
        if (optimizeFlag == 1) {
            sol2.optimize(maptag + "opti_log_dese", maxIter, xtarget, ytarget);
            sol2.outputConfigSet(maptag + "configdata_dese.txt");
            sol2.outputPolicy(maptag + "policydata_dese.txt");
        }

        if (simulateFlag == 1) {
            sol2.simulate(maptag + "policydata_dese.txt", maptag + "probSolution", x0, y0, nstep, outputFreq);
        }
        if (calFirstPassageTimeFlag){
        sol2.calFirstPassageTime(maptag + "policydata_dese.txt", x0, y0, xtarget, ytarget, target_width, nstep, FirstPassageTimeOpt, outputFreq);
        }
        
    }
    //	system("pause");
    std::cout << "finish!" << std::endl;
    return (0);
}



