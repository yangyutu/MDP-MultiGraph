#ifndef SOLUTION_H
#define SOLUTION_H


#include<vector>
#include<iostream>
#include<map>
#include<unordered_set>
#include<unordered_map>
#include<string>
#include "utilities.h"
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
class Solution{
	
public:
    
        enum actionMode{ optimal=0, diffusion=1, slow=2, fast=3};

	std::vector<GNode *> configSet;
	std::vector<int> finishSet;
	double *costSet;
	double *distTo;
	Edge *edgeTo;
	std::unordered_set<CoorPair,CoorPairHash,CoorPairEqual> obstacleSet;
	static double PI;
	std::unordered_map <Config, int, ConfigHash, ConfigEqual> configIndex;
public:
	double delx;
	double dely;
	double delphi;
	int **bitmap;
	int rowBM;
	int colBM;
	int rodLen;
/** number of angle state */
	int phiNMax;
	int numV;
	int numE;
/** number of obstacles
    */
	int numObs;
/** number of actuation type*/
	int numActuation;

	double probThresh;
	double gamma;
	int maxIter;
	int xtarget;
	int ytarget;
	int mapXMax;
	int mapYMax;
        
        double defaultBigCost;
	std::map<int,double> controlMap;
	std::unordered_map<int,GNode *> mapV;
        std::unordered_set<int> *isolationSet;
	std::vector< std::vector<Edge> > *connectTo;
	std::vector< std::vector<Edge> > *connectBy;
	std::vector< JumpEvent > *jumpEventArr;
	
	
	   friend class boost::serialization::access;
 
 // only Serialize selected members for optimization purpose
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & mapV;
      ar & numV;
//      ar & configSet;
      ar & gamma; 
       ar & numActuation;
      for(int i=0; i< numActuation; i++){
      ar & connectTo[i];
      ar & connectBy[i];
      }
     
      ar & phiNMax;
      ar & configIndex;   
    }

public:
 Solution(){};
	void readJumpMatrix(std::string tag);
	void outputConfigSet(std::string filename);
	void outputPolicy(std::string filename);
	void outputGraph(std::string);
	void readBitmap(std::string s);
	void constructObtacle();
	void initialize(int rodLen,int phiNMax,double gamma, int numActuation, double probThresh, double defaultBigCost);
	void constructConfigSet();
	bool isOverlapObstacle(int x0,int y0,int phi);
	void constructGraph();
	void addEdges(int i,int m);
	bool isControlNeighbor(GNode & g1,GNode& g2, double move);
	bool isPathIntersectObstacle(int x, int y, int phi, int newx, int newy, int newphi);
	bool containsEdge(int i,int j,int m);
	int  mapPos(int x0,int y0,int phi0);
	void calPath(int s);
	void relax(int v);
	void printPath(int w,int s);
	void pirntPathHelper(int w,int s, std::ostream &os);
	double round(double);
	void optimize(std::string,int maxIter0, int xtarget0, int ytarget0);
	bool converged(double gamma, double &totalValue);
	void initialCost();
	bool outOfBoundary(int x, int y);
        void refineJump();
        void outputNodeInfo(int index);

	void simulate(std::string policyname,std::string, int x0, int y0, int nstep, int outputfreq);
	void outputProbDist(std::string outputfile,int count, const std::vector<double> &newSol);
        void calFirstPassageTime(std::string policyname, int x0, int y0, int xtarget, int ytarget, int width, int nstep, actionMode act,int);
        bool inAdsorbingRegion(GNode *g, int xtarget, int ytarget, int width);
        
};

#endif
