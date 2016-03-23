#include<iostream>
#include<fstream>
#include<sstream>
#include<queue>
#include<cassert>
#include<cstdio>
#include<cstdlib>
//#include"pq.h"
#include"Solution_est.h"
#include"utilities.h"


using namespace std;

double Solution::PI = 3.1415926;

double Solution::round(double num) {
    return (num > 0.0) ? floor(num + 0.5) : ceil(num - 0.5);
}

void Solution::outputConfigSet(string filename) {

    ofstream os;
    os.open(filename);
    for (auto& g : configSet) {
        os << g->pos.x << "\t" << g->pos.y << "\t" << g->phi << "\t" << g->phi * 360 / phiNMax << endl;
    }
    os.close();
}

void Solution::outputGraph(string filename) {
    ofstream os;
    os.open(filename);
    for (auto &it : mapV) {
        GNode *g = it.second;
        for (int m = 0; m < numActuation; m++) {
            os << g->index << "\t";
            for (auto &e : connectTo[m][g->index]) {
                assert(g->index == e.from);
                os << e.to << "\t";
            }
        }
        os << endl;
    }
    os.close();
}

void Solution::outputNodeInfo(int index) {
    ofstream os;
    char buffer[50];
    sprintf(buffer, "%d", index);
    string filename(buffer);
    os.open(filename);
    GNode *g = mapV[index];

    os << g->index << "\t";
    os << g->pos.x << "\t";
    os << g->pos.y << "\t";
    os << g->phi << "\t";
    os << g->OptControl << "\t";
    os << g->JFunc << "\t";
    os << g->isolation << "\t";
    os << endl;
    double tempJFunc[3] = {0.0, 0.0, 0.0};
    for (int m = 0; m < numActuation; m++) {
        for (auto &e : connectTo[m][g->index]) {
            tempJFunc[m] += gamma * e.transProb * (mapV[e.to]->JFunc);
        }
    }
    os << tempJFunc[0] + costSet[0] << "\t";
    os << tempJFunc[1] + costSet[1] << "\t";
    os << tempJFunc[2] + costSet[2] << "\t";

    for (int m = 0; m < numActuation; m++) {
        for (auto &e : connectTo[m][g->index]) {
            GNode* g2 = (mapV[e.to]);
            os << m << "\t";
            os << g2->index << "\t";
            os << g2->pos.x << "\t";
            os << g2->pos.y << "\t";
            os << g2->phi << "\t";
            os << g2->OptControl << "\t";
            os << g2->JFunc << "\t";
            os << g2->isolation << "\t";
            os << e.transProb << "\t";
            os << connectTo[m][g->index].size() << "\t";
            os << endl;
        }
        tempJFunc[m] += costSet[m];
    }
    os.close();
}

void Solution::outputPolicy(string filename) {
    ofstream os;
    os.open(filename);

    for (auto &it : mapV) {
        GNode *g = it.second;
        if (g->isolation == 1) {
            if (g->OptControl != 0) cout << "wrong: isolation node has non-zero control option. the node is" << g->index << endl;
        }
        os << g->index << "\t";
        os << g->pos.x << "\t";
        os << g->pos.y << "\t";
        os << g->phi << "\t";
        os << g->OptControl << "\t";
        os << g->JFunc << "\t";
        os << g->isolation << "\t";
        os << endl;
    }
    os.close();
}

void Solution::readBitmap(string s) {

    ifstream is;
    is.open(s);
    string line;
    string data;
    int lineNum = 0;
    int row = 0;
    int temp;
    while (getline(is, line)) {
        lineNum++;
        stringstream linestream(line);

        if (lineNum == 1) {
            linestream>>rowBM;
            mapYMax = rowBM;
            cout << "The number of rows of the matrix are: " << rowBM << endl;
        } else if (lineNum == 2) {
            linestream>>colBM;
            mapXMax = colBM;
            cout << "The number of cols of the matrix are: " << colBM << endl;
            bitmap = new int*[rowBM];
            for (int i = 0; i < rowBM; i++) {
                bitmap[i] = new int[colBM];
            }
        } else {


            for (int i = 0; i < colBM; i++) {
                linestream>>temp;
                bitmap[row][i] = temp;
            }
            row++;
        }
    }
    ofstream os;
    os.open(s + "mapdata.txt");
    for (int i = 0; i < rowBM; i++) {
        for (int j = 0; j < colBM; j++) {
            if (bitmap[i][j] == 1) os << j << "\t" << i << endl;
            ;
        }
    }
    os.close();
    is.close();


    os.open(s + "walldata.txt");
    int count = 0;
    for (double i = -1; i <= rowBM; i += 0.5) {
        for (double j = -1; j <= colBM; j += 0.5) {
            if ((int) round(i) >= 0 && (int) round(i) < rowBM && (int) round(j) >= 0 && (int) round(j) < colBM) {
                if (bitmap[(int) round(i)][(int) round(j)] == 0) {
                    count++;
                    os << count << "\t" << j << "\t" << i << endl;
                }
            } else if (abs(i + 1) < 0.25 || abs(i - rowBM) < 0.25 ||
                    abs(j + 1) < 0.25 || abs(j - colBM) < 0.25) {
                count++;
                os << count << "\t" << j << "\t" << i << endl;
            }
        }

    }
    os.close();



}

void Solution::initialize(int rodLen0, int phiNMax0, double gamma0, int numActuation0, double probThresh0, double defaultBigCost0) {
    rodLen = rodLen0;
    phiNMax = phiNMax0;
    gamma = gamma0;
    numActuation = numActuation0;
    probThresh = probThresh0;
    defaultBigCost = defaultBigCost0;
    costSet = new double[numActuation];
    for (int i = 0; i < numActuation; i++) {
        costSet[i] = 1.0;
    }
    delx = 1;
    dely = 1;
    delphi = 1;

    connectTo = new vector< vector<Edge> >[numActuation];
    connectBy = new vector< vector<Edge> >[numActuation];
    jumpEventArr = new vector<JumpEvent>[numActuation];
    isolationSet = new unordered_set<int>[numActuation];
}

void Solution::readJumpMatrix(string tag) {
    ifstream file;
    string line;
    char fileend[10];
    string filename;
    for (int i = 0; i < numActuation; i++) {
        sprintf(fileend, "%d", i);
        filename = "diffuse_5s.txt";
        file.open(filename);
        while (getline(file, line)) {
            JumpEvent je;
            stringstream linestream(line);
            linestream >> je.x;
            linestream >> je.y;
            linestream >> je.phi;
            linestream >> je.prob;
            jumpEventArr[i].push_back(je);
        }
        file.close();
    }
}

void Solution::constructObtacle() {
    for (int i = 0; i < rowBM; i++)
        for (int j = 0; j < colBM; j++) {
            if (bitmap[i][j] == 0) {
                obstacleSet.insert(CoorPair(j, i));
            }
        }
    numObs = obstacleSet.size();
}

void Solution::constructConfigSet() {
    int count = 0;
    for (int i = 0; i < rowBM; i++)
        for (int j = 0; j < colBM; j++) {
            // if there is not obstacle
            if (bitmap[i][j] == 1) {
                for (int k = 0; k < phiNMax; k++) {
                    int phi = k;
                    if (!isOverlapObstacle(j, i, phi)) {
                        configSet.push_back(new GNode(j, i, phi));
                        configIndex[Config(j, i, phi)] = count;
                        count++;
                    }
                }
            }
        }
    numV = configSet.size();
}

bool Solution::outOfBoundary(int x, int y) {
    if (x < 0) return true;
    if (y < 0) return true;
    if (x >= mapXMax) return true;
    if (y >= mapYMax) return true;
    return false;
}

// here since we are estimating point mass diffusion, we are not considering the geometry effect
bool Solution::isOverlapObstacle(int x0, int y0, int phi) {
    int xtemp;
    int ytemp;
    
    xtemp = x0;
    ytemp = y0;

    for (int i = -rodLen / 2; i <= rodLen / 2; i++) {
        for (int j = -rodLen / 2; j <= rodLen / 2; j++) {
            xtemp = x0 + i;
            ytemp = y0 + j;
        if (outOfBoundary(xtemp, ytemp)) return true;
        if (obstacleSet.find(CoorPair(xtemp, ytemp)) != obstacleSet.end()) return true;
        }
    }

    if (outOfBoundary(xtemp, ytemp)) return true;
    if (obstacleSet.find(CoorPair(xtemp, ytemp)) != obstacleSet.end()) return true;

    return false;
}

void Solution::constructGraph() {
    int count = 0;
    for (auto n : configSet) {
        mapV[count] = n;
        n->index = count;
        count++;
    }
    numV = count;

    for (int i = 0; i < numActuation; i++) {
        connectTo[i].reserve(numV);
        connectBy[i].reserve(numV);
        for (int j = 0; j < numV; j++) {
            vector<Edge> v1;
            vector<Edge> v2;
            connectTo[i].push_back(v1);
            connectBy[i].push_back(v2);
        }
    }
    for (int i = 0; i < numV; i++) {
        cout << "working on node: " << i << endl;
        for (int m = 0; m < numActuation; m++) {
            //adding edges due to actuation type m
            addEdges(i, m);
        }
    }
}

void Solution::addEdges(int i, int m) {
    double sum_prob = 0;
    GNode &g = *(mapV[i]);
    int newx;
    int newy;
    int newphi;
 /*   
    for (auto &je : jumpEventArr[m]) {
        newx = g.pos.x + (int) round(je.x * cos(g.phi * PI * 2.0 / (double) phiNMax) -
                je.y * sin(g.phi * PI * 2.0 / (double) phiNMax));
        newy = g.pos.y + (int) round(je.x * sin(g.phi * PI * 2.0 / (double) phiNMax) +
                je.y * cos(g.phi * PI * 2.0 / (double) phiNMax));
        newphi = (g.phi + je.phi + phiNMax) % phiNMax;
        if (!isPathIntersectObstacle(g.pos.x, g.pos.y, g.phi, newx, newy, newphi)) {
            auto index = configIndex.find(Config(newx, newy, newphi));
            if (index != configIndex.end()) {
                sum_prob += je.prob;
            }
        }
    }
    // if the probability from one node to all other nodes1 is larger than the thresh and it is not diffusion mode
    if (sum_prob > this->probThresh) {
        for (auto &je : jumpEventArr[m]) {
            newx = g.pos.x + (int) round(je.x * cos(g.phi * PI * 2.0 / (double) phiNMax) -
                    je.y * sin(g.phi * PI * 2.0 / (double) phiNMax));
            newy = g.pos.y + (int) round(je.x * sin(g.phi * PI * 2.0 / (double) phiNMax) +
                    je.y * cos(g.phi * PI * 2.0 / (double) phiNMax));
            newphi = (g.phi + je.phi + phiNMax) % phiNMax;
            if (!isPathIntersectObstacle(g.pos.x, g.pos.y, g.phi, newx, newy, newphi)) {
                auto index = configIndex.find(Config(newx, newy, newphi));
                if (index != configIndex.end()) {
                    int j = index->second;
                    Edge e(i, j, costSet[m], m, je.prob / sum_prob);
                    connectBy[m][j].push_back(e);
                    connectTo[m][i].push_back(e);
                    numE++;
                }
            }
        }
// if the control option is diffusion then no matter how small the effective of 
//        transition probablity is, we allow the existence of such control
    } else if (m == 0) {
        for (auto &je : jumpEventArr[m]) {
            newx = g.pos.x + (int) round(je.x * cos(g.phi * PI * 2.0 / (double) phiNMax) -
                    je.y * sin(g.phi * PI * 2.0 / (double) phiNMax));
            newy = g.pos.y + (int) round(je.x * sin(g.phi * PI * 2.0 / (double) phiNMax) +
                    je.y * cos(g.phi * PI * 2.0 / (double) phiNMax));
            newphi = (g.phi + je.phi + phiNMax) % phiNMax;
            if (!isPathIntersectObstacle(g.pos.x, g.pos.y, g.phi, newx, newy, newphi)) {
                auto index = configIndex.find(Config(newx, newy, newphi));
                if (index != configIndex.end()) {
                    int j = index->second;
                    Edge e(i, j, costSet[m], m, je.prob / sum_prob);
                    connectBy[m][j].push_back(e);
                    connectTo[m][i].push_back(e);
                    numE++;
                }
            }
        }
    }
   */ 
    
    if (m == 0) {
        for (auto &je : jumpEventArr[m]) {
            newx = g.pos.x + je.x;
            newy = g.pos.y + je.y;
            newphi = 0;
            if (!isPathIntersectObstacle(g.pos.x, g.pos.y, g.phi, newx, newy, newphi)) {
                auto index = configIndex.find(Config(newx, newy, newphi));
                if (index != configIndex.end()) {
                    sum_prob += je.prob;
                }
            }
        }
        
        for (auto &je : jumpEventArr[m]) {
            
//            newx = g.pos.x + (int) round(je.x * cos(g.phi * PI * 2.0 / (double) phiNMax) -
//                    je.y * sin(g.phi * PI * 2.0 / (double) phiNMax));
//            newy = g.pos.y + (int) round(je.x * sin(g.phi * PI * 2.0 / (double) phiNMax) +
//                    je.y * cos(g.phi * PI * 2.0 / (double) phiNMax));
//            newphi = (g.phi + je.phi + phiNMax) % phiNMax;
            newx = g.pos.x + je.x;
            newy = g.pos.y + je.y;
            newphi = 0;
            
            if (!isPathIntersectObstacle(g.pos.x, g.pos.y, g.phi, newx, newy, newphi)) {
                auto index = configIndex.find(Config(newx, newy, newphi));
                if (index != configIndex.end()) {
                    int j = index->second;
                    Edge e(i, j, costSet[m], m, je.prob / sum_prob);
                    connectBy[m][j].push_back(e);
                    connectTo[m][i].push_back(e);
                    numE++;
                }
            }
        }
    }
}

bool Solution::isPathIntersectObstacle(int x, int y, int phi, int newx, int newy, int newphi) {

    int xtemp;
    int ytemp;
    if (x == newx && y == newy) {
        if (isOverlapObstacle(x, y, newphi)) return true;
    } else if (x == newx) {
        for (int i = 0; i <= abs(newy - y); i++) {
            xtemp = x;
            ytemp = (int) round(y + (newy - y) / abs(newy - y) * i);
            if (isOverlapObstacle(xtemp, ytemp, phi)) return true;
            if (isOverlapObstacle(xtemp, ytemp, newphi)) return true;
            //				if(obstacleSet.find(CoorPair(xtemp,ytemp))!=obstacleSet.end()) return true;
        }
    } else if (y == newy) {
        for (int i = 0; i <= abs(newx - x); i++) {
            xtemp = (int) round(x + (newx - x) / abs(newx - x) * i);
            ytemp = y;
            if (isOverlapObstacle(xtemp, ytemp, phi)) return true;
            if (isOverlapObstacle(xtemp, ytemp, newphi)) return true;
            //				if(obstacleSet.find(CoorPair(xtemp,ytemp))!=obstacleSet.end()) return true;
        }
    } else {
        //		double slope=((double)newy-(double)y)/((double)newx-(double)x);
        double len = sqrt((x - newx)*(x - newx)+(y - newy)*(y - newy));
        for (int i = 0; i <= len; i++) {
            xtemp = (int) round(x + i * ((double) newx - (double) x) / len);
            ytemp = (int) round(y + i * ((double) newy - (double) y) / len);

            if (isOverlapObstacle(xtemp, ytemp, phi)) return true;
            if (isOverlapObstacle(xtemp, ytemp, newphi)) return true;
            //			if(obstacleSet.find(CoorPair(xtemp,ytemp)) != obstacleSet.end()) return true;
        }
    }

    return false;
}

void Solution::initialCost() {
    //	double costToGo=0.0;
    for (auto &it : mapV) {
        GNode *g = it.second;
        g->OptControl_old = 0;
        g->OptControl = 0;
// at the beginning, the cost of every node is set to a big number
        g->JFunc = defaultBigCost;
        g->JFunc_old = g->JFunc;
        g->isolation = 0;
    }
}

void Solution::optimize(string logfile, int maxIter0, int xtarget0, int ytarget0) {
    ofstream ofs;
    ofs.open(logfile);

    int maxIter = maxIter0;
    int xtarget = xtarget0;
    int ytarget = ytarget0;

    queue<int> Q;
    int GNodeIndex;
    double totalValue;
    int *marked = new int[numV];
    double *tempJFunc = new double[numActuation];

    initialCost();
    refineJump();
    for (int iter = 0; iter < maxIter; iter++) {
        cout << iter << endl;
        for (int i = 0; i < numV; i++) marked[i] = 0;
        //		Q.		

        for (int phi = 0; phi < phiNMax; phi++) {
// first find the target node
            auto index = configIndex.find(Config(xtarget, ytarget, phi));
            if (index != configIndex.end()) {
                Q.push(index->second);
                marked[index->second] = 1;
                GNode &g = *(mapV[index->second]);

// if we want fast convergence and sacrifice the theoretical soundness, uncomment the following                
                //			g.JFunc = 0.0;
                //			g.OptControl=0;
                //			g.JFunc_old = 0.0;
                //			g.OptControl_old=0;
                //                        g.isolation = 0;
            }
        }

        while (!Q.empty()) {
            GNodeIndex = Q.front();
            Q.pop();		
            GNode& g = *(mapV[GNodeIndex]);
            if (g.isolation == 1) {
                cout << "wrong isolated node in queue: " << g.index << endl;
            }

            for (int m = 0; m < numActuation; m++) {
                tempJFunc[m] = 0.0;
                for (auto &e : connectTo[m][GNodeIndex]) {
                    tempJFunc[m] += gamma * e.transProb * (mapV[e.to]->JFunc);
                }
                //		if Gnode is not connect to any other nodes via actutaion m
                if (connectTo[m][GNodeIndex].size() == 0.0) {
                    tempJFunc[m] = defaultBigCost;
                }
                tempJFunc[m] += costSet[m];
                if (g.pos.x == xtarget && g.pos.y == ytarget) tempJFunc[m] -= costSet[m];
            }

            int minIndex = 0;
            double minVal = g.JFunc_old;
            for (int m = 0; m < numActuation; m++) {
                if (tempJFunc[m] < minVal) {
                    minIndex = m;
                    minVal = tempJFunc[m];
                }

            }
// take newer better value and decision            
            g.OptControl_old = g.OptControl;
            g.OptControl = minIndex;
            g.JFunc_old = g.JFunc;
            g.JFunc = minVal;

            for (int m = 0; m < numActuation; m++) {
                for (auto &e : connectBy[m][GNodeIndex]) {
                    assert(e.to == GNodeIndex);
                    if (!marked[e.from]) {
                        Q.push(e.from);
                        marked[e.from] = 1;
                    }

                }
            }
        }

        if (converged(gamma, totalValue)) break;
        ofs << iter << "\t" << totalValue << endl;

    }

    ofs.close();


    double diffVal = 0.0;
    int opt;
    for (auto &it : mapV) {
        GNode *g = it.second;
        if (!(g->isolation)) {
            opt = g->OptControl;
            diffVal = g->JFunc;
            for (auto &e : connectTo[opt][g->index])
                diffVal -= gamma * e.transProb * (mapV[e.to]->JFunc);
            diffVal -= costSet[opt];
        }
        if (diffVal > 1e-6) {
            cout << g->index << "\t" << diffVal << endl;
        }
    }
}

bool Solution::converged(double gamma, double &totalValue) {
    double epi = 1e-7;
    int diff = 0;
    double totalVal = 0.0;
    double diffVal = 0.0;
    for (auto &it : mapV) {
        GNode *g = it.second;
        diff += abs(g->OptControl - g->OptControl_old);
        diffVal += abs(g->JFunc - g->JFunc_old);
        if (!(g->isolation)) totalVal += g->JFunc;
    }
    totalValue = totalVal / numV;
    diffVal /= numV;
    // here the value iteration condition is from approximate dynamic programming	
    return diff == 0 && (diffVal < (epi * (1 - gamma) / gamma / 2.0));
}

/** Here I refine the graph, specifically the jump events
 *  1) if a node cannot jump to another nodes with probability greater than the 
 * thresh at any actuation, then the node is considered isolated. that means 
 * connectTo[m][idx] has size 0 for all m;
 *  2) if a node is isolated, then I set it to be reflective. That is other
 *  nodes should not be able to jump it. Then I will break the edge.  
    3) how to add penalty if they go to the isolated node?
 *  4) go to isolated node is different from hitting the obstacle
 **/
void Solution::refineJump() {
    int isolateFlag;
    // if the node is not connecting to anything, then it is isolated    
    for (int i = 0; i < numV; i++) {
        isolateFlag = 1;
        for (int m = 0; m < numActuation; m++) {
            if (connectTo[m][i].size() != 0) isolateFlag = 0;
        }
        if (isolateFlag == 1) mapV[i]->isolation = 1;
    }

    // here I set the transition probability to the isolated node to be 0.0 and
    //   and renormalize the transition probability
    /*   double normFactor;
       for(int m=0;m<numActuation;m++){
     
         for(int i=0;i < numV; i++){  
                normFactor = 0.0;
               for(auto &e: connectTo[m][i]){
                   if(mapV[e.to]->isolation == 1) e.transProb = 0.0;            
                   normFactor += e.transProb;
               }
            
               for(auto &e: connectTo[m][i]){           
                   e.transProb /= normFactor;
               }
        
           } 
       }  
     */
}

void Solution::simulate(string policyname, string outputfile, int x0, int y0, int nstep, int outputfreq) {
    ifstream is;
    string line;


/*    is.open(policyname);
    int dum;
    int index;
    for (int i = 0; i < numV; i++) {
        getline(is, line);
        stringstream linestream(line);
        linestream >> index;
        linestream >> dum;
        linestream >> dum;
        linestream >> dum;
        linestream >> (mapV[index]->OptControl);
        linestream >> (mapV[index]->JFunc);
        linestream >> (mapV[index]->isolation);
    }

*/

    int step = 0;

    map<int, double> oldSol, newSol;

    for (int i = 0; i < numV; i++) {
        oldSol[i] = 0.0;
        newSol[i] = 0.0;

    }

    auto it = configIndex.find(Config(x0, y0, 0));
    if (it == configIndex.end()) {
        cout << "initial configuration not admissible" << endl;
        exit(100);
    }

    int initialIndex = it->second;
    oldSol[initialIndex] = 1.0;
    newSol[initialIndex] = 1.0;
    int count = 0;
    int countedge;
    double sumtemp;
    while (step < nstep) {
        //  update the probability of every node
        //  probability will flow into 
        // and flow out        
        sumtemp = 0.0;
        countedge = 0;
        for (int i = 0; i < numV; i++) {
//            int option = mapV[i]->OptControl;
//      0 is the diffusion mode
            int option = 0;
            for (Edge &e : connectTo[option][i]) {
                newSol[i] -= oldSol[i] * e.transProb;
                newSol[e.to] += oldSol[i] * e.transProb;

            }

        }
        cout << "simulate step: " << step << endl;
        if (step == 0 || (step + 1) % outputfreq == 0) {
            count++;
            outputProbDist(outputfile, count, newSol);

        }
        for (int i = 0; i < numV; i++) {
            oldSol[i] = newSol[i];
        }
        step++;
    }
}

void Solution::outputProbDist(string outputfile, int count, const map<int, double> &newSol) {
    std::stringstream fileend;
    string filename;
    ofstream os;
    fileend << count;
    filename = outputfile + fileend.str() + ".txt";
    os.open(filename);

    for (auto &it : mapV) {

        GNode *g = it.second;
	
        os << g->index << "\t";
        os << g->pos.x << "\t";
        os << g->pos.y << "\t";
        os << g->phi << "\t";
        os << newSol.at(g->index) << "\t";
        os << g->isolation << "\t";
        os << endl;
    }
    os.close();
}

void Solution::calFirstPassageTime(string policyname, int x0, int y0,  int xtarget, int ytarget, int width, int nstep, actionMode act) {
    ifstream is;
    string line;
    is.open(policyname);
    int dum;
    int index;
    for (int i = 0; i < numV; i++) {
        getline(is, line);
        stringstream linestream(line);
        linestream >> index;
        linestream >> dum;
        linestream >> dum;
        linestream >> dum;
        linestream >> (mapV[index]->OptControl);
        linestream >> (mapV[index]->JFunc);
        linestream >> (mapV[index]->isolation);
    }



    int step = 0;

    map<int, double> oldSol, newSol;

    for (int i = 0; i < numV; i++) {
        oldSol[i] = 0.0;
        newSol[i] = 0.0;

    }

    auto it = configIndex.find(Config(x0, y0, 0));
    if (it == configIndex.end()) {
        cout << "initial configuration not admissible" << endl;
        exit(100);
    }

    int initialIndex = it->second;
    oldSol[initialIndex] = 1.0;
    newSol[initialIndex] = 1.0;
    int count = 0;
    int countedge;
    double prob_mass_sum ;
    int option;

    ofstream os;
    os.open("first_pass_time.txt");
    while (step < nstep) {
        //  update the probability of every node
        //  probability will flow into 
        // and flow out        
        prob_mass_sum = 0.0;
        countedge = 0;
        for (int i = 0; i < numV; i++) {
            
            if(oldSol[i] < 1e-12){
                continue;
            }
            if(act == diffusion) {
                option = 0;
            }
            // for nodes in the adsorbing region or nodes have too low probabiliy mass, we last neglect it.
            if (!inAdsorbingRegion(mapV[i], xtarget, ytarget, width)) {
                for (Edge &e : connectTo[option][i]) {
                    newSol[i] -= oldSol[i] * e.transProb;
                    newSol[e.to] += oldSol[i] * e.transProb;

                }
            }

        }
        
        cout << "simulate step in first passage time calculation: " << step << endl;

        // we need to sum up all the survival probability mass
        for (int i = 0; i < numV; i++) {
            // for nodes in the adsorbing region or nodes have too low probabiliy mass, we last neglect it.
            if (!inAdsorbingRegion(mapV[i], xtarget, ytarget, width) ){
                prob_mass_sum += newSol[i];
            }
        }
        
        os << step << "\t";
        os << prob_mass_sum << endl;

        for (int i = 0; i < numV; i++) {
            oldSol[i] = newSol[i];
        }
        step++;
        
    }
}



bool Solution::inAdsorbingRegion(GNode *g, int xtarget, int ytarget, int width){

	if ( abs(g->pos.x - xtarget) < width && abs(g->pos.y - ytarget) < width) {
		return true;
	}

	return false;
}