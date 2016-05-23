#include<cmath>
#include<armadillo>
#include<vector>
#include<sstream>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include "omp.h"
int main(int argc, char** argv) {

//    double v_set0={0.125，0.25，0.5，1，1.5，2，3，4，6，10};
//    double t_set0={0.1，0.5，1，2，4};
//    double Dt_set0={0.2, 1, 5};
#ifdef OPENMP
    std::cout << "Number of processors available = " << omp_get_num_procs() << "\n";

    std::cout << "Number of max threads =              " << omp_get_max_threads() << "\n";
    
   
    int nthread;
    
    if (argc == 2){
        nthread = boost::lexical_cast<int>(argv[1]);
        nthread = std::min(nthread, omp_get_num_procs());
    } else {
        nthread = omp_get_num_procs();
    }
    omp_set_num_threads(nthread);
    std::cout << "threads used              " << nthread << "\n";
#endif
    
    
    std::vector<double> v_set {0.125,0.25,0.5,1,1.5,2,3,4,6,10};
    std::vector<double> t_set {0.5,1,2,4};
    std::vector<double> Dt_set {0.2,1,5,15};

    double v_ratio;
    double Dt_ratio;
    double t_ratio;

    int count = 0;
    
#ifdef OPENMP
#pragma omp parallel 
#pragma omp for
#endif    
    
    for (int v_idx = 0; v_idx < v_set.size(); v_idx++) {
        for (int t_idx = 0; t_idx < t_set.size(); t_idx++) {
            for (int Dt_idx = 0; Dt_idx < Dt_set.size(); Dt_idx++) {
#ifdef OPENMP 
#pragma omp critical
#endif             
                std::cout << count++ << std::endl;
                v_ratio = v_set[v_idx];
                t_ratio = t_set[t_idx];
                Dt_ratio = Dt_set[Dt_idx];
                double radius = 1000e-9; //radius of spherical particle, unit m
                double kb = 1.38e-23; // Boltzmann constant
                double T = 293.0; //temperature unit K
                double vis = 1e-3; // viscosity
                double Dt1 = 5.13e-13; // m^2/s
                double Dt2 = 4.02e-13;
                double Dr = 0.55;
                double control_t = 1; //s
                double vmax = 4.5e-6;
                Dt1 = Dt1*Dt_ratio;
                Dt2 = Dt2*Dt_ratio;
                Dr = Dr;
                control_t = control_t*t_ratio;
                vmax = vmax*v_ratio;
                int control_step;

                double v;
                double dist_thresh = vmax*control_t;


//                double sample_t = 100;

                long long nsample = (long) (10000 / control_t);
                double dt = control_t / 100;
                control_step = (int) (control_t / dt);
                double x = 0.0;
                double y = 0.0;
                double theta = 0.0;

                arma::vec xset, yset, thetaset;
                xset.zeros(nsample);
                yset.zeros(nsample);
                thetaset.zeros(nsample);

                for (int traj = 0; traj < nsample; traj++) {
                    xset(traj) = x;
                    yset(traj) = y;
                    thetaset(traj) = theta;
                    // v should equal to the projection
                    double proj = -(cos(theta) * x + sin(theta) * y);
                    if (proj < 0.5 * dist_thresh) {
                        v = 0.0;
                    } else {
                        v = vmax;
                    }
                    arma::mat randomdisp(5, control_step);
                    randomdisp.randn();
                    for (int i = 0; i < control_step; i++) {
                        x = x + sqrt(2 * Dt1 * dt) * randomdisp(0, i) * cos(theta) - sqrt(2 * Dt2 * dt) * randomdisp(1, i) * sin(theta) + v * cos(theta) * dt;
                        y = y + sqrt(2 * Dt1 * dt) * randomdisp(2, i) * sin(theta) + sqrt(2 * Dt2 * dt) * randomdisp(3, i) * cos(theta) + v * sin(theta) * dt;
                        theta = theta + sqrt(2 * Dr * dt) * randomdisp(4, i);
                    }
                }

                xset = xset / radius;
                yset = yset / radius;
#ifdef OPENMP 
#pragma omp critical
#endif             
                {
                std::stringstream ss;
                ss <<  v_ratio << "V" << t_ratio << "T" << Dt_ratio << "Dt";
                xset.save("xset" + ss.str() + ".txt", arma::raw_ascii);
                yset.save("yset" + ss.str() + ".txt", arma::raw_ascii);
            }
            }
        }
    }
    return 0;
}


