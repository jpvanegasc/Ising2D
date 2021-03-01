/**
 * Monte Carlo Ising 2D Metropolis
 */
#include<iostream>
#include<cmath>

#include"random64.h"
#include"constants.h"

class SpinSystem{
    private:
        int s[L][L]; int E, M;
    public:
        void initialize_down(void);
        void metropolis_step(double beta, CRandom &ran);
        double get_E(void);
        double get_M(void);
};

inline double SpinSystem::get_E(void){return (double)E;}
inline double SpinSystem::get_M(void){return (double)std::abs(M);}

void SpinSystem::initialize_down(void){
    for(int i=0; i<L; i++)
        for(int j=0; j<L; j++)
            s[i][j] = -1;
    M = -L2; E = -2*L2;
}

void SpinSystem::metropolis_step(double beta, CRandom &ran){
    int n = (int)L2*ran.r(); int i = n%L, j = n/L; // Just one random number to save computing time
    int dE = 2*s[i][j]*(s[(i+1)%L][j] + s[(L+i-1)%L][j] + s[i][(j+1)%L] + s[i][(L+j-1)%L]);

    if(dE <= 0){
        s[i][j] *= -1; E += dE; M += 2*s[i][j];
        return;
    }
    else if(ran.r() < std::exp(-beta*dE)){
        s[i][j] *= -1; E += dE; M += 2*s[i][j];
        return;
    }
}

