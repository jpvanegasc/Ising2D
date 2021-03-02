/**
 * Monte Carlo Ising 2D Metropolis
 */
#ifndef __ISING2D_H
#define __ISING2D_H

#include<iostream>
#include<fstream>
#include<cmath>

#include"random64.h"
#include"constants.h"

class Ising2D{
    private:
        CRandom *ran = NULL;

        int *system = NULL;

        int E, M;
    public:
        double beta = 0.0;

        Ising2D(unsigned long long seed);
        ~Ising2D();

        void initialize_down(void);
        void metropolis_step(void);

        void save(std::string filename, bool gnuplot);

        double get_E(void){return (double)E;}
        double get_M(void){return (double)std::abs(M);}
};

// 2D to 1D
#define SIZE (Lx*Ly)
#define X_MULT Ly

/**
 * Transform from 2D notation to 1D notation 
 * @return 1D macro-coordinate on array
 */
#define get_1D(ix, iy) (( (ix)*X_MULT) + (iy))


// 1D to 2D
/**
 * Transform from 1D notation to 2D notation
 * @return x coordinate
 */
#define get_ix(index) ((index)/X_MULT)
/**
 * Transform from 1D notation to 2D notation
 * @return y coordinate
 */
#define get_iy(index) ((index)%Ly)


#endif // __ISING2D_H