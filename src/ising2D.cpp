#include"ising2D.h"

Ising2D::Ising2D(unsigned long long seed){
    ran = new CRandom(seed);

    system = new int[SIZE];
}

Ising2D::~Ising2D(){
    delete ran;
    delete[] system;
}

void Ising2D::initialize_down(void){
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++)
            system[get_1D(ix, iy)] = -1;
    M = -SIZE; E = -2*SIZE;
}

void Ising2D::metropolis_step(void){
    int i = (int)SIZE*ran->r();
    int ix = get_ix(i), iy = get_iy(i);

    int dE = 2*system[get_1D(ix, iy)]*(
        system[get_1D( (ix+1)%Lx, iy )] + system[get_1D( (Lx+ix-1)%Lx, iy )] +
        system[get_1D( ix, (iy+1)%Ly )] + system[get_1D( ix, (Ly+iy-1)%Ly )]
    );

    if(dE <= 0){
        system[get_1D(ix, iy)] *= -1;
        M += 2*system[get_1D(ix, iy)]; E += dE;
        return;
    }
    else if(ran->r() < std::exp(-beta*dE)){
        system[get_1D(ix, iy)] *= -1;
        M += 2*system[get_1D(ix, iy)]; E += dE;
        return;
    }
}

void Ising2D::save(std::string filename, bool gnuplot){
    std::ofstream file(filename);

    for(int ix=0; ix<Lx; ix++){
        for(int iy=0; iy<Ly; iy++){
            if(gnuplot) file << ix << ',' << iy << ',' << system[get_1D(ix, iy)] << '\n';
        }
        if(gnuplot) file << '\n';
    }

    file.close();
}

