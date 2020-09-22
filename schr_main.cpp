#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <complex.h>
#include <chrono> //time meansurement

#include <math.h> //-lm
#include <armadillo> //-larmadillo
#include <omp.h> //-lopenm, for loops and armadillo
#include <fftw3.h> //-lfftw3 armadillo uses 
#include <cmath>
#include <lapacke.h> //-llapacke for armadillo
#include <cblas.h> //-lblas for ...

#define ARMA_USE_OPENMP

#define BASE 2
#define hbar 1
#define mass 1

using namespace arma;
using namespace std;

/*
You need to create a "out" folder for files and "figs" for python plots;

g++ schr_main.cpp -larmadillo -fopenmp -lm -lfftw3 -llapacke -Ofast -lblas -o BINARY  --> to generate the binary with name "BINARY"
./BINARY --> to run with the current configuration


Hi, this code use the split operator method for the potential "Potential" where potential isn't time dependent, it use mainly the armadillo library with lapacke, fttw3 and blas in the backend. Some consideratios:
-> x and y are a square lattice;
-> to improve speed the code use paralelism for all possible "for" statements (I'm thinking you'll use a huge system, for small system N<2^12 serial for is faster):
    --> for N = 2^12 with openMP the code take 10.7s without it takes 15.3s (for each time step);
    --> for N = 2^13 with openMP the code take 56.31s without it takes 46.85s (for each time step);
    ///-------------MY PC HAS ONLY 4 CORES--------------///
-> to avoid paralel for's just comment all ''#pragma...'' lines
-> dont change the BASE value, it's really necessary for the code speed
-> the plots.py program generate a video with seaborn and ffmpeg;
-> some lines can be hard to understand because I avoid to use "if" statement;

FOR ANY RECOMENDATION, COMENTARY OR RECOMENDATION  => joaocassianox7x@gmail.com

*/


void Initial_Condition(cx_cube& IC,int N, int i, int j, double a, double delta)
{

    double width = 10.; //width of the wave packet 2-d
    double center_x = 15.; //center of gaussian package in x axis
    double center_y = 300.; //center of gaussian package in y axis

    double momentum_x = +35.0; //momentum of plane wave in x axis, in quantum mechanics kx
    double momentum_y = +0.0; //momentum of plane wave in y axis, in quantum mechanics ky

    double aux_x, aux_y; //auxiliar variables to define x_i and y_u

    #pragma omp parallel for shared(IC, width, center_x,center_y,momentum_x, momentum_y) private(i,j)
    for(i = 0; i<N; i++)
    {
        #pragma omp parallel for shared(IC, width, center_x,center_y,momentum_x, momentum_y,i) private(j)
        for(j = 0; j<N; j++)
        {
            aux_x = a+i*delta - center_x; 
            aux_y = a+j*delta - center_y;
            IC(i,j,0) = 10*(width*width/M_PIl)*cexp(-(0.5/width/width)*(aux_x*aux_x+aux_y*aux_y)+1I*(-aux_x*momentum_x/hbar-aux_y*momentum_y/hbar)); //initial condition, Gaussian package
        }
    } //end of parallel for :)//

}

void Potential(cx_mat& POT,int N, int i, int j, double delta)
{ 
    #pragma omp parallel for shared(POT) private(i,j) //parallel for to fill the potential matrix faster 
    for(i = 0; i<N; i++)
    {
        #pragma omp parallel for shared(POT,i) private(j)
        for(j = 0; j<N; j++)
        {
            POT(i,j) = 1000*(delta*i>=450 && delta*i<=455) + 1000*(delta*i>=250 && delta*i<=280) + 1000*(delta*i>=150 && delta*i<=170); //potential, three diferent almost infinity walls
        }
    } //end of parallel for :)//

}

void Laplace(cx_mat& DER, int N, int i, int j, double delta)
{
    double k = 2*M_PIl/delta;
    double deltak = k/(double)N;

    double k_x, k_y;
    #pragma omp parallel for shared(DER,deltak) private(i,j) //parallel for to fil the laplacian faster
    for(i = 0; i<N; i++)
    {
        k_x = i*deltak - 1*(i>=N/2)*N*deltak;
        #pragma omp parallel for shared(DER,deltak,i) private(j) //parallel for to fil the laplacian faster
        for(j = 0; j<N; j++)
        {
            k_y = j*deltak - 1*(j>=N/2)*N*deltak;
            DER(i,j) = k_x*k_x+k_y*k_y; //matrix of kx and ky in reciprocal space
        }
    } //end of parallel for :)//
}


int main()
{
    auto start = std::chrono::high_resolution_clock::now();

    //grid values
    double a = 0.; // "x" and "y" are in [a,b] interval
    double b = +500.; // ...
    int N = pow(BASE,14); //number of points, BASE = 2, so the points are a power of 2 for the FFTW3 work faster with divide and conquest
    double delta = (b-a)/(double)N; //deltax = deltay = delta


    //time values
    double t0 = 0.; //initial time
    double tf = 10.; //final time
    double dt = 5; // delta t for our problem
    int N_Ts = 1;//(int)(tf-t0)/dt; // number of time steps

    int i,j;

    //matrix and tensor necessary for psi_{ij,t} and v_{ij}
    cx_cube IC(N,N,2,fill::zeros); //vector for initial condition //zeros for memory allocation
    cx_mat POT(N,N,fill::ones); //vector for potential //zeros for memory allocation, more in http://arma.sourceforge.net/docs.html#Cube
    cx_mat DER(N,N,fill::zeros); //laplace operator in reciprocal space

    Laplace(DER,N, i, j,delta);

    
    Potential(POT,N, i, j, delta);
    Initial_Condition(IC,N, i, j,a,delta);


    FILE *fc = fopen("potencial.dat","w"); //write the potential for python program be able to polot it
    for(int i = 0; i<N; i++)
    {
            for(int j = 0; j<N; j++)
            {
            fprintf(fc,"%lf ",norm(POT(i,j)));
            }
            fprintf(fc,"\n");
    }
    fclose(fc);


    //method by his own
    DER = exp(-1I*DER*hbar*dt/(1.0*hbar));
    POT = exp(-1I*POT*hbar*dt/(2.0*hbar));

    for(int ts = 1; ts<N_Ts; ts++) //time steps
    {
        char nome1[64]; 
        snprintf(nome1, sizeof(char) * 64, "out/data%i.txt", ts - 1);

        FILE *fc = fopen(nome1,"w");

        for(int i = 0; i<N; i++)
        {
                for(int j = 0; j<N; j++)
                {
                    fprintf(fc,"%lf ",pow(norm(IC(i,j,0)),2));
                }
                fprintf(fc,"\n");
        }
        fclose(fc);

        cx_mat part1 = fft2(POT%IC.slice(0));
        cx_mat part2 = ifft2(DER%part1);
        cx_mat part3 = POT%part2;
        IC.slice(1) = part3;
        IC.slice(0) = IC.slice(1);
    }
    
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    return 0;
}









if(A>0=)
{
    x = 1;
}else
{
    x = -1;
}

//Is the samething of:
 
x = 1*(A>=0) - 1*(A<0);