
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>

void Fuerza_verlet(int N, double x[], double F[], double K){
    for(int i=0; i<N; i++){
        F[i] = -K * x[i];
    }
}

void Fuerza_euler(int N, double x[], double p[], double F[], double K,double eta, double m){ //Emplearemos esta función tanto para Euler como para Runge-Kutta
    for(int i=0; i<N; i++){
        F[i] = -K * x[i] - eta * p[i]/m;
    }
}

double Energia_cinetica_instantanea(int N, double v[], double m){
    double K=0;
    for(int i=0;i<N;i++){
        K=K+0.5*m*v[i]*v[i];
    }
    return K;
}

double Energia_potencial_instantanea(int N, double x[], double m, double K){
    double U=0;
    for(int i=0;i<N;i++){
        U=U+0.5*K*x[i]*x[i];
    }
    return U;
}

double Energia_total_instantanea(int N, double x[], double v[], double m, double K){
    double E=0;
    E=Energia_cinetica_instantanea(N,v,m)+Energia_potencial_instantanea(N,x,m,K);
    return E;
}





