
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