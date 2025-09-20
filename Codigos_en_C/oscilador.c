#include <stdio.h>
#include "integracion.h"
#include "funciones_oscilador.h"

int main(){
    int FLAG=1; // FLAG=1 para Verlet, FLAG=2 para Euler-Maruyama y FLAG=3 para Runge-Kutta

    if(FLAG==1){
        double kb=1.0;
        double Temperatura=1.0;
        double alfa=0;
        int N=1;
        double dt=0.01;
        double m=1.0;
        int pasos=10000;
        double x_0[N];
        double v_0[N];
        double K=1.0;

        x_0[0]=1.0; 
        v_0[0]=0.0;

        Verlet(K, kb, Temperatura, alfa, N, dt, m, pasos, Fuerza_verlet, x_0, v_0);
    }else{
        if(FLAG==2){
            // Implementación de Euler-Maruyama
    }else{
        if(FLAG==3){
            // Implementación de Runge-Kutta
        }else{
            printf("FLAG no válido. Use 1 para Verlet, 2 para Euler-Maruyama o 3 para Runge-Kutta.\n");
        }
    }
    }

}