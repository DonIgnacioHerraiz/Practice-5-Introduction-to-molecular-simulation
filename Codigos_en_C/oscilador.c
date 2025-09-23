#include <stdio.h>
#include "integracion.h"
#include "funciones_oscilador.h"
#include "random.h"

int main(){
    int FLAG=2; // FLAG=1 para Verlet, FLAG=2 para Euler-Maruyama y FLAG=3 para Runge-Kutta
    
    inicializa_PR(123456); // Inicializa el generador con una semilla

    if(FLAG==1){
        double kb=1.0;
        double Temperatura=1.0;
        double alfa=10;
        int N=1;
        double dt=0.05;
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
        double kb=1.0;
        double Temperatura=1.0;
        double eta=10;
        int N=1;
        double h=0.001;
        double m=1.0;
        int pasos=10000;
        double x_0[N];
        double p_0[N];
        double K=1.0;

        x_0[0]=1.0; 
        p_0[0]=0.0;

        EulerMaruyama(K, kb, Temperatura, eta, N, h, m, pasos, Fuerza_euler, x_0, p_0);
    }else{
        if(FLAG==3){
            // Implementación de Runge-Kutta
        }else{
            printf("FLAG no válido. Use 1 para Verlet, 2 para Euler-Maruyama o 3 para Runge-Kutta.\n");
        }
    }
    }

}