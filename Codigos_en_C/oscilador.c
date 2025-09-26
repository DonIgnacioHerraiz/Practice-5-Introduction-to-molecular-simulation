#include <stdio.h>
#include "integracion.h"
#include "funciones_oscilador.h"
#include "random.h"

int main(){
    int FLAG=4; // FLAG=1 para Verlet, FLAG=2 para Euler-Maruyama y FLAG=3 para Runge-Kutta, FLAG=4 para Histogramas, FLAG 5 PARA EQUIPARTICION
    
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
    }   else{
            if(FLAG==3){
              // Implementación de Runge-Kutta
                double kb=1.0;
                double Temperatura=1.0;
                double eta=0.0;
                int N=1;
                double h=0.001;
                double m=1.0;
                int pasos=10000;
                double x_0[N];
                double p_0[N];
                double K=1.0;

                x_0[0]=1.0; 
                p_0[0]=0.0;

                RungeKutta2(K, kb, Temperatura, eta, N, h, m, pasos, Fuerza_euler, x_0, p_0);
            }else{
                if(FLAG==4){
                    generar_histogramas(
        "Resultados_simulacion/OSCILADOR/RUNGE-KUTTA",                      // carpeta de entrada
        "Resultados_simulacion/OSCILADOR/RUNGE-KUTTA/HISTOGRAMAS/VELOCIDADES", // carpeta salida velocidades
        "Resultados_simulacion/OSCILADOR/RUNGE-KUTTA/HISTOGRAMAS/POSICIONES",  // carpeta salida posiciones
        "R-K",   // prefijo de los ficheros de entrada (R-K_0.txt, R-K_1.txt, …)
        50     // número de bins del histograma
    );
                
            }else{
                if(FLAG==5){
                        procesar_archivos(
        "Resultados_simulacion/OSCILADOR/RUNGE-KUTTA",
        "Resultados_simulacion/OSCILADOR/RUNGE-KUTTA/EQUIPARTICION.txt",
        "R-K"
    );
                }
        }
    }
    return 0;
}}}