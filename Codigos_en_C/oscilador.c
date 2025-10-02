#include <stdio.h>
#include "integracion.h"
#include "funciones_oscilador.h"
#include "random.h"

int main(){
    int FLAG=4; // FLAG=1 para Verlet, FLAG=2 para Euler-Maruyama y FLAG=3 para Runge-Kutta, FLAG=4 para Histogramas, FLAG 5 PARA EQUIPARTICION
    int FICHEROS=2;// FLAG=1 para Verlet, FLAG=2 para Euler-Maruyama y FLAG=3 para Runge-Kutta
    inicializa_PR(123456); // Inicializa el generador con una semilla

    if(FLAG==1){
        double kb=1.0;
        double Temperatura=1.0;
        double alfa=0.0;
        int N=1;
        double dt=0.05;
        double m=1.0;
        int pasos=10000000;
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
        double eta=0.0;
        int N=1;
        double h=0.001;
        double m=1.0;
        int pasos=10000000;
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
                double h=0.01;
                double m=1.0;
                int pasos=10000000;
                double x_0[N];
                double p_0[N];
                double K=1.0;
                int decision=0;
                x_0[0]=1.0; 
                p_0[0]=0.0;
                if (decision==0){
                    RungeKutta2(K, kb, Temperatura, eta, N, h, m, pasos, Fuerza_euler, x_0, p_0);
                }
                else{
                    RungeKutta2_histograma(1);
                }
            }else{
                if(FLAG==4){
                    switch(FICHEROS){
                        case 1:
                            generar_histogramas(
                            "Resultados_simulacion/OSCILADOR/VERLET",                      // carpeta de entrada
                            "Resultados_simulacion/OSCILADOR/VERLET/HISTOGRAMAS/VELOCIDADES", // carpeta salida velocidades
                            "Resultados_simulacion/OSCILADOR/VERLET/HISTOGRAMAS/POSICIONES",  // carpeta salida posiciones
                            "V",   // prefijo de los ficheros de entrada (R-K_0.txt, R-K_1.txt, …)
                            50     // número de bins del histograma
                            );break;
                        case 2:
                            generar_histogramas(
                            "Resultados_simulacion/OSCILADOR/EULER-MARUYAMA",                      // carpeta de entrada
                            "Resultados_simulacion/OSCILADOR/EULER-MARUYAMA/HISTOGRAMAS/VELOCIDADES", // carpeta salida velocidades
                            "Resultados_simulacion/OSCILADOR/EULER-MARUYAMA/HISTOGRAMAS/POSICIONES",  // carpeta salida posiciones
                            "E-M",   // prefijo de los ficheros de entrada (R-K_0.txt, R-K_1.txt, …)
                            50     // número de bins del histograma
                            );break;
                        case 3:
                            generar_histogramas(
                            "Resultados_simulacion/OSCILADOR/RUNGE-KUTTA",                      // carpeta de entrada
                            "Resultados_simulacion/OSCILADOR/RUNGE-KUTTA/HISTOGRAMAS/VELOCIDADES", // carpeta salida velocidades
                            "Resultados_simulacion/OSCILADOR/RUNGE-KUTTA/HISTOGRAMAS/POSICIONES",  // carpeta salida posiciones
                            "R-K",   // prefijo de los ficheros de entrada (R-K_0.txt, R-K_1.txt, …)
                            50     // número de bins del histograma
                            );break;}
            }else{
                if(FLAG==5){
                    switch(FICHEROS){
                        case 1:
                            procesar_archivos(
                            "Resultados_simulacion/OSCILADOR/VERLET",
                            "Resultados_simulacion/OSCILADOR/VERLET/EQUIPARTICION.txt",
                            "V"
                            );break;
                        case 2:
                            procesar_archivos(
                            "Resultados_simulacion/OSCILADOR/EULER-MARUYAMA",
                            "Resultados_simulacion/OSCILADOR/EULER-MARUYAMA/EQUIPARTICION.txt",
                            "E-M"
                            );break;
                        case 3:
                            procesar_archivos(
                            "Resultados_simulacion/OSCILADOR/RUNGE-KUTTA",
                            "Resultados_simulacion/OSCILADOR/RUNGE-KUTTA/EQUIPARTICION.txt",
                            "R-K"
                            );break;
                }
        }
    }
    return 0;
}}}}