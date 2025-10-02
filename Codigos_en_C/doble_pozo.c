#include <stdio.h>
#include "integracion.h"
#include "funciones_oscilador.h"
#include "random.h"

int main(){
    int FLAG=3; // FLAG=3 para Runge-Kutta, FLAG=4 para Histogramas, FLAG 5 PARA EQUIPARTICION
    inicializa_PR(123456); // Inicializa el generador con una semilla

    if(FLAG==3){
        // Implementación de Runge-Kutta
        double kb=0.2;
        double Temperatura=1.0;
        double eta=0.0;
        int N=1;
        double h=0.01;
        double m=1.0;
        int pasos=10000000;
        double x_0[N];
        double p_0[N];
        double A=1.0;
        int decision=0;
        x_0[0]=0.5; 
        p_0[0]=0.0;
        if (decision==0){
            RungeKutta2(A, kb, Temperatura, eta, N, h, m, pasos, Fuerza_euler, x_0, p_0);
        }
        else{
            RungeKutta2_histograma(1);
        }
        }else{
        if(FLAG==4){                
            generar_histogramas(
            "Resultados_simulacion/DOBLE_POZO",                      // carpeta de entrada
            "Resultados_simulacion/DOBLE_POZO/HISTOGRAMAS/VELOCIDADES", // carpeta salida velocidades
            "Resultados_simulacion/DOBLE_POZO/HISTOGRAMAS/POSICIONES",  // carpeta salida posiciones
            "R-K",   // prefijo de los ficheros de entrada (R-K_0.txt, R-K_1.txt, …)
            50     // número de bins del histograma
            );
        }else{
        if(FLAG==5){
            procesar_archivos(
            "Resultados_simulacion/DOBLE_POZO",
            "Resultados_simulacion/DOBLE_POZO/EQUIPARTICION.txt",
            "R-K"
            );
        }
        }
    }
    return 0;
}