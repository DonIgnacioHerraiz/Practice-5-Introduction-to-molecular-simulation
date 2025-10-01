#pragma once

#include "random.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
// Un paso del integrador de Verlet con ruido térmico
void un_paso_verlet(double betta[], double b, double a, int N,
                    double x_antiguo[], double x_nuevo[],
                    double v_antiguo[], double v_nuevo[],
                    double F_antiguo[], double F_nuevo[],
                    double dt, double m,
                    void (*Fuerza)(int, double[], double[], double),
                    double K);

// Integra la trayectoria completa con Verlet y guarda en archivo
void verlet_trayectoria(char* filename_input,
                        double kb, double Temperatura, double alfa,
                        int N, double dt, double m, int pasos,
                        void (*Fuerza)(int, double[], double[], double),
                        char* filename_output,
                        double x_0[], double v_0[], double K);

// Escribe archivo de parámetros de entrada para una simulación de Verlet
void escribe_input_verlet(double kb, double Temperatura, double alfa,
                          int N, double dt, double m, int pasos,
                          double x_0[], double v_0[],
                          char filename[], double K);

// Ejecuta simulación completa de Verlet (crea input y trayectoria)
void Verlet(double K, double kb, double Temperatura, double alfa,
            int N, double dt, double m, int pasos,
            void (*Fuerza)(int, double[], double[], double),
            double x_0[], double v_0[]);


// Un paso del integrador de Euler-Maruyama con ruido térmico
 void un_paso_euler(double gamma[], int N, double x_antiguo[], double x_nuevo[],
                    double p_antiguo[], double p_nuevo[],
                    double h, double m,
                    void (*Fuerza)(int, double[], double[], double [],double,double,double),
                    double Fuerzas[], double K, double eta);

// Integra la trayectoria completa con Euler-Maruyama y guarda en archivo
void euler_trayectoria(char* filename_input,
                        double kb, double Temperatura, double eta,
                        int N, double h, double m, int pasos,
                        void (*Fuerza)(int, double[], double[], double[],double,double,double),
                        char* filename_output,
                        double x_0[], double p_0[], double K);

// Escribe archivo de parámetros de entrada para una simulación de Euler-Maruyama
void escribe_input_euler(double kb, double Temperatura, double eta, 
                          int N, double h, double m, int pasos,
                          double x_0[], double p_0[],
                          char filename[], double K);

// Ejecuta simulación completa de Euler-Maruyama (crea input y trayectoria)
void EulerMaruyama(double K, double kb, double Temperatura, double eta,
            int N, double h, double m, int pasos,
            void (*Fuerza)(int, double[], double[], double [],double,double,double),
            double x_0[], double p_0[]);

// Un paso del integrador de Runge-Kutta de orden 2
 void paso_RungeKutta2(double Z[], int N, double x_antiguo[], double x_nuevo[],double p_antiguo[], double p_nuevo[], 
    double h, double m, void (*Fuerza)(int, double[],double[],double [],double,double,double), double K, double eta);

// Integra la trayectoria completa con Runge-Kutta y guarda en archivo
void rungeKutta_trayectoria(char* filename_input,double kb, double Temperatura,double eta,int N,double h, double m, int pasos,
     void (*Fuerza)(int, double[],double[],double[],double,double,double), 
     char*filename_output, double x_0[], double p_0[], double K);

// Escribe archivo de parámetros de entrada para una simulación de Runge-Kutta
void escribe_input_RungeKutta(double kb, double Temperatura, double eta, int N, double h, double m, int pasos,
      double x_0[], double p_0[], char filename[], double K);

// Ejecuta simulación completa de Runge-Kutta (crea input y trayectoria)
void RungeKutta2(double K,double kb, double Temperatura,double eta,int N,double h, double m, int pasos, 
    void (*Fuerza)(int, double[],double[],double [],double,double,double), double x_0[], double p_0[]);
            
//Funcion que me saque directamente las posiciones y momentos 
void RungeKutta2_pos_mom(double x_final[], double p_final[], int k, int num_data);

//Funcion que me cree el histograma de las posiciones y momentos finales con la funcion Histograma2D y me cree archivo con el propio histograma
void RungeKutta2_histograma(int k);