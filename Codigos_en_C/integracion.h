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