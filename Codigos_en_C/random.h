#pragma once

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

#define PI 3.14159265358979323846

// Devuelve un número aleatorio uniforme en (0,1)
double fran(void);

// Inicializa el generador Parisi-Rapuano con una semilla
void inicializa_PR(int SEMILLA);

// Genera un número aleatorio N(0,1) con Box-Muller
double gaussian(void);

// Función histograma 1D
void histogram (double *H, int N, double *data, int Thist, double *max, double *min, double *delta);

void histogram2D (double **H, int N, double *data1, double *data2, int Thist1, int Thist2, double *max1, double *min1, double *delta1,
    double *max2, double *min2, double *delta2, char* filename);

// Función principal de procesamiento
void generar_histogramas(const char *carpeta_in, 
                         const char *carpeta_out_vel, 
                         const char *carpeta_out_pos, 
                         const char *prefijo, 
                         int bins);