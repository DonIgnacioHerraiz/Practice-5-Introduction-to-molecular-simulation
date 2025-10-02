#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>

// #define OSCILADOR
#define DOBLE_POZO

void Fuerza_verlet(int N, double x[], double F[], double K);

void Fuerza_euler(int N, double x[], double p[], double F[], double K,double eta, double m);

double Energia_cinetica_instantanea(int N, double v[], double m);

double Energia_potencial_instantanea(int N, double x[], double m, double K);

double Energia_total_instantanea(int N, double x[], double v[], double m, double K);

int es_archivo_valido(const char *nombre, const char *prefijo);

int extraer_indice(const char *nombre, const char *prefijo);

void procesar_archivos(const char *carpeta, const char *archivo_salida, const char *c_input);



