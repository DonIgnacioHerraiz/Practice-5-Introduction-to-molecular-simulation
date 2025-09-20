#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>

void Fuerza_verlet(int N, double x[], double F[], double K);

double Energia_cinetica_instantanea(int N, double v[], double m);

double Energia_potencial_instantanea(int N, double x[], double m, double K);

double Energia_total_instantanea(int N, double x[], double v[], double m, double K);

