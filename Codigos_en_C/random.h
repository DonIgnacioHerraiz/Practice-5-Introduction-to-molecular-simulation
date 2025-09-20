#pragma once

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

#define PI 3.14159265358979323846

// Devuelve un número aleatorio uniforme en (0,1)
double fran(void);

// Genera un número aleatorio N(0,1) con Box-Muller
double gaussian(void);

