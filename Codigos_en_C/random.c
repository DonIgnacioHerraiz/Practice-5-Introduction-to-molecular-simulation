#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#define PI 3.14159265358979323846

//JOEL PON EN ESTA FUNCION PARISI-RAPUANO, LA QUE HE PUESTO AHORA ES PANCHITA
//Esta función devuelve un numero aleatorio uniforme en (0,1)
double fran(){
    return rand()/(RAND_MAX+1.0); 
}


// Función que genera un número aleatorio N(0,1) con Box-Muller
double gaussian() {
    double u1, u2;
    do {
        u1 = fran();
    } while (u1 <= 1e-10); // evitamos log(0)
    u2 = fran();

    return sqrt(-2.0 * log(u1)) * cos(2.0 *PI * u2);
}