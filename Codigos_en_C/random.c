#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#define PI 3.14159265358979323846



// Variables que hay que definir para Parisi-Rapuano
#define NormRANu (2.3283063671E-10F)
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;


//Esta función devuelve un numero aleatorio uniforme en (0,1)
double fran(void) // Cambiado a double
{
    double r; // Cambiado a double
    ig1 = (ind_ran - 24) & 255; // Asegura que el índice esté en [0, 255]
    ig2 = (ind_ran - 55) & 255; // Asegura que el índice esté en [0, 255]
    ig3 = (ind_ran - 61) & 255; // Asegura que el índice esté en [0, 255]
    irr[ind_ran] = irr[ig1] + irr[ig2];
    ir1 = (irr[ind_ran] ^ irr[ig3]);
    ind_ran = (ind_ran + 1) & 255; // Incrementa y asegura que esté en [0, 255]
    r = ir1 * (double)NormRANu; // Asegura que el cálculo sea en double
    return r;
}
void inicializa_PR(int SEMILLA)
{
    int INI,FACTOR,SUM,i;

    srand(SEMILLA);

    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;

    for(i=0;i<256;i++)
    {
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;
    }
    ind_ran=ig1=ig2=ig3=0;
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

//Función histograma 1D
/** 
 * @param H     Puntero al array donde se almacenará el histograma (debe tener tamaño Thist).
 * @param N     Número de datos en el array data.
 * @param data  Puntero al array de datos a analizar.
 * @param Thist Tamaño del histograma.
 * @param max   Puntero a variable donde se almacenará el valor máximo encontrado en data.
 * @param min   Puntero a variable donde se almacenará el valor mínimo encontrado en data.
 * @param delta Puntero a variable donde se almacenará el tamaño del intervalo del histograma.
*/

void histogram (int *H, int N, double *data, int Thist, double *max, double *min, double *delta){
    int j;
    double normalization_factor;
    // Inicializamos el histograma a cero
    for(int i=0; i<Thist; i++){
        H[i]=0;
    }
    // Encontramos el valor máximo y mínimo en los datos
    *max = data[0];
    *min = data[0]; 
    for(int i=1; i<N; i++){
        if(data[i]>*max) *max=data[i];
        if(data[i]<*min) *min=data[i];
    }
    // Calculamos el tamaño del intervalo
    *delta = (*max - *min) / Thist;
    // Llenamos el histograma
    for(int i=0; i<N; i++){
        j = (int)((data[i] - *min) / *delta);
        if(j >= 0 && j <= Thist) {
            if (j == Thist) j--; // Aseguramos que el valor máximo caiga en el último bin
            H[j]++;
        }
    }
    //Normalizamos el histograma
    normalization_factor = 1.0/(N * (*delta));
    for(int i=0; i<Thist; i++){
        H[i]=H[i]* normalization_factor;
    }
}

//Función histograma 2D
/** 
 * @param H     Puntero al array 2D donde se almacenará el histograma (debe tener tamaño Thist1 x Thist2).
 * @param N     Número de datos en los arrays data1 y data2.
 * @param data1 Puntero al array de datos a analizar en la primera dimensión.   
 * @param data2 Puntero al array de datos a analizar en la segunda dimensión.
 * @param Thist1 Tamaño del histograma en la primera dimensión.
 * @param Thist2 Tamaño del histograma en la segunda dimensión.
 * @param max1  Puntero a variable donde se almacenará el valor máximo encontrado en data1.
 * @param min1  Puntero a variable donde se almacenará el valor mínimo encontrado en data1.
 * @param delta1 Puntero a variable donde se almacenará el tamaño del intervalo del histograma en la primera dimensión.
 * @param max2  Puntero a variable donde se almacenará el valor máximo encontrado en data2.
 * @param min2  Puntero a variable donde se almacenará el valor mínimo encontrado en data2.
 * @param delta2 Puntero a variable donde se almacenará el tamaño del intervalo del histograma en la segunda dimensión.
*/
void histogram2D (int **H, int N, double *data1, double *data2, int Thist1, int Thist2, double *max1, double *min1, double *delta1,
    double *max2, double *min2, double *delta2){
        int j1,j2;
        double normalization_factor;
        // Inicializamos el histograma a cero
        for (int i=0; i<Thist1; i++){
            for(int k=0; k<Thist2; k++){
                H[i][k]=0;
            }
        }
        // Encontramos el valor máximo y mínimo en los datos
        *max1 = data1[0];
        *min1 = data1[0];
        *max2 = data2[0];
        *min2 = data2[0];
        for(int i=1; i<N; i++){
            if(data1[i]>*max1) *max1=data1[i];
            if(data1[i]<*min1) *min1=data1[i];
            if(data2[i]>*max2) *max2=data2[i];
            if(data2[i]<*min2) *min2=data2[i];
        }
        // Calculamos el tamaño del intervalo
        *delta1 = (*max1 - *min1) / Thist1;
        *delta2 = (*max2 - *min2) / Thist2;
        // Llenamos el histograma
        for(int i=0; i<N; i++){
            j1 = (int)((data1[i] - *min1) / *delta1);
            j2 = (int)((data2[i] - *min2) / *delta2);
            if(j1 >= 0 && j1 <= Thist1 && j2 >= 0 && j2 <= Thist2) {
                if (j1 == Thist1) j1--; // Aseguramos que el valor máximo caiga en el último bin
                if (j2 == Thist2) j2--; // Aseguramos que el valor máximo caiga en el último bin
                H[j1][j2]++;
            }
        }
        //Normalizamos el histograma
        normalization_factor = 1.0/(N * (*delta1) * (*delta2));
        for (int i=0; i<Thist1;i++){
            for(int k=0; k<Thist2; k++){
                H[i][k]=H[i][k]* normalization_factor;
            }
        }
    }

