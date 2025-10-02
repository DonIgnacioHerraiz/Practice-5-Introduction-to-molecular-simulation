#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <dirent.h>
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

void histogram (double *H, int N, double *data, int Thist, double *max, double *min, double *delta){
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
    normalization_factor = 1.0/(N)/(*delta);
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
 * @param filename Nombre del archivo donde se guardará el histograma.
*/
void histogram2D (double **H, int N, double *data1, double *data2, int Thist1, int Thist2, double *max1, double *min1, double *delta1,
    double *max2, double *min2, double *delta2, char* filename){
        int j1,j2;
        double normalization_factor;
        // Asignar memoria para las filas
        H = (double **)malloc(Thist1 * sizeof(double *));
        if (H == NULL) {
            printf("Error: No se pudo asignar memoria para las filas.\n");
            return;
        }

        // Asignar memoria para cada columna de cada fila
        for (int i = 0; i < Thist1; i++) {
            H[i] = (double *)calloc(Thist2, sizeof(double));
            if (H[i] == NULL) {
                printf("Error: No se pudo asignar memoria para la fila %d.\n", i);
                return;
            }
        }
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
        //Escribimos el histograma en un archivo
        FILE *file=fopen(filename, "w");
        if(file==NULL){
            printf("Error al abrir el archivo %s\n", filename);
            return;
        }
        for (int i=0; i<Thist1; i++){
            for (int k=0; k<Thist2; k++){
                fprintf(file, "%lf", H[i][k]);
            }
            fprintf(file, "\n");
        }
        fclose(file);
        // Liberamos la memoria dinámica
        for (int i = 0; i < Thist1; i++) {
            free(H[i]);
        }
    }

// Función principal de procesamiento
void generar_histogramas(const char *carpeta_in, 
                         const char *carpeta_out_vel, 
                         const char *carpeta_out_pos, 
                         const char *prefijo, 
                         int bins) {
    DIR *dir;
    struct dirent *ent;

    dir = opendir(carpeta_in);
    if (!dir) {
        perror("No se pudo abrir la carpeta de entrada");
        return;
    }

    while ((ent = readdir(dir)) != NULL) {
        if (strncmp(ent->d_name, prefijo, strlen(prefijo)) == 0 && strstr(ent->d_name, ".txt")) {
            // Ruta al archivo de entrada
            char ruta_in[512];
            snprintf(ruta_in, sizeof(ruta_in), "%s/%s", carpeta_in, ent->d_name);

            FILE *fin = fopen(ruta_in, "r");
            if (!fin) {
                perror("No se pudo abrir archivo de entrada");
                continue;
            }

            char linea[512];
            // saltar la primera línea
            if (!fgets(linea, sizeof(linea), fin)) {
                fclose(fin);
                continue;
            }

            // reservar memoria dinámica para datos
            int capacidad = 1000, N = 0;
            double *posiciones = malloc(capacidad * sizeof(double));
            double *velocidades = malloc(capacidad * sizeof(double));

            while (fgets(linea, sizeof(linea), fin)) {
                double c1, pos, vel, c4, c5, c6;
                if (sscanf(linea, "%lf %lf %lf %lf %lf %lf", &c1, &pos, &vel, &c4, &c5, &c6) == 6) {
                    if (N >= capacidad) {
                        capacidad *= 2;
                        posiciones = realloc(posiciones, capacidad * sizeof(double));
                        velocidades = realloc(velocidades, capacidad * sizeof(double));
                    }
                    posiciones[N] = pos;
                    velocidades[N] = vel;
                    N++;
                }
            }
            fclose(fin);

            // arrays para histogramas
            double *H_pos = malloc(bins * sizeof(double));
            double *H_vel = malloc(bins * sizeof(double));

            double min_p, max_p, delta_p;
            double min_v, max_v, delta_v;

            histogram(H_pos, N, posiciones, bins, &max_p, &min_p, &delta_p);
            histogram(H_vel, N, velocidades, bins, &max_v, &min_v, &delta_v);

            // crear archivos de salida
            char ruta_pos[512], ruta_vel[512];
            snprintf(ruta_pos, sizeof(ruta_pos), "%s/%s", carpeta_out_pos, ent->d_name);
            snprintf(ruta_vel, sizeof(ruta_vel), "%s/%s", carpeta_out_vel, ent->d_name);

            FILE *fpos = fopen(ruta_pos, "w");
            FILE *fvel = fopen(ruta_vel, "w");

            if (fpos) {
                for (int i = 0; i < bins; i++) {
                    double x_centro = min_p + (i + 0.5) * delta_p;
                    fprintf(fpos, "%lf %lf\n", x_centro, H_pos[i]);
                }
                fclose(fpos);
            }

            if (fvel) {
                for (int i = 0; i < bins; i++) {
                    double v_centro = min_v + (i + 0.5) * delta_v;
                    fprintf(fvel, "%lf %lf\n", v_centro, H_vel[i]);
                }
                fclose(fvel);
            }


            free(posiciones);
            free(velocidades);
            free(H_pos);
            free(H_vel);
        }
    }
    closedir(dir);
}