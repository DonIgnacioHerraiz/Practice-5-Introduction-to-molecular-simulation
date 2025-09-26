
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>

void Fuerza_verlet(int N, double x[], double F[], double K){
    for(int i=0; i<N; i++){
        F[i] = -K * x[i];
    }
}

void Fuerza_euler(int N, double x[], double p[], double F[], double K,double eta, double m){ //Emplearemos esta función tanto para Euler como para Runge-Kutta
    for(int i=0; i<N; i++){
        F[i] = -K * x[i] - eta * p[i]/m;
    }
}

double Energia_cinetica_instantanea(int N, double v[], double m){
    double K=0;
    for(int i=0;i<N;i++){
        K=K+0.5*m*v[i]*v[i];
    }
    return K;
}

double Energia_potencial_instantanea(int N, double x[], double m, double K){
    double U=0;
    for(int i=0;i<N;i++){
        U=U+0.5*K*x[i]*x[i];
    }
    return U;
}

double Energia_total_instantanea(int N, double x[], double v[], double m, double K){
    double E=0;
    E=Energia_cinetica_instantanea(N,v,m)+Energia_potencial_instantanea(N,x,m,K);
    return E;
}

int es_archivo_valido(const char *nombre, const char *prefijo) {
    size_t len = strlen(prefijo);
    // comprobar si empieza con prefijo + "_"
    if (strncmp(nombre, prefijo, len) == 0 && nombre[len] == '_' && strstr(nombre, ".txt")) {
        return 1;
    }
    return 0;
}

int extraer_indice(const char *nombre, const char *prefijo) {
    int idx;
    char formato[64];
    snprintf(formato, sizeof(formato), "%s_%%d.txt", prefijo);
    sscanf(nombre, formato, &idx);
    return idx;
}

void procesar_archivos(const char *carpeta, const char *archivo_salida, const char *c_input) {
    DIR *dir;
    struct dirent *ent;
    FILE *fout;
    int LINEA_MAX = 512;

    dir = opendir(carpeta);
    if (!dir) {
        perror("No se pudo abrir la carpeta");
        return;
    }

    fout = fopen(archivo_salida, "w"); // sobrescribe el archivo si existe
    if (!fout) {
        perror("No se pudo abrir archivo de salida");
        closedir(dir);
        return;
    }

    char *archivos[1024];
    int n = 0;

    while ((ent = readdir(dir)) != NULL) {
        if (es_archivo_valido(ent->d_name, c_input)) {
            archivos[n] = strdup(ent->d_name);
            n++;
        }
    }
    closedir(dir);

    // ordenar archivos por índice k
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (extraer_indice(archivos[i], c_input) > extraer_indice(archivos[j], c_input)) {
                char *tmp = archivos[i];
                archivos[i] = archivos[j];
                archivos[j] = tmp;
            }
        }
    }

    // procesar archivos
    for (int i = 0; i < n; i++) {
        char ruta[512];
        snprintf(ruta, sizeof(ruta), "%s/%s", carpeta, archivos[i]);

        FILE *fin = fopen(ruta, "r");
        if (!fin) {
            perror("No se pudo abrir archivo de entrada");
            free(archivos[i]);
            continue;
        }

        char linea[LINEA_MAX];
        if (!fgets(linea, LINEA_MAX, fin)) { // saltar primera fila
            fclose(fin);
            free(archivos[i]);
            continue;
        }

        double col4, col5, col6;
        double suma4 = 0.0, suma5 = 0.0, suma6 = 0.0;
        int count = 0;

        while (fgets(linea, LINEA_MAX, fin)) {
            double c1, c2, c3;
            if (sscanf(linea, "%lf %lf %lf %lf %lf %lf",
                       &c1, &c2, &c3, &col4, &col5, &col6) == 6) {
                suma4 += col4;
                suma5 += col5;
                suma6 += col6;
                count++;
            }
        }
        fclose(fin);

        if (count > 0) {
            fprintf(fout, "------- %s -------\n\n", archivos[i]);
            fprintf(fout, "Ek: %lf\n", suma4 / count);
            fprintf(fout, "Ep: %lf\n", suma5 / count);
            fprintf(fout, "E: %lf\n\n", suma6 / count);
        }

        free(archivos[i]);
    }

    fclose(fout);
}






