#include "random.h"
#include "funciones_oscilador.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>


/**
 * Realiza un paso en la integración del movimiento usando el método de Verlet.
 * @param betta       Array con términos aleatorios para el ruido térmico.
 * @param b           Coeficiente dependiente de alfa y dt.
 * @param a           Coeficiente dependiente de alfa y dt.
 * @param N           Número de partículas o elementos.
 * @param x_antiguo   Array con las posiciones en el paso de tiempo anterior.
 * @param x_nuevo     Array donde se almacenarán las nuevas posiciones.
 * @param v_antiguo   Array con las velocidades en el paso de tiempo anterior.
 * @param v_nuevo     Array donde se almacenarán las nuevas posiciones.
 * @param dt          Paso de tiempo.
 * @param m           Masa de las partículas (asumida igual para todas).
 * @param Fuerza      Puntero a función que calcula las fuerzas; recibe N, las posiciones las velocidades y un real.
 * @param F_antiguo   Array donde se almacenarán las fuerzas en el paso de tiempo anterior.
 * @param F_nuevo     Array donde se almacenarán las nuevas fuerzas.
 * @param K           Constante del oscilador armónico.
 */

 void un_paso_verlet(double betta[],double b, double a, int N, double x_antiguo[], double x_nuevo[],double v_antiguo[], double v_nuevo[],double F_antiguo[], double F_nuevo[], double dt, double m, void (*Fuerza)(int, double[],double[],double), double K) {
   
    Fuerza(N, x_antiguo, F_antiguo,K);


    for(int i=0; i<N; i++) {
        x_nuevo[i] = x_antiguo[i] + v_antiguo[i]*dt*b+F_antiguo[i]*dt*dt*b/(2*m)+b*dt*betta[i];
    }


    Fuerza(N, x_nuevo, F_nuevo,K);
   
    for(int i=0; i<N; i++) {
        v_nuevo[i] = a*v_antiguo[i] + (a*F_antiguo[i]+F_nuevo[i])*dt/(2*m)+b*betta[i]/m;
    }


}


/**
    * Realiza la integración del movimiento usando el método de Verlet durante un número dado de pasos.
    * Los resultados se guardan en un archivo.
    * @param kb              Constante de Boltzmann.
    * @param Temperatura     Temperatura del sistema.
    * @param alfa           Coeficiente de fricción.
    * @param N               Número de partículas o elementos.
    * @param dt              Paso de tiempo.  
    * @param m               Masa de las partículas (asumida igual para todas).
    * @param pasos           Número de pasos de tiempo a simular.
    * @param Fuerza          Puntero a función que calcula las fuerzas; recibe N y un array de posiciones.
    * @param filename_output Nombre del archivo donde se guardarán los resultados.
    * @param x_0            Array con las posiciones iniciales.
    * @param v_0            Array con las velocidades iniciales.
 */

void verlet_trayectoria(char* filename_input,double kb, double Temperatura,double alfa,int N,double dt, double m, int pasos, void (*Fuerza)(int, double[],double[],double), char*filename_output, double x_0[], double v_0[], double K) {
   
    FILE *archivo = fopen(filename_output, "w");
    if(archivo == NULL) {
        printf("Error al abrir el archivo %s\n", filename_output);
        return;
    }


    double a = (1.0 - alfa* dt / (2.0 * m)) / (1.0 + alfa * dt / (2.0 * m));
    double b = 1.0 / (1.0 + alfa * dt / (2.0 * m));
   
    fprintf(archivo, "%.6f %d\t%s\n", dt, pasos, filename_input);


    double x_antiguo[N];
    double x_nuevo[N];
    double v_antiguo[N];
    double v_nuevo[N];
    double F_antiguo[N];
    double F_nuevo[N];
    double betta[N];
    double Ek;
    double Ep;
    double Et;

    for(int i=0; i<N; i++) {
        x_antiguo[i] = x_0[i];
        v_antiguo[i] = v_0[i];
        x_nuevo[i] = 0;
        v_nuevo[i] = 0;
        F_antiguo[i] = 0;
        F_nuevo[i] = 0;
    }


    for(int paso=0; paso<pasos; paso++) {
       
        for(int i=0; i<N; i++) {
            betta[i] = gaussian()*sqrt(2*alfa*Temperatura*kb*dt);
        }


        un_paso_verlet(betta,b,a,N,x_antiguo,x_nuevo,v_antiguo,v_nuevo,F_antiguo,F_nuevo,dt,m,Fuerza,K);
       
        fprintf(archivo, "%.6f", paso * dt); // Tiempo en la primera columna


         for(int i = 0; i < N; i++) {
              fprintf(archivo, " %.6f", x_nuevo[i]); // Posiciones
        }
         for(int i = 0; i < N; i++) {
              fprintf(archivo, " %.6f", v_nuevo[i]); // Velocidades
        }
        Ek=Energia_cinetica_instantanea(N,v_nuevo,m);
        Ep=Energia_potencial_instantanea(N,x_nuevo,m,K);
        Et=Energia_total_instantanea(N,x_nuevo,v_nuevo,m,K);
        fprintf(archivo, " %.6f", Ek); // Cinetica
        fprintf(archivo, " %.6f", Ep); // Potencial
        fprintf(archivo, " %.6f", Et); // Total
        fprintf(archivo, "\n");
       
        for(int i=0; i<N; i++) {
            x_antiguo[i] = x_nuevo[i];
            v_antiguo[i] = v_nuevo[i];
        }
    }
    fclose(archivo);
}

/**
    * Escribe en un fichero los parámetros de la simulación de Verlet. Lo hace en la carpeta PARAMETROS\OSCILADOR con el formato V_i, siendo i el primer número natural tal que no existe un archivo con ese nombre en la carpeta.
    * @param kb              Constante de Boltzmann.
    * @param Temperatura     Temperatura del sistema.
    * @param alfa           Coeficiente de fricción.
    * @param N               Número de partículas o elementos.
    * @param dt              Paso de tiempo.  
    * @param m               Masa de las partículas (asumida igual para todas).
    * @param pasos           Número de pasos de tiempo a simular.
    * @param x_0            Array con las posiciones iniciales.
    * @param v_0            Array con las velocidades iniciales.
    * @param filename       Devuelve el nombre del archivo creado.
 */

void escribe_input_verlet(double kb, double Temperatura, double alfa, int N, double dt, double m, int pasos,
                          double x_0[], double v_0[], char filename[], double K) {
    const char* folder = "PARAMETROS/OSCILADOR/VERLET";
    FILE* file;
    int k = 0;


    // Buscamos el primer V_k.txt que no exista
    while (1) {
        snprintf(filename, 256, "%s/V_%d.txt", folder, k);
        file = fopen(filename, "r");
        if (file) {
            // El archivo existe, cerramos y probamos el siguiente
            fclose(file);
            k++;
        } else {
            // No existe, podemos usar este k
            break;
        }
    }


    // Creamos el archivo para escritura
    file = fopen(filename, "w");
    if (!file) {
        printf("No se pudo crear el archivo %s\n", filename);
        filename[0] = '\0'; // indicamos fallo
        return;
    }


    // Escribimos los parámetros
    fprintf(file, "K %g\n", K);
    fprintf(file, "kb %g\n", kb);
    fprintf(file, "Temperatura %g\n", Temperatura);
    fprintf(file, "alfa %g\n", alfa);
    fprintf(file, "N %d\n", N);
    fprintf(file, "dt %g\n", dt);
    fprintf(file, "m %g\n", m);
    fprintf(file, "pasos %d\n", pasos);


    // Escribimos los vectores x_0 y v_0
    for (int i = 0; i < N; i++) {
        fprintf(file, "x_0_%d %g\n", i, x_0[i]);
    }
    for (int i = 0; i < N; i++) {
        fprintf(file, "v_0_%d %g\n", i, v_0[i]);
    }


    fclose(file);


    printf("Archivo creado: %s\n", filename);
}

/**
    * Lleva a cabo la simulacion completa de Verlet. Los parametros estan en la carpeta PARAMETROS\OSCILADOR, mientras que la trayectoria está en Resultados_simulacion\OSCILADOR\VERLET
    * @param kb              Constante de Boltzmann.
    * @param Temperatura     Temperatura del sistema.
    * @param alfa           Coeficiente de fricción.
    * @param N               Número de partículas o elementos.
    * @param dt              Paso de tiempo.  
    * @param m               Masa de las partículas (asumida igual para todas).
    * @param pasos           Número de pasos de tiempo a simular.
    * @param x_0            Array con las posiciones iniciales.
    * @param v_0            Array con las velocidades iniciales.
    * @param Fuerza       Puntero a función que calcula las fuerzas; recibe N y un array de posiciones.
 */

void Verlet(double K,double kb, double Temperatura,double alfa,int N,double dt, double m, int pasos, void (*Fuerza)(int, double[],double[],double), double x_0[], double v_0[]){
    char filename_input[256];
    escribe_input_verlet(kb, Temperatura, alfa, N, dt, m, pasos, x_0, v_0, filename_input,K);
   
   
    const char* folder = "Resultados_simulacion/OSCILADOR/VERLET";
    char filename_output[256];
    FILE* file;
    int i = 0;


    // Buscamos el primer V_i.txt que no exista
    while (1) {
        snprintf(filename_output, sizeof(filename_output), "%s/V_%d.txt", folder, i);
        file = fopen(filename_output, "r");
        if (file) {
            // El archivo existe, cerramos y probamos el siguiente
            fclose(file);
            i++;
        } else {
            // No existe, podemos usar este i
            break;
        }
    }


    // Creamos el archivo
    file = fopen(filename_output, "w");
    if (!file) {
        printf("No se pudo crear el archivo %s\n", filename_output);
        return;
    }
    fclose(file);


    verlet_trayectoria(filename_input,kb,Temperatura,alfa,N,dt,m,pasos,Fuerza,filename_output,x_0,v_0,K);


}



//********************************EULER-MARUYAMA**********************************************



/**
 * Realiza un paso en la integración del movimiento usando el método de Verlet.
 * @param gamma       Array con términos aleatorios para el ruido térmico.
 * @param N           Número de partículas o elementos.
 * @param x_antiguo   Array con las posiciones en el paso de tiempo anterior.
 * @param x_nuevo     Array donde se almacenarán las nuevas posiciones.
 * @param p_antiguo   Array con las velocidades en el paso de tiempo anterior.
 * @param p_nuevo     Array donde se almacenarán las nuevas velocidades.
 * @param h          Paso de tiempo.
 * @param m           Masa de las partículas (asumida igual para todas).
 * @param Fuerza      Puntero a función que calcula las fuerzas; recibe N, las posiciones las velocidades y un real.
 * @param Fuerzas    Array donde se almacenarán las fuerzas en el paso de tiempo anterior.
 * @param K           Constante del oscilador armónico.
 * @param coef_damping Coeficiente de amortiguamiento.
 */

 void un_paso_euler(double gamma[], int N, double x_antiguo[], double x_nuevo[],double p_antiguo[], double p_nuevo[], double h, double m, void (*Fuerza)(int, double[],double[],double [],double,double,double), double Fuerzas[], double K, double eta) {

    Fuerza(N, x_antiguo,p_antiguo,Fuerzas,K,eta,m);


    for(int i=0; i<N; i++) {
        x_nuevo[i] = x_antiguo[i] + p_antiguo[i]/m*h;
    }
   
    for(int i=0; i<N; i++) {
        p_nuevo[i] = p_antiguo[i] + Fuerzas[i]*h + gamma[i];
    }


}


/**
    * Realiza la integración del movimiento usando el método de Verlet durante un número dado de pasos.
    * Los resultados se guardan en un archivo.
    * @param kb              Constante de Boltzmann.
    * @param Temperatura     Temperatura del sistema.
    * @param eta           Coeficiente de fricción.
    * @param N               Número de partículas o elementos.
    * @param h              Paso de tiempo.  
    * @param m               Masa de las partículas (asumida igual para todas).
    * @param pasos           Número de pasos de tiempo a simular.
    * @param Fuerza          Puntero a función que calcula las fuerzas; recibe N y un array de posiciones.
    * @param filename_output Nombre del archivo donde se guardarán los resultados.
    * @param x_0            Array con las posiciones iniciales.
    * @param p_0            Array con las velocidades iniciales.
 */

void euler_trayectoria(char* filename_input,double kb, double Temperatura,double eta,int N,double h, double m, int pasos, void (*Fuerza)(int, double[],double[],double[],double,double,double), char*filename_output, double x_0[], double p_0[], double K) {

    FILE *archivo = fopen(filename_output, "w");
    if(archivo == NULL) {
        printf("Error al abrir el archivo %s\n", filename_output);
        return;
    }
   
    fprintf(archivo, "%.6f %d\t%s\n", h, pasos, filename_input);


    double x_antiguo[N];
    double x_nuevo[N];
    double p_antiguo[N];
    double p_nuevo[N];
    double v_antiguo[N]; //estas velocidades son solo para calcular la energía cinética
    double v_nuevo[N];
    double Fuerzas[N];
    double gamma[N];
    double Ek;
    double Ep;
    double Et;

    for(int i=0; i<N; i++) {
        x_antiguo[i] = x_0[i];
        p_antiguo[i] = p_0[i];
        x_nuevo[i] = 0;
        p_nuevo[i] = 0;
        Fuerzas[i] = 0;
    }


    for(int paso=0; paso<pasos; paso++) {
       
        for(int i=0; i<N; i++) {
            gamma[i] = gaussian()*sqrt(2*eta*Temperatura*kb*h);
        }


        un_paso_euler(gamma,N,x_antiguo,x_nuevo,p_antiguo,p_nuevo,h,m,Fuerza,Fuerzas,K,eta);
       
        

        for(int i=0; i<N; i++) {
            v_antiguo[i]=p_antiguo[i]/m;
            v_nuevo[i]=p_nuevo[i]/m;
        }
        
        Ek=Energia_cinetica_instantanea(N,v_nuevo,m);
        Ep=Energia_potencial_instantanea(N,x_nuevo,m,K);
        Et=Energia_total_instantanea(N,x_nuevo,v_nuevo,m,K);

        if((h*paso*10-((int)(h*paso*10)))<1e-10) {   // Guardamos cada 0.1 unidades de tiempo
            fprintf(archivo, "%.6f", paso * h); // Tiempo en la primera columna


            for(int i = 0; i < N; i++) {
              fprintf(archivo, " %.6f", x_nuevo[i]); // Posiciones
            }
            for(int i = 0; i < N; i++) {
              fprintf(archivo, " %.6f", v_nuevo[i]); // Velocidades
            }


            fprintf(archivo, " %.6f", Ek); // Cinetica
            fprintf(archivo, " %.6f", Ep); // Potencial
            fprintf(archivo, " %.6f", Et); // Total
            fprintf(archivo, "\n");

        }

        for(int i=0; i<N; i++) {
            x_antiguo[i] = x_nuevo[i];
            p_antiguo[i] = p_nuevo[i];
        }
    }
    fclose(archivo);
}

/**
    * Escribe en un fichero los parámetros de la simulación de Verlet. Lo hace en la carpeta PARAMETROS\OSCILADOR con el formato V_i, siendo i el primer número natural tal que no existe un archivo con ese nombre en la carpeta.
    * @param kb              Constante de Boltzmann.
    * @param Temperatura     Temperatura del sistema.
    * @param eta           Coeficiente de fricción.
    * @param N               Número de partículas o elementos.
    * @param h              Paso de tiempo.  
    * @param m               Masa de las partículas (asumida igual para todas).
    * @param pasos           Número de pasos de tiempo a simular.
    * @param x_0            Array con las posiciones iniciales.
    * @param p_0            Array con los momentos iniciales.
    * @param filename       Devuelve el nombre del archivo creado.
 */

void escribe_input_euler(double kb, double Temperatura, double eta, int N, double h, double m, int pasos,
                          double x_0[], double p_0[], char filename[], double K) {
    const char* folder = "PARAMETROS/OSCILADOR/EULER-MARUYAMA";
    FILE* file;
    int k = 0;


    // Buscamos el primer E-M_k.txt que no exista
    while (1) {
        snprintf(filename, 256, "%s/E-M_%d.txt", folder, k);
        file = fopen(filename, "r");
        if (file) {
            // El archivo existe, cerramos y probamos el siguiente
            fclose(file);
            k++;
        } else {
            // No existe, podemos usar este k
            break;
        }
    }


    // Creamos el archivo para escritura
    file = fopen(filename, "w");
    if (!file) {
        printf("No se pudo crear el archivo %s\n", filename);
        filename[0] = '\0'; // indicamos fallo
        return;
    }


    // Escribimos los parámetros
    fprintf(file, "K %g\n", K);
    fprintf(file, "kb %g\n", kb);
    fprintf(file, "Temperatura %g\n", Temperatura);
    fprintf(file, "eta %g\n", eta);
    fprintf(file, "N %d\n", N);
    fprintf(file, "h %g\n", h);
    fprintf(file, "m %g\n", m);
    fprintf(file, "pasos %d\n", pasos);


    // Escribimos los vectores x_0 y p_0
    for (int i = 0; i < N; i++) {
        fprintf(file, "x_0_%d %g\n", i, x_0[i]);
    }
    for (int i = 0; i < N; i++) {
        fprintf(file, "p_0_%d %g\n", i, p_0[i]);
    }


    fclose(file);


    printf("Archivo creado: %s\n", filename);
}

/**
    * Lleva a cabo la simulacion completa de Verlet. Los parametros estan en la carpeta PARAMETROS\OSCILADOR, mientras que la trayectoria está en Resultados_simulacion\OSCILADOR\VERLET
    * @param kb              Constante de Boltzmann.
    * @param Temperatura     Temperatura del sistema.
    * @param eta           Coeficiente de fricción.
    * @param N               Número de partículas o elementos.
    * @param h              Paso de tiempo.  
    * @param m               Masa de las partículas (asumida igual para todas).
    * @param pasos           Número de pasos de tiempo a simular.
    * @param x_0            Array con las posiciones iniciales.
    * @param p_0            Array con los momentos iniciales.
    * @param Fuerza       Puntero a función que calcula las fuerzas; recibe N y un array de posiciones.
 */

void EulerMaruyama(double K,double kb, double Temperatura,double eta,int N,double h, double m, int pasos, void (*Fuerza)(int, double[],double[],double [],double,double,double), double x_0[], double p_0[]){
    char filename_input[256];
    escribe_input_euler(kb, Temperatura, eta, N, h, m, pasos, x_0, p_0, filename_input,K);


    const char* folder = "Resultados_simulacion/OSCILADOR/EULER-MARUYAMA";
    char filename_output[256];
    FILE* file;
    int i = 0;


    // Buscamos el primer E-M_i.txt que no exista
    while (1) {
        snprintf(filename_output, sizeof(filename_output), "%s/E-M_%d.txt", folder, i);
        file = fopen(filename_output, "r");
        if (file) {
            // El archivo existe, cerramos y probamos el siguiente
            fclose(file);
            i++;
        } else {
            // No existe, podemos usar este i
            break;
        }
    }


    // Creamos el archivo
    file = fopen(filename_output, "w");
    if (!file) {
        printf("No se pudo crear el archivo %s\n", filename_output);
        return;
    }
    fclose(file);


    euler_trayectoria(filename_input,kb,Temperatura,eta,N,h,m,pasos,Fuerza,filename_output,x_0,p_0,K);


}

//****************RUNGE KUTTA ESTOCASTICO DE ORDEN 2 ******************************
/**
 * Realiza un paso en la integración del movimiento usando el método de Verlet.
 * @param Z           Array con términos aleatorios para el ruido térmico.
 * @param N           Número de partículas o elementos.
 * @param x_antiguo   Array con las posiciones en el paso de tiempo anterior.
 * @param x_nuevo     Array donde se almacenarán las nuevas posiciones.
 * @param p_antiguo   Array con las velocidades en el paso de tiempo anterior.
 * @param p_nuevo     Array donde se almacenarán las nuevas velocidades.
 * @param h          Paso de tiempo.
 * @param m           Masa de las partículas (asumida igual para todas).
 * @param Fuerza      Puntero a función que calcula las fuerzas; recibe N, las posiciones las velocidades y un real.
 * @param K           Constante del oscilador armónico.
 * @param eta         Coeficiente de amortiguamiento.
 */
void paso_RungeKutta2(double Z[], int N, double x_antiguo[], double x_nuevo[],double p_antiguo[], double p_nuevo[], double h, double m, void (*Fuerza)(int, double[],double[],double [],double,double,double), double K, double eta) {
    double x_intermedio[N];
    double p_intermedio[N];
    //Funciones primera etapa
    double fx1[N];
    double gp1[N];
    //Funciones segunda etapa
    double fx2[N];
    double gp2[N];
    //Calculo de la fuerza para la primera etapa
    for (int i=0; i<N; i++){
        p_intermedio[i]=Z[i]+p_antiguo[i];// pn+Z
    }
    Fuerza(N, x_antiguo,p_intermedio,gp1,K,eta,m);//g(xn, pn+Z)
    //Primera etapa
    for(int i=0; i<N; i++) {
        fx1[i] = p_intermedio[i]/m;
    }
    //Calculo de la fuerza para la segunda etapa
    for (int i=0; i<N; i++){
        x_intermedio[i]=x_antiguo[i]+h*fx1[i];// xn+h*f(xn, pn+Z)
        p_intermedio[i]=p_antiguo[i]+h*gp1[i];// pn+h*g(xn, pn+Z)
    }
    Fuerza(N, x_intermedio,p_intermedio,gp2,K,eta,m);//g(xn+h*f(xn, pn+Z), pn+h*g(xn, pn+Z))
    //Segunda etapa
    for(int i=0; i<N; i++) {
        fx2[i] = p_intermedio[i]/m;
    }
    //Actualización de las variables
    for(int i=0; i<N; i++) {
        x_nuevo[i] = x_antiguo[i] + (h/2.0)*(fx1[i]+fx2[i]);
        p_nuevo[i] = p_antiguo[i] + (h/2.0)*(gp1[i]+gp2[i]) + Z[i];
    }
}

/**
    * Realiza la integración del movimiento usando el método de Verlet durante un número dado de pasos.
    * Los resultados se guardan en un archivo.
    * @param kb              Constante de Boltzmann.
    * @param Temperatura     Temperatura del sistema.
    * @param eta           Coeficiente de fricción.
    * @param N               Número de partículas o elementos.
    * @param h              Paso de tiempo.  
    * @param m               Masa de las partículas (asumida igual para todas).
    * @param pasos           Número de pasos de tiempo a simular.
    * @param Fuerza          Puntero a función que calcula las fuerzas; recibe N y un array de posiciones.
    * @param filename_output Nombre del archivo donde se guardarán los resultados.
    * @param x_0            Array con las posiciones iniciales.
    * @param p_0            Array con las velocidades iniciales.
    * @param K           Constante del oscilador armónico.
 */

void rungeKutta_trayectoria(char* filename_input,double kb, double Temperatura,double eta,int N,double h, double m, int pasos, void (*Fuerza)(int, double[],double[],double[],double,double,double),  char*filename_output, double x_0[], double p_0[], double K) {

    FILE *archivo = fopen(filename_output, "w");
    if(archivo == NULL) {
        printf("Error al abrir el archivo %s\n", filename_output);
        return;
    }
   
    fprintf(archivo, "%.6f %d\t%s\n", h, pasos, filename_input);


    double x_antiguo[N];
    double x_nuevo[N];
    double p_antiguo[N];
    double p_nuevo[N];
    double v_antiguo[N]; //estas velocidades son solo para calcular la energía cinética
    double v_nuevo[N];
    double gamma[N];
    double Ek;
    double Ep;
    double Et;

    //Inicialización de arrays a cero para los que serán editados y a sus valores iniciales para los que no
    for(int i=0; i<N; i++) {
        x_antiguo[i] = x_0[i];
        p_antiguo[i] = p_0[i];
        x_nuevo[i] = 0;
        p_nuevo[i] = 0;
    }


    for(int paso=0; paso<pasos; paso++) {
        for(int i=0; i<N; i++) {
            gamma[i] = gaussian()*sqrt(2*eta*Temperatura*kb*h);
        }


        paso_RungeKutta2(gamma,N,x_antiguo,x_nuevo,p_antiguo,p_nuevo,h,m,Fuerza,K,eta);
       
        fprintf(archivo, "%.6f", paso * h); // Tiempo en la primera columna


         for(int i = 0; i < N; i++) {
              fprintf(archivo, " %.6f", x_nuevo[i]); // Posiciones
        }
         for(int i = 0; i < N; i++) {
              fprintf(archivo, " %.6f", p_nuevo[i]); // Velocidades
        }

        for(int i=0; i<N; i++) {
            v_antiguo[i]=p_antiguo[i]/m;
            v_nuevo[i]=p_nuevo[i]/m;
        }
        
        Ek=Energia_cinetica_instantanea(N,v_nuevo,m);
        Ep=Energia_potencial_instantanea(N,x_nuevo,m,K);
        Et=Energia_total_instantanea(N,x_nuevo,v_nuevo,m,K);
        fprintf(archivo, " %.6f", Ek); // Cinetica
        fprintf(archivo, " %.6f", Ep); // Potencial
        fprintf(archivo, " %.6f", Et); // Total
        fprintf(archivo, "\n");
       
        for(int i=0; i<N; i++) {
            x_antiguo[i] = x_nuevo[i];
            p_antiguo[i] = p_nuevo[i];
        }
    }
    fclose(archivo);
}

/**
    * Escribe en un fichero los parámetros de la simulación de Verlet. Lo hace en la carpeta PARAMETROS\OSCILADOR con el formato V_i, siendo i el primer número natural tal que no existe un archivo con ese nombre en la carpeta.
    * @param kb              Constante de Boltzmann.
    * @param Temperatura     Temperatura del sistema.
    * @param eta           Coeficiente de fricción.
    * @param N               Número de partículas o elementos.
    * @param h              Paso de tiempo.  
    * @param m               Masa de las partículas (asumida igual para todas).
    * @param pasos           Número de pasos de tiempo a simular.
    * @param x_0            Array con las posiciones iniciales.
    * @param p_0            Array con los momentos iniciales.
    * @param filename       Devuelve el nombre del archivo creado.
    * @param K           Constante del oscilador armónico.
 */

void escribe_input_RungeKutta(double kb, double Temperatura, double eta, int N, double h, double m, int pasos,double x_0[], double p_0[], char filename[], double K) {
    const char* folder = "PARAMETROS/OSCILADOR/RUNGE-KUTTA";
    FILE* file;
    int k = 0;


    // Buscamos el primer R-K_k.txt que no exista
    while (1) {
        snprintf(filename, 256, "%s/R-K_%d.txt", folder, k);
        file = fopen(filename, "r");
        if (file) {
            // El archivo existe, cerramos y probamos el siguiente
            fclose(file);
            k++;
        } else {
            // No existe, podemos usar este k
            break;
        }
    }


    // Creamos el archivo para escritura
    file = fopen(filename, "w");
    if (!file) {
        printf("No se pudo crear el archivo %s\n", filename);
        filename[0] = '\0'; // indicamos fallo
        return;
    }


    // Escribimos los parámetros
    fprintf(file, "K %g\n", K);
    fprintf(file, "kb %g\n", kb);
    fprintf(file, "Temperatura %g\n", Temperatura);
    fprintf(file, "eta %g\n", eta);
    fprintf(file, "N %d\n", N);
    fprintf(file, "h %g\n", h);
    fprintf(file, "m %g\n", m);
    fprintf(file, "pasos %d\n", pasos);


    // Escribimos los vectores x_0 y p_0
    for (int i = 0; i < N; i++) {
        fprintf(file, "x_0_%d %g\n", i, x_0[i]);
    }
    for (int i = 0; i < N; i++) {
        fprintf(file, "p_0_%d %g\n", i, p_0[i]);
    }


    fclose(file);


    printf("Archivo creado: %s\n", filename);
}

/**
    * Lleva a cabo la simulacion completa de Verlet. Los parametros estan en la carpeta PARAMETROS\OSCILADOR, mientras que la trayectoria está en Resultados_simulacion\OSCILADOR\VERLET
    * @param kb              Constante de Boltzmann.
    * @param Temperatura     Temperatura del sistema.
    * @param eta           Coeficiente de fricción.
    * @param N               Número de partículas o elementos.
    * @param h              Paso de tiempo.  
    * @param m               Masa de las partículas (asumida igual para todas).
    * @param pasos           Número de pasos de tiempo a simular.
    * @param x_0            Array con las posiciones iniciales.
    * @param p_0            Array con los momentos iniciales.
    * @param Fuerza       Puntero a función que calcula las fuerzas; recibe N y un array de posiciones.
    * @param K           Constante del oscilador armónico.
 */

void RungeKutta2(double K,double kb, double Temperatura,double eta,int N,double h, double m, int pasos, 
    void (*Fuerza)(int, double[],double[],double [],double,double,double), double x_0[], double p_0[]){
    char filename_input[256];
    escribe_input_RungeKutta(kb, Temperatura, eta, N, h, m, pasos, x_0, p_0, filename_input,K);


    const char* folder = "Resultados_simulacion/OSCILADOR/RUNGE-KUTTA";
    char filename_output[256];
    FILE* file;
    int i = 0;


    // Buscamos el primer R-K_i.txt que no exista
    while (1) {
        snprintf(filename_output, sizeof(filename_output), "%s/R-K_%d.txt", folder, i);
        file = fopen(filename_output, "r");
        if (file) {
            // El archivo existe, cerramos y probamos el siguiente
            fclose(file);
            i++;
        } else {
            // No existe, podemos usar este i
            break;
        }
    }


    // Creamos el archivo
    file = fopen(filename_output, "w");
    if (!file) {
        printf("No se pudo crear el archivo %s\n", filename_output);
        return;
    }
    fclose(file);


    rungeKutta_trayectoria(filename_input,kb,Temperatura,eta,N,h,m,pasos,Fuerza,filename_output,x_0,p_0,K);


}
/** 
    @param x_final: array donde se guardan las posiciones finales
    @param p_final: array donde se guardan los momentos finales
    @param k: número del archivo R-K_k.txt del que se quieren obtener los datos
*/

void RungeKutta2_pos_mom(double x_final[], double p_final[], int k, int num_data){
    FILE *file;
    const char* folder = "Resultados_simulacion/OSCILADOR/RUNGE-KUTTA";
    char filename[256];
    // Buscamos el primer R-K_k.txt que queramos dentro de los que hayamos creado previamente
    snprintf(filename, 256, "%s/R-K_%d.txt", folder,k);
    file = fopen(filename, "r");
    if (!file)
    {
        printf("No se pudo abrir el archivo %s\n", filename);
        return;
    }
    //leer archivo para obtener las posiciones y momentos finales
    int contador_lineas=0;
    while (!feof(file))
    {
        double a1,a2,a3,a4,a5;
        if (contador_lineas==1){
            fscanf(file, "%lf %lf %lf %lf %lf",&a1,&a2,&a3,&a4,&a5); // Saltar el tiempo
            x_final[contador_lineas-1]=a2;
            p_final[contador_lineas-1]=a3;  
        }

        contador_lineas++;

    }    
    num_data=contador_lineas; 
    fclose(file);
} 
/**
 * @param k: número del archivo R-K_k.txt del que se quieren obtener los datos
 * @param num_data: número de datos que hay en el archivo
 * @param num_bins: número de bins que se quieren en el histograma
 */

//Funcion que me cree el histograma de las posiciones y momentos finales con la funcion Histograma2D y me cree archivo con el propio histograma
void RungeKutta2_histograma(int k){
    int num_data;
    int num_bins=50;
    double **H;//Donde guardaremos el histograma 2D
    double max1, min1, delta1;
    double max2, min2, delta2;
    char filename[256];
    double x_final[num_data];
    double p_final[num_data];
    //Leemos el archivo R-K_k.txt para obtener las posiciones y momentos finales
    RungeKutta2_pos_mom(x_final, p_final, k, num_data);
    //Creamos el histograma 2D
    const char* folder = "Resultados_simulacion/OSCILADOR/RUNGE-KUTTA/HISTOGRAMA/2D";
    snprintf(filename, 256, "%s/Histograma_R-K-2D_%d.txt", folder,k);
    histogram2D(H,num_data,x_final,p_final,num_bins,num_bins,&max1,&min1,&delta1,&max2,&min2,&delta2,filename);
}



