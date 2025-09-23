#include <stdio.h>
#include "random.h"
int main() {
    inicializa_PR(123456); // Inicializa el generador con una semilla
    FILE *f = fopen("TESTS/Box-Muller/datos_box_muller.txt", "w");
    int contador = 0;
    double media = 0.0;
    double varianza = 0.0;

    // Generar y mostrar 10 números aleatorios con la función Box-Muller
    for (int i = 0; i < 100000; i++) {
        double num = gaussian();
        fprintf(f, "%f\n", num);
        media += num;
        contador++;
    }

    media /= contador;
    fclose(f);
    // Calcular la varianza
    f = fopen("TESTS/Box-Muller/datos_box_muller.txt", "r");
    for (int i = 0; i < 100000; i++) {
        double num;
        fscanf(f, "%lf", &num);
        varianza += (num - media) * (num - media);
    }
    varianza /= (contador-1);

    printf("Media: %f\n", media);
    printf("Varianza: %f\n", varianza);
    


    fclose(f);
    return 0;
}