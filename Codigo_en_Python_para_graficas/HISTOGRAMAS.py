import os
import matplotlib.pyplot as plt
import numpy as np

def plot_histograms(prefijo, carpeta_entrada, carpeta_salida):
    tipos = ['POSICIONES', 'VELOCIDADES']

    for tipo in tipos:
        carpeta_in = os.path.join(carpeta_entrada, tipo)
        carpeta_out = os.path.join(carpeta_salida, tipo)

        os.makedirs(carpeta_out, exist_ok=True)  # crear carpeta si no existe

        # listar todos los archivos que empiecen con prefijo y terminen en .txt
        archivos = sorted([f for f in os.listdir(carpeta_in) if f.startswith(prefijo) and f.endswith('.txt')])

        for archivo in archivos:
            ruta_in = os.path.join(carpeta_in, archivo)

            # leer datos
            x_vals = []
            y_vals = []
            with open(ruta_in, 'r') as f:
                for linea in f:
                    parts = linea.strip().split()
                    if len(parts) == 2:
                        x, y = map(float, parts)
                        x_vals.append(x)
                        y_vals.append(y)

            x_vals = np.array(x_vals)
            y_vals = np.array(y_vals)

            # calcular anchos de barra dinámicos
            if len(x_vals) > 1:
                # distancia entre cada centro de bin
                widths = np.diff(x_vals)
                # para el último bin, usamos el mismo ancho que el penúltimo
                widths = np.append(widths, widths[-1])
            else:
                widths = np.array([1.0])  # un único bin

            # crear figura
            plt.figure(figsize=(8,5))
            plt.bar(x_vals, y_vals, width=widths, align='center', color='skyblue', edgecolor='black')
            plt.xlabel('x')
            plt.ylabel('Frecuencia')
            plt.title(f'Histograma {tipo} - {archivo}')
            plt.tight_layout()

            # guardar figura
            ruta_out = os.path.join(carpeta_out, archivo.replace('.txt', '.png'))
            plt.savefig(ruta_out)
            plt.close()
            print(f"Guardado: {ruta_out}")

# --- Ejemplo de uso ---
prefijo = "R-K"
carpeta_entrada = r"Resultados_simulacion\OSCILADOR\RUNGE-KUTTA\HISTOGRAMAS"
carpeta_salida  = r"Graficas\OSCILADOR\RUNGE-KUTTA\HISTOGRAMAS"

plot_histograms(prefijo, carpeta_entrada, carpeta_salida)
