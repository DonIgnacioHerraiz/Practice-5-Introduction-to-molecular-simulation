import os
import numpy as np
import matplotlib.pyplot as plt

# ──────────────────────────────
# Configuración de rutas relativas
# ──────────────────────────────
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Nombre del archivo
filename = "V_4.txt"
filename_no_ext = os.path.splitext(filename)[0]  

# Rutas de datos y parámetros
data_file = os.path.join(base_dir, "Resultados_simulacion", "OSCILADOR", "Verlet", filename)
param_file = os.path.join(base_dir, "PARAMETROS", "OSCILADOR", filename)
output_dir = os.path.join(base_dir, "Graficas", "OSCILADOR", "TRAYECTORIAS")
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, f"{filename_no_ext}.png")  # guardará V_1.png

# ──────────────────────────────
# Leer parámetros del archivo
# ──────────────────────────────
params = {}
with open(param_file, 'r') as f:
    for line in f:
        key, value = line.split()
        params[key] = float(value)

dt = params["dt"]
pasos = int(params["pasos"])
x_max = pasos * dt  # límite superior para el eje X

# ──────────────────────────────
# Cargar datos numéricos ignorando filas con texto
# ──────────────────────────────
data = np.genfromtxt(data_file, usecols=(0,1), invalid_raise=False)
data = data[~np.isnan(data).any(axis=1)]  # filtrar filas con NaN

tiempo = data[:, 0]
posicion = data[:, 1]

# ──────────────────────────────
# Crear gráfica
# ──────────────────────────────
plt.figure(figsize=(8,5))
plt.plot(tiempo, posicion, label=f"Oscilador armónico ({filename})", color="blue", marker='o', markersize=2, linestyle='-')
plt.xlabel("Tiempo")
plt.ylabel("Posición")
plt.title(f"Trayectoria de un oscilador armónico ({filename})")
plt.grid(True)
plt.legend()
plt.tight_layout()

# Ajustar límites de los ejes
plt.xlim(0, x_max)      # eje X basado en pasos*dt
plt.ylim(3.5, -3.5)     # invertido para que 1.5 esté arriba y -1.5 abajo

# Guardar la figura
plt.savefig(output_file, dpi=300)
print(f"Gráfica guardada en: {output_file}")

# Mostrar la gráfica
plt.show()
