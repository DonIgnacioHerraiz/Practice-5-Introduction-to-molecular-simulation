import os
import numpy as np
import matplotlib.pyplot as plt

# ──────────────────────────────
# Configuración de rutas relativas
# ──────────────────────────────
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# Nombre del archivo
filename = "V_1.txt"
filename_no_ext = os.path.splitext(filename)[0]  

# Rutas de datos y parámetros
data_file = os.path.join(base_dir, "Resultados_simulacion", "OSCILADOR", "VERLET", filename)
param_file = os.path.join(base_dir, "PARAMETROS", "OSCILADOR", "VERLET", filename)
output_dir = os.path.join(base_dir, "Graficas", "OSCILADOR", "VERLET", "TRAYECTORIAS_ENERGIAS")
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
# Columnas: 0=tiempo, 3=E_potencial, 4=E_cinetica, 5=E_total
data = np.genfromtxt(data_file, usecols=(0,3,4,5), invalid_raise=False)
data = data[~np.isnan(data).any(axis=1)]  # filtrar filas con NaN

tiempo = data[:, 0]
E_potencial = data[:, 1]
E_cinetica = data[:, 2]
E_total = data[:, 3]

# ──────────────────────────────
# Crear gráfica solo con puntos
# ──────────────────────────────
plt.figure(figsize=(8,5))
plt.plot(tiempo, E_potencial, 'o', label="E_potencial", color="red", markersize=2)
plt.plot(tiempo, E_cinetica, 'o', label="E_cinetica", color="green", markersize=2)
plt.plot(tiempo, E_total, 'o', label="E_total", color="blue", markersize=2)

plt.xlabel("Tiempo")
plt.ylabel("Energía")
plt.title(f"Energías del oscilador ({filename})")
plt.grid(True)
plt.legend()
plt.tight_layout()

# Ajustar límites de los ejes
plt.xlim(0, x_max)
plt.ylim(min(E_potencial.min(), E_cinetica.min(), E_total.min())*1.1,
         max(E_potencial.max(), E_cinetica.max(), E_total.max())*1.1)

# Guardar la figura
plt.savefig(output_file, dpi=300)
print(f"Gráfica guardada en: {output_file}")

# Mostrar la gráfica
plt.show()
