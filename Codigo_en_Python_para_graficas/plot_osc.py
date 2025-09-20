import os
import numpy as np
import matplotlib.pyplot as plt

# ──────────────────────────────
# Configuración de rutas relativas
# ──────────────────────────────
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
data_file = os.path.join(base_dir, "Resultados_simulacion", "OSCILADOR", "Verlet", "V_1.txt")
output_dir = os.path.join(base_dir, "Graficas")
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, "trayectoria_oscilador.png")

# ──────────────────────────────
# Cargar datos numéricos ignorando filas con texto
# ──────────────────────────────
data = np.genfromtxt(data_file, usecols=(0,1), invalid_raise=False)

# Filtrar filas que tengan NaN (convertidas de texto)
data = data[~np.isnan(data).any(axis=1)]

tiempo = data[:, 0]
posicion = data[:, 1]

# ──────────────────────────────
# Crear gráfica
# ──────────────────────────────
plt.figure(figsize=(8,5))
plt.plot(tiempo, posicion, label="Oscilador armónico", color="blue", marker='o', markersize=2, linestyle='-')
plt.xlabel("Tiempo")
plt.ylabel("Posición")
plt.title("Trayectoria de un oscilador armónico")
plt.grid(True)
plt.legend()
plt.tight_layout()

# Ajustar límites de los ejes
plt.xlim(0, 100)
plt.ylim(1.5, -1.5)  # invertido para que 1.5 esté arriba y -1.5 abajo

# Guardar la figura
plt.savefig(output_file, dpi=300)
print(f"Gráfica guardada en: {output_file}")

# Mostrar la gráfica
plt.show()
