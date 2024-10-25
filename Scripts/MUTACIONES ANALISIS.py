#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 23:11:35 2024

@author: ics2
"""

import matplotlib.pyplot as plt
import pandas as pd

# Cargar el archivo CSV de mutaciones
df = pd.read_csv('/home/ics2/mutaciones_resultado.csv')

# Convertir las regiones a posiciones de aminoácidos
REGIONES = [
    {"start": int((r["start"] - 1) / 3) + 1, "end": int(r["end"] / 3), "gene": r["gene"], "product": r["product"]}
    for r in [
        {"start": 95, "end": 436, "gene": "poly_", "product": "ancC:"},
        {"start": 95, "end": 394, "gene": "poly_", "product": "C:"},
        {"start": 437, "end": 934, "gene": "poly_", "product": "prM:"},
        {"start": 437, "end": 709, "gene": "poly_", "product": "pr:"},
        {"start": 710, "end": 934, "gene": "poly_", "product": "M:"},
        {"start": 935, "end": 2413, "gene": "poly_", "product": "E:"},
        {"start": 2414, "end": 3469, "gene": "poly_", "product": "NS1:"},
        {"start": 3470, "end": 4123, "gene": "poly_", "product": "NS2A:"},
        {"start": 4124, "end": 4513, "gene": "poly_", "product": "NS2B:"},
        {"start": 4514, "end": 6370, "gene": "poly_", "product": "NS3:"},
        {"start": 6371, "end": 6751, "gene": "poly_", "product": "NS4A:"},
        {"start": 6752, "end": 6820, "gene": "poly_", "product": "2K:"},
        {"start": 6821, "end": 7564, "gene": "poly_", "product": "NS4B:"},
        {"start": 7565, "end": 10264, "gene": "poly_", "product": "NS5:"},
    ]
]

# Extraer año
df['año'] = df['nombre'].str.extract(r'/(\d{4})$')

# Convertir la columna 'año' a enteros para facilitar la comparación
df['año'] = df['año'].astype(int)

# Dividir las mutaciones en dos grupos según el año
df_mayor_2021 = df[df['año'] > 2021]
df_menor_2021 = df[df['año'] <= 2021]

def extract_mutations_aa(mutations_str):
    mutations = mutations_str.split(';')
    valid_mutations = []
    for m in mutations:
        parts = m.split(':')
        if len(parts) == 2:  # Verifica que haya exactamente dos partes
            valid_mutations.append(parts[1])
        else:
            print(f"Formato inesperado en la mutación: '{m}'")  # Mensaje de error o registro
    return valid_mutations

# Extraer y acumular todas las mutaciones para ambos grupos
all_mutations_mayor_2021 = []
for index, row in df_mayor_2021.iterrows():
    all_mutations_mayor_2021.extend(extract_mutations_aa(row.get('mutaciones', '')))

all_mutations_menor_2021 = []
for index, row in df_menor_2021.iterrows():
    all_mutations_menor_2021.extend(extract_mutations_aa(row.get('mutaciones', '')))

# Eliminar duplicados
unique_mutations_mayor_2021 = set(all_mutations_mayor_2021)
unique_mutations_menor_2021 = set(all_mutations_menor_2021)

# Determinar las mutaciones únicas y remanentes
mutaciones_unicas = unique_mutations_mayor_2021 - unique_mutations_menor_2021
mutaciones_remanentes = unique_mutations_mayor_2021 & unique_mutations_menor_2021

# Visualizar el genoma
fig, ax = plt.subplots(figsize=(13, 3), dpi=300)

# Dibujar las regiones del genoma en aminoácidos
for region in REGIONES:
    ax.plot([region["start"], region["end"]], [0, 0], label=region["product"], linewidth=15)

# Marcar mutaciones únicas (rojo arriba)
mutation_y_pos_unicas = 0.08  # Elevar ligeramente los puntos
for mutation in mutaciones_unicas:
    pos = int(mutation[1:-1])  # Extraer la posición de la mutación
    ax.plot(pos, mutation_y_pos_unicas, 'ro', markersize=2, zorder=5)

# Marcar mutaciones remanentes (azul abajo)
mutation_y_pos_remanentes = -0.08  # Bajar ligeramente los puntos
for mutation in mutaciones_remanentes:
    pos = int(mutation[1:-1])  # Extraer la posición de la mutación
    ax.plot(pos, mutation_y_pos_remanentes, 'bo', markersize=2, zorder=5)

# Añadir líneas verticales grises, opacas, entrecortadas
for pos in [312, 364, 433, 496]:
    ax.axvline(x=pos, color='gray', linestyle='--', alpha=0.5)

# Ajustar los límites del eje y para evitar que los puntos vayan demasiado arriba o abajo
ax.set_ylim(-0.5, 0.5)

ax.set_yticks([])
ax.set_xlabel('Posición de aminoácidos', fontsize=14)
ax.set_xlim(0, int(10264 / 3))  # Dividir el máximo de nucleótidos por 3
ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1.05))
plt.title('Mapa del Genoma del Dengue 3 (Post-pandemia) con Mutaciones Únicas (Rojo) y Remanentes (Azul)', fontsize=16)
plt.show()



#


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 23:11:35 2024

@author: ics2
"""

import matplotlib.pyplot as plt
import pandas as pd

# Cargar el archivo CSV de mutaciones
df = pd.read_csv('/home/ics2/mutaciones_resultado.csv')

# Convertir las regiones a posiciones de aminoácidos
REGIONES = [
    {"start": int((r["start"] - 1) / 3) + 1, "end": int(r["end"] / 3), "gene": r["gene"], "product": r["product"]}
    for r in [
        {"start": 95, "end": 436, "gene": "poly_", "product": "ancC:"},
        {"start": 95, "end": 394, "gene": "poly_", "product": "C:"},
        {"start": 437, "end": 934, "gene": "poly_", "product": "prM:"},
        {"start": 437, "end": 709, "gene": "poly_", "product": "pr:"},
        {"start": 710, "end": 934, "gene": "poly_", "product": "M:"},
        {"start": 935, "end": 2413, "gene": "poly_", "product": "E:"},
        {"start": 2414, "end": 3469, "gene": "poly_", "product": "NS1:"},
        {"start": 3470, "end": 4123, "gene": "poly_", "product": "NS2A:"},
        {"start": 4124, "end": 4513, "gene": "poly_", "product": "NS2B:"},
        {"start": 4514, "end": 6370, "gene": "poly_", "product": "NS3:"},
        {"start": 6371, "end": 6751, "gene": "poly_", "product": "NS4A:"},
        {"start": 6752, "end": 6820, "gene": "poly_", "product": "2K:"},
        {"start": 6821, "end": 7564, "gene": "poly_", "product": "NS4B:"},
        {"start": 7565, "end": 10264, "gene": "poly_", "product": "NS5:"},
    ]
]

# Extraer año
df['año'] = df['nombre'].str.extract(r'/(\d{4})$')

# Convertir la columna 'año' a enteros para facilitar la comparación
df['año'] = df['año'].astype(int)

# Dividir las mutaciones en dos grupos según el año
df_mayor_2021 = df[df['año'] > 2021]
df_menor_2021 = df[df['año'] <= 2021]

def extract_mutations_aa(mutations_str):
    mutations = mutations_str.split(';')
    valid_mutations = []
    for m in mutations:
        parts = m.split(':')
        if len(parts) == 2:  # Verifica que haya exactamente dos partes
            valid_mutations.append(parts[1])
        else:
            print(f"Formato inesperado en la mutación: '{m}'")  # Mensaje de error o registro
    return valid_mutations

# Extraer y acumular todas las mutaciones para ambos grupos
all_mutations_mayor_2021 = []
for index, row in df_mayor_2021.iterrows():
    all_mutations_mayor_2021.extend(extract_mutations_aa(row.get('mutaciones', '')))

all_mutations_menor_2021 = []
for index, row in df_menor_2021.iterrows():
    all_mutations_menor_2021.extend(extract_mutations_aa(row.get('mutaciones', '')))

# Eliminar duplicados
unique_mutations_mayor_2021 = set(all_mutations_mayor_2021)
unique_mutations_menor_2021 = set(all_mutations_menor_2021)

# Determinar las mutaciones únicas y remanentes
mutaciones_unicas = unique_mutations_mayor_2021 - unique_mutations_menor_2021
mutaciones_remanentes = unique_mutations_mayor_2021 & unique_mutations_menor_2021

# Visualizar el genoma
fig, ax = plt.subplots(figsize=(13, 3), dpi=300)

# Dibujar las regiones del genoma en aminoácidos
for region in REGIONES:
    ax.plot([region["start"], region["end"]], [0, 0], label=region["product"], linewidth=15)

# Marcar mutaciones únicas (rojo arriba)
mutation_y_pos_unicas = 0.08  # Elevar ligeramente los puntos
for mutation in mutaciones_unicas:
    pos = int(mutation[1:-1])  # Extraer la posición de la mutación
    ax.plot(pos, mutation_y_pos_unicas, 'ro', markersize=2, zorder=5)

# Marcar mutaciones remanentes (azul abajo)
mutation_y_pos_remanentes = -0.08  # Bajar ligeramente los puntos
for mutation in mutaciones_remanentes:
    pos = int(mutation[1:-1])  # Extraer la posición de la mutación
    ax.plot(pos, mutation_y_pos_remanentes, 'bo', markersize=2, zorder=5)

# Marcar mutaciones pre-pandemia (lila más abajo)
mutation_y_pos_pre_pandemia = -0.15  # Posicionar más abajo
for mutation in unique_mutations_menor_2021:
    pos = int(mutation[1:-1])  # Extraer la posición de la mutación
    ax.plot(pos, mutation_y_pos_pre_pandemia, 'mo', markersize=2, zorder=5) # 'mo' es para color lila/morado

# Añadir líneas verticales grises, opacas, entrecortadas
for pos in [312, 364, 433, 496]:
    ax.axvline(x=pos, color='gray', linestyle='--', alpha=0.5)

# Ajustar los límites del eje y para evitar que los puntos vayan demasiado arriba o abajo
ax.set_ylim(-0.5, 0.5)

ax.set_yticks([])
ax.set_xlabel('Posición de aminoácidos', fontsize=14)
ax.set_xlim(0, int(10264 / 3))  # Dividir el máximo de nucleótidos por 3
ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1.05))
plt.title('Mapa del Genoma del Dengue 3 (Post-pandemia) con Mutaciones: Únicas (Rojo), Remanentes (Azul), Pre-pandemia (Lila)', fontsize=16)
plt.show()