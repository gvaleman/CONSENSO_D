#!/usr/bin/env python3
import os
import re
import shutil
import argparse

def organizar_fastq(directorio=None):
    """
    Organiza archivos FASTQ en carpetas según el nombre de la muestra.
    """
    # Usar directorio actual si no se especifica
    if directorio is None:
        directorio = os.getcwd()
    
    # Cambiar al directorio especificado
    os.chdir(directorio)
    
    # Patrón para identificar archivos R1
    patron = re.compile(r'(.+?)_S\d+_L001_R1_001\.fastq\.gz')
    
    # Contador para estadísticas
    muestras_procesadas = 0
    archivos_movidos = 0
    
    print(f"Organizando archivos FASTQ en: {directorio}")
    
    # Recorrer todos los archivos en el directorio
    for archivo in os.listdir('.'):
        # Buscar archivos R1
        coincidencia = patron.match(archivo)
        if coincidencia:
            nombre_muestra = coincidencia.group(1)
            archivo_r1 = archivo
            archivo_r2 = archivo.replace('_R1_', '_R2_')
            
            # Verificar si existe el archivo R2 correspondiente
            if os.path.isfile(archivo_r2):
                print(f"Procesando muestra: {nombre_muestra}")
                
                # Crear directorio si no existe
                if not os.path.exists(nombre_muestra):
                    os.makedirs(nombre_muestra)
                
                # Mover archivos R1 y R2 al directorio
                shutil.move(archivo_r1, os.path.join(nombre_muestra, archivo_r1))
                shutil.move(archivo_r2, os.path.join(nombre_muestra, archivo_r2))
                
                print(f"  → Archivos movidos a {nombre_muestra}/")
                muestras_procesadas += 1
                archivos_movidos += 2
            else:
                print(f"ADVERTENCIA: No se encontró archivo R2 para {archivo_r1}")
    
    print(f"Organización completada. {muestras_procesadas} muestras procesadas, {archivos_movidos} archivos movidos.")
    return muestras_procesadas, archivos_movidos

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Organiza archivos FASTQ en carpetas por muestra")
    parser.add_argument("--dir", help="Directorio donde se encuentran los archivos FASTQ", default=None)
    args = parser.parse_args()
    
    organizar_fastq(args.dir)
