#!/usr/bin/env python3
import subprocess
import os
import argparse
from pathlib import Path

def run_glue_genotyping(input_fasta, output_dir, serotype):
    """
    Automatiza el proceso de genotipificación usando GLUE.
    
    Args:
        input_fasta (str): Ruta al archivo fasta de entrada
        output_dir (str): Directorio donde se guardarán los resultados
        serotype (int): Serotipo de dengue (1-4)
    """
    # Validar serotipo
    if serotype not in [1, 2, 3, 4]:
        raise ValueError("Serotipo debe ser 1, 2, 3 o 4")
    
    # Crear directorio de salida si no existe
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Preparar el nombre del archivo de salida
    output_file = output_dir / f"resultados_genotipo_denv{serotype}.txt"
    
    # Preparar los comandos GLUE
    glue_commands = [
        "console set cmd-output-file-format tab",
        f"console set next-cmd-output-file {output_file}",
        f"module denv{serotype}MaxLikelihoodGenotyper genotype file -f {input_fasta}",
        "exit"
    ]
    
    # Crear archivo temporal con comandos
    temp_command_file = output_dir / "temp_commands.txt"
    with open(temp_command_file, 'w') as f:
        f.write('\n'.join(glue_commands))
    
    try:
        # Ejecutar GLUE con los comandos
        process = subprocess.Popen(
            ['gluetools.sh'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Leer y ejecutar comandos
        with open(temp_command_file, 'r') as f:
            stdout, stderr = process.communicate(f.read())
        
        # Verificar si hubo errores
        if process.returncode != 0:
            print(f"Error al ejecutar GLUE: {stderr}")
            return False
        
        print(f"Genotipificación completada. Resultados guardados en: {output_file}")
        return True
        
    except Exception as e:
        print(f"Error durante la ejecución: {str(e)}")
        return False
    
    finally:
        # Limpiar archivo temporal
        if temp_command_file.exists():
            temp_command_file.unlink()

def main():
    parser = argparse.ArgumentParser(description='Automatización de genotipificación GLUE para dengue')
    parser.add_argument('-i', '--input', required=True, help='Archivo FASTA de entrada')
    parser.add_argument('-o', '--output', required=True, help='Directorio de salida')
    parser.add_argument('-s', '--serotype', type=int, required=True, choices=[1,2,3,4],
                      help='Serotipo de dengue (1-4)')
    
    args = parser.parse_args()
    
    # Ejecutar el proceso
    run_glue_genotyping(args.input, args.output, args.serotype)

if __name__ == "__main__":
    main()
