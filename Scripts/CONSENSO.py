import argparse
import sys
from args_parser import parsear_argumentos
from menu_consenso import main as mostrar_menu
# from assemble import ensamblar
import subprocess

def parsear_argumentos():
    parser = argparse.ArgumentParser(description='Ensamblaje y genotipificación de genomas virales')
    
    # Argumentos del ensamblaje
    parser.add_argument('--tec_seq', choices=['NANO', 'ILLUMINA'], help='Tipo de secuenciación: Nanopore o Illumina', required=True)
    parser.add_argument('--seq_ref', help='Secuencia de referencia')
    parser.add_argument('--primers', help='Primers usados durante la amplificación (PCR)')
    parser.add_argument('--data', help='Metadatos: Nombre de la muestra vs barcode')
    parser.add_argument('--gui', action='store_true', help='Interfaz gráfica de Usuario')
    parser.add_argument('--menu', action='store_true', help='Menú para ingresar datos')
    parser.add_argument('virus', help='Virus a ensamblar')
    parser.add_argument('path', help='Ruta de los datos de secuenciación')

    # Argumentos para el script de genotipificación
    parser.add_argument('archivo_fasta', help='Archivo FASTA con las secuencias')
    parser.add_argument('directorio_salida', help='Directorio donde se guardarán los resultados')
    parser.add_argument('serotipo', help='Serotipo del dengue (DENV_1, DENV_2, DENV_3, DENV_4)')

    args = parser.parse_args()

    # Validación de serotipo
    if not args.serotipo.startswith('DENV_') or not args.serotipo.split('_')[1].isdigit() or int(args.serotipo.split('_')[1]) not in range(1, 5):
        parser.error("El serotipo debe ser uno de: DENV_1, DENV_2, DENV_3, DENV_4")

    return args

def ejecutar_script_bash(fasta_file, output_dir, serotype):
    result = subprocess.run(
        ["bash", "tu_script_bash.sh", fasta_file, output_dir, serotype],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print(f"Error ejecutando el script bash: {result.stderr}")
    else:
        print(f"Resultado del script bash: {result.stdout}")

def main():
    args = parsear_argumentos()

    if len(sys.argv) == 1:
        # No se proporcionaron argumentos, mostrar el menú
        mostrar_menu()
    else:
        if args.gui:
            import gui
        elif args.menu:
            args = manejar_menu(args)
        
        # Verificar si hay información suficiente para ejecutar el ensamblaje
        if args.virus and args.path:
            ensamblar(args.virus, args.path, args.tec_seq, args.seq_ref, args.primers, args.data)
            
            # Ejecutar el script bash para genotipificación
            ejecutar_script_bash(args.archivo_fasta, args.directorio_salida, args.serotipo)
        else:
            print("Falta información para ejecutar el ensamblaje. Usa --menu o --gui para ingresar los datos.")

if __name__ == "__main__":
    main()