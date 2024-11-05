import argparse
import subprocess
import sys
from pathlib import Path

def script_path(script_name):
    # Obtener el directorio del script actual y navegar a Scripts
    base_dir = Path(__file__).resolve().parent.parent
    return base_dir / 'Scripts' / script_name

def parsear_argumentos():
    parser = argparse.ArgumentParser(description='Ensamblaje y genotipificación de genomas virales')
    subparsers = parser.add_subparsers(dest='command', required=True, help='Sub-commands: ASSEMBLY, LINEAGE, MUTATION')

    # Subcomando ASSEMBLY
    parser_assembly = subparsers.add_parser('ASSEMBLY', help='Ejecución del ensamblaje')
    parser_assembly.add_argument('--tec_seq', choices=['NANO', 'ILLUMINA'], required=True, help='Tipo de secuenciación: Nanopore o Illumina')
    parser_assembly.add_argument('serotipo', choices=['DENV_1', 'DENV_2', 'DENV_3', 'DENV_4', 'SARS_COV_2', 'RABV'], help='Virus o serotipo a ensamblar')
    parser_assembly.add_argument('path', help='Ruta de los datos de secuenciación')

    # Subcomando LINEAGE
    parser_lineage = subparsers.add_parser('LINEAGE', help='Ejecución de genotipificación')
    parser_lineage.add_argument('archivo_fasta', help='Archivo FASTA con las secuencias')
    parser_lineage.add_argument('directorio_salida', help='Directorio donde se guardarán los resultados')
    parser_lineage.add_argument('serotipo', choices=['DENV_1', 'DENV_2', 'DENV_3', 'DENV_4'], help='Serotipo del dengue')

    # Subcomando MUTATION
    parser_mutation = subparsers.add_parser('MUTATION', help='Análisis de mutaciones')
    parser_mutation.add_argument('serotipo', choices=['DENV_1', 'DENV_2', 'DENV_3', 'DENV_4', 'RABV'], help='Serotipo o Virus')
    parser_mutation.add_argument('archivo_fasta', help='Archivo FASTA con las secuencias')

    return parser.parse_args()

def ejecutar_script(script_name, *args):
    script_file = script_path(script_name)
    subprocess.run(['bash', str(script_file)] + list(args), check=True)

def main():
    if len(sys.argv) == 1:
        # No se proporcionaron argumentos, mostrar el menú
        menu_script = script_path('menu_consenso.py')
        subprocess.run(['python3', str(menu_script)])
        return

    args = parsear_argumentos()

    try:
        if args.command == 'ASSEMBLY':
            ejecutar_script('ASSEMBLER.sh', args.tec_seq, args.serotipo, args.path)
        elif args.command == 'LINEAGE':
            ejecutar_script('GLUE_lineage.sh', args.archivo_fasta, args.directorio_salida, args.serotipo)
        elif args.command == 'MUTATION':
            ejecutar_script('Variant_calling.sh', args.serotipo, args.archivo_fasta)
        else:
            print("Comando no reconocido. Use ASSEMBLY, LINEAGE o MUTATION.")
    except subprocess.CalledProcessError as e:
        print(f"Error en la ejecución del script {args.command}: {e}")

if __name__ == "__main__":
    main()



"""
import argparse
import sys
import subprocess
import os

def parsear_argumentos():
    parser = argparse.ArgumentParser(description='Ensamblaje y genotipificación de genomas virales')
    parser.add_argument('--tec_seq', choices=['NANO', 'ILLUMINA'], help='Tipo de secuenciación: Nanopore o Illumina')
    parser.add_argument('--seq_ref', help='Secuencia de referencia')
    parser.add_argument('--primers', help='Primers usados durante la amplificación (PCR)')
    parser.add_argument('--data', help='Metadatos: Nombre de la muestra vs barcode')
    parser.add_argument('--gui', action='store_true', help='Interfaz gráfica de Usuario')
    parser.add_argument('--menu', action='store_true', help='Menú para ingresar datos')
    parser.add_argument('virus', nargs='?', help='Virus a ensamblar')
    parser.add_argument('path', nargs='?', help='Ruta de los datos de secuenciación')
    parser.add_argument('archivo_fasta', nargs='?', help='Archivo FASTA con las secuencias')
    parser.add_argument('directorio_salida', nargs='?', help='Directorio donde se guardarán los resultados')
    parser.add_argument('serotipo', nargs='?', help='Serotipo del dengue (DENV_1, DENV_2, DENV_3, DENV_4)')

    args = parser.parse_args()

    # Validación de serotipo si se proporciona
    if args.serotipo and (not args.serotipo.startswith('DENV_') or not args.serotipo.split('_')[1].isdigit() or int(args.serotipo.split('_')[1]) not in range(1, 5)):
        parser.error("El serotipo debe ser uno de: DENV_1, DENV_2, DENV_3, DENV_4")

    return args

def main():
    if len(sys.argv) == 1:
        # No se proporcionaron argumentos, mostrar el menú
        subprocess.run(["python3", "menu_consenso.py"])
        return

    if len(sys.argv) == 2:
        arg = sys.argv[1].lower()
        if arg == "menu":
            subprocess.run(["python3", "menu_consenso.py"])
            return
        elif arg == "gui":
            subprocess.run(["python3", "GUI CONSENSO_D.py"])
            return

    args = parsear_argumentos()

    if args.gui:
        subprocess.run(["python3", "GUI CONSENSO_D.py"])
    elif args.menu:
        subprocess.run(["python3", "menu_consenso.py"])
    
    # Verificar si hay información suficiente para ejecutar el ensamblaje
    if args.virus and args.path:
        # Asumiendo que la función ensamblar está definida en algún lugar.
        ensamblar(args.virus, args.path, args.tec_seq, args.seq_ref, args.primers, args.data)
        
        # Ejecutar el script bash para genotipificación
        if args.archivo_fasta and args.directorio_salida and args.serotipo:
            ejecutar_script_bash(args.archivo_fasta, args.directorio_salida, args.serotipo)
        else:
            print("Falta información para ejecutar el script de genotipificación.")
    else:
        print("Falta información para ejecutar el ensamblaje. Usa --menu o --gui para ingresar los datos.")

if __name__ == "__main__":
    main()
    """