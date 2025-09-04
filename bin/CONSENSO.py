import argparse
import subprocess
import sys
from pathlib import Path

def script_path(script_name):
    # Obtener el directorio del script actual y navegar a Scripts
    base_dir = Path(__file__).resolve().parent.parent
    return base_dir / 'Scripts' / script_name

def parsear_argumentos():
    parser = argparse.ArgumentParser(description='CONSENSO_D: Ensamblaje y análisis de genomas virales')
    
    # Argumentos globales primero
    parser.add_argument('--gui', action='store_true', help='Interfaz gráfica de Usuario')
    parser.add_argument('--menu', action='store_true', help='Menú para ingresar datos')
    parser.add_argument('--version', action='version', version='CONSENSO_D v1.0')
    parser.add_argument('--help-all', action='store_true', help='Mostrar ayuda detallada de todos los comandos')
    parser.add_argument('--log', help='Guardar la salida en un archivo de log (por defecto: comando_timestamp.log)')
    
    # Subcomandos (opcional si se usa --gui o --menu)
    subparsers = parser.add_subparsers(dest='command', required=False, help='Funcionalidades disponibles')
    
    # Subcomando ASSEMBLY
    parser_assembly = subparsers.add_parser('ASSEMBLY', help='Ensamblaje de genomas')
    parser_assembly.add_argument('--tec_seq', choices=['NANO', 'ILLUMINA'], required=True, help='Tipo de secuenciación: Nanopore o Illumina')
    parser_assembly.add_argument('--seq_ref', help='Secuencia de referencia')
    parser_assembly.add_argument('--primers', help='Primers usados durante la amplificación (PCR)')
    parser_assembly.add_argument('--data', help='Metadatos: Nombre de la muestra vs barcode')
    parser_assembly.add_argument('virus', choices=['DENV_1', 'DENV_2', 'DENV_3', 'DENV_4', 'SARS_COV_2', 'RABV', 'N_RABV'], help='Virus a ensamblar')
    parser_assembly.add_argument('path', help='Ruta de los datos de secuenciación')
    
    # Subcomando LINEAGE
    parser_lineage = subparsers.add_parser('LINEAGE', help='Determinación de linajes')
    parser_lineage.add_argument('archivo_fasta', help='Archivo FASTA con las secuencias')
    parser_lineage.add_argument('directorio_salida', help='Directorio donde se guardarán los resultados')
    parser_lineage.add_argument('serotipo', choices=['DENV_1', 'DENV_2', 'DENV_3', 'DENV_4'], help='Serotipo del dengue')
    
    # Subcomando VARIANT (nuevo)
    parser_variant = subparsers.add_parser('VARIANT', help='Llamado de variantes')
    parser_variant.add_argument('serotipo', choices=['DENV_1', 'DENV_2', 'DENV_3', 'DENV_4', 'RABV'], help='Serotipo/virus')
    parser_variant.add_argument('archivo_fasta', help='Archivo multiFASTA con secuencias')
    
    # Subcomando MUTATION (mantener compatibilidad)
    parser_mutation = subparsers.add_parser('MUTATION', help='Análisis de mutaciones')
    parser_mutation.add_argument('serotipo', choices=['DENV_1', 'DENV_2', 'DENV_3', 'DENV_4', 'RABV'], help='Serotipo o Virus')
    parser_mutation.add_argument('archivo_fasta', help='Archivo FASTA con las secuencias')
    
    # Nuevo subcomando ORGANIZE_FASTQ
    parser_organize = subparsers.add_parser('ORGANIZE_FASTQ', help='Organizar archivos FASTQ en carpetas')
    parser_organize.add_argument('directorio', help='Directorio con archivos FASTQ para organizar')
    
    # Manejar casos especiales
    if len(sys.argv) == 1:
        # Sin argumentos, mostrar menú por defecto
        return type('Args', (), {'gui': False, 'menu': True, 'command': None})()
    
    args = parser.parse_args()
    
    # Si se pide GUI o MENU, no requiere subcomando
    if args.gui or args.menu:
        return args
    
    # Si no hay subcomando y no es GUI/MENU, mostrar ayuda
    if not args.command:
        parser.print_help()
        sys.exit(1)
        
    return args

def ejecutar_script(script_name, *args, log_file=None):
    script_file = script_path(script_name)
    
    if log_file:
        # Redirigir salida a un archivo de log
        with open(log_file, 'w') as f:
            result = subprocess.run(
                ['bash', str(script_file)] + list(args),
                stdout=f,
                stderr=subprocess.STDOUT,  # Combina stderr con stdout
                check=False  # No lanzar excepción para poder capturar el código de retorno
            )
        return result.returncode
    else:
        # Comportamiento normal
        return subprocess.run(['bash', str(script_file)] + list(args), check=True)

def ejecutar_script_python(script_name, *args, log_file=None):
    script_file = script_path(script_name)
    
    if log_file:
        # Redirigir salida a un archivo de log
        with open(log_file, 'w') as f:
            result = subprocess.run(
                ['python3', str(script_file)] + list(args),
                stdout=f,
                stderr=subprocess.STDOUT,
                check=False
            )
        return result.returncode
    else:
        # Comportamiento normal
        return subprocess.run(['python3', str(script_file)] + list(args), check=True)

def main():
    args = parsear_argumentos()
    
    try:
        # Determinar archivo de log si se especificó
        log_file = None
        if hasattr(args, 'log') and args.log:
            if args.log == 'auto':
                # Generar nombre automático basado en comando y timestamp
                from datetime import datetime
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                log_file = f"{args.command if hasattr(args, 'command') else 'consenso'}_{timestamp}.log"
            else:
                log_file = args.log
        
        # Manejar GUI
        if hasattr(args, 'gui') and args.gui:
            gui_script = script_path('GUI CONSENSO_D.py')
            subprocess.run(['python3', str(gui_script)])
            return
        
        # Manejar MENU
        if hasattr(args, 'menu') and args.menu:
            menu_script = script_path('menu_consenso.py')
            subprocess.run(['python3', str(menu_script)])
            return
        
        # Manejar subcomandos
        if args.command == 'ASSEMBLY':
            # Usar virus en lugar de serotipo para mantener compatibilidad
            return_code = ejecutar_script('ASSEMBLER.sh', args.tec_seq, args.virus, args.path, log_file=log_file)
            if log_file:
                if return_code == 0:
                    print(f"Ensamblaje completado exitosamente. Log guardado en: {log_file}")
                else:
                    print(f"Error en el ensamblaje. Revise el log en: {log_file}")
        elif args.command == 'LINEAGE':
            return_code = ejecutar_script('GLUE_lineage.sh', args.archivo_fasta, args.directorio_salida, args.serotipo, log_file=log_file)
            if log_file:
                print(f"Análisis de linaje completado. Log guardado en: {log_file}")
        elif args.command == 'VARIANT':
            return_code = ejecutar_script('Variant_caller.sh', args.serotipo, args.archivo_fasta, log_file=log_file)
            if log_file:
                print(f"Análisis de variantes completado. Log guardado en: {log_file}")
        elif args.command == 'MUTATION':
            return_code = ejecutar_script('Variant_calling.sh', args.serotipo, args.archivo_fasta, log_file=log_file)
            if log_file:
                print(f"Análisis de mutaciones completado. Log guardado en: {log_file}")
        elif args.command == 'ORGANIZE_FASTQ':
            # Usar el nuevo script Python para organizar FASTQ
            return_code = ejecutar_script_python('organizar_fastq.py', '--dir', args.directorio, log_file=log_file)
            if log_file:
                print(f"Organización de archivos FASTQ completada. Log guardado en: {log_file}")
        else:
            print("Comando no reconocido.")
            
    except subprocess.CalledProcessError as e:
        print(f"Error en la ejecución del script {args.command}: {e}")

if __name__ == "__main__":
    main()
