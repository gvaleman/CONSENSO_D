import argparse
import sys

def parsear_argumentos():
    parser = argparse.ArgumentParser(description='CONSENSO_D: Ensamblaje y análisis de genomas virales')
    
    # Argumentos globales primero
    parser.add_argument('--gui', action='store_true', help='Interfaz gráfica de Usuario')
    parser.add_argument('--menu', action='store_true', help='Menú para ingresar datos')
    parser.add_argument('--version', action='version', version='CONSENSO_D v1.0')
    parser.add_argument('--help-all', action='store_true', help='Mostrar ayuda detallada de todos los comandos')
    
    # Subcomandos (opcional si se usa --gui o --menu)
    subparsers = parser.add_subparsers(dest='command', required=False, help='Funcionalidades disponibles')

    # Subcomando ASSEMBLY (tu versión original)
    parser_assembly = subparsers.add_parser('ASSEMBLY', help='Ensamblaje de genomas')
    parser_assembly.add_argument('--tec_seq', choices=['NANO', 'ILLUMINA'], required=True, help='Tipo de secuenciación: Nanopore o Illumina')
    parser_assembly.add_argument('--seq_ref', help='Secuencia de referencia')
    parser_assembly.add_argument('--primers', help='Primers usados durante la amplificación (PCR)')
    parser_assembly.add_argument('--data', help='Metadatos: Nombre de la muestra vs barcode')
    parser_assembly.add_argument('virus', choices=['DENV_1', 'DENV_2', 'DENV_3', 'DENV_4', 'SARS_COV_2', 'RABV', 'N_RABV'], help='Virus a ensamblar')
    parser_assembly.add_argument('path', help='Ruta de los datos de secuenciación')

    # Subcomando VARIANT (nuevo, simple)
    parser_variant = subparsers.add_parser('VARIANT', help='Llamado de variantes')
    parser_variant.add_argument('serotipo', choices=['DENV_1', 'DENV_2', 'DENV_3', 'DENV_4', 'RABV'], help='Serotipo/virus')
    parser_variant.add_argument('archivo_fasta', help='Archivo multiFASTA con secuencias')

    # Subcomando LINEAGE (existente)
    parser_lineage = subparsers.add_parser('LINEAGE', help='Determinación de linajes')
    parser_lineage.add_argument('archivo_fasta', help='Archivo FASTA con las secuencias')
    parser_lineage.add_argument('directorio_salida', help='Directorio donde se guardarán los resultados')
    parser_lineage.add_argument('serotipo', choices=['DENV_1', 'DENV_2', 'DENV_3', 'DENV_4'], help='Serotipo del dengue')

    # Argumentos globales al final (para compatibilidad)
    # Ya están definidos arriba
    
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
