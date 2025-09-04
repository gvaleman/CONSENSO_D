import argparse

def parsear_argumentos():
    parser = argparse.ArgumentParser(description='Ensamblaje de genomas virales')
    parser.add_argument('--tec_seq', choices=['NANO', 'ILLUMINA'], help='Tipo de secuenciación: Nanopore o Illumina', required=True)
    parser.add_argument('--seq_ref', help='Secuencia de referencia')
    parser.add_argument('--primers', help='Primers usados durante la amplificación (PCR)')
    parser.add_argument('--data', help='Metadatos: Nombre de la muestra vs barcode')
    parser.add_argument('--gui', action='store_true', help='Interfaz gráfica de Usuario')
    parser.add_argument('--menu', action='store_true', help='Menú para ingresar datos')
    parser.add_argument('virus', help='Virus a ensamblar')
    parser.add_argument('path', help='Ruta de los datos de secuenciación')
    
    return parser.parse_args()