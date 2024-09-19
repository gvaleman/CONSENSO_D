import argparse

def parsear_argumentos():
    parser = argparse.ArgumentParser(description='Descripción de tu script.')
    parser.add_argument('--tec_seq', help='Descripción de tec_seq')
    parser.add_argument('--seq_ref', help='Descripción de seq_ref')
    parser.add_argument('--primers', help='Descripción de primers')
    parser.add_argument('--data', help='Descripción de data')
    parser.add_argument('--gui', action='store_true', help='Usa GUI para ingresar datos')
    parser.add_argument('--menu', action='store_true', help='Muestra el menú para ingresar datos')
    parser.add_argument('virus', nargs='?', help='Virus a ensamblar')
    parser.add_argument('path', nargs='?', help='Ruta de los datos')
    
    return parser.parse_args()