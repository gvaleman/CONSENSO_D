import sys
from args_parser import parsear_argumentos
from menu_consenso import main as mostrar_menu
#from assemble import ensamblar

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
        
        # VerificaR si hay información suficiente para ejecutar el ensamblaje
        if args.virus and args.path:
            ensamblar(args.virus, args.path, args.tec_seq, args.seq_ref, args.primers, args.data)
        else:
            print("Falta información para ejecutar el ensamblaje. Usa --menu o --gui para ingresar los datos.")

if __name__ == "__main__":
    main()