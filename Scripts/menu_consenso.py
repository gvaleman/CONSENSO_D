import time
import os

def welcome_message():
    print("""
=============================================
=         Bienvenido a CONSENSO_D           =
=      Ensamblador de genomas virales       =
=                                           =
=     developed by / desarollado por:       =
=   Sustainable Sciences Institute Team     =
=============================================""")
    time.sleep(2)
    print()  # Adding a new line for better readability

def main_menu():
    print("""
=================================================
=                 MENÚPRINCIPAL                =
=                   MAIN MENU                   =
=================================================
   |  1. Usar CONSENSO_D en consola           |
   |      (Use CONSENSO_D in console)         |
   |                                          |
   |  2. Usar CONSENSO_D de manera guiada     |
   |      (Use CONSENSO_D in guided console)  |
   |                                          |
   |  3. Usar CONSENSO_D en GUI               |
   |       (Use CONSENSO_D in GUI)            |
   |                                          |
   |  4. Determinar Linaje con GLU CLI        |
   |                                          |
   |  5. Salir  / Exit                        |
   |                                          |
=================================================""")
    
def use_console():
    print("""
***********************************************
*          Usar CONSENSO_D en consola         *
***********************************************

ESPAÑOL:
A continuación puedes usar la consola.
    Parámetros:
        CONSENSO_D: Nombre del script (es una constante)
        VIRUS: Nombre del virus (DENV_1, DENV_2, DENV_3, DENV_4, RABV)
        path: Ubicación del directorio donde están los archivos fastq crudos

    Ejemplo:
        CONSENSO_D DENV_1 "~/Documentos/fast_files"
          

ENGLISH:
You can now use the console.
    Parameters:
        CONSENSO_D: Name of the script (constant)
        VIRUS: Name of the virus (DENV_1, DENV_2, DENV_3, DENV_4, RABV)
        path: Directory location where the raw fastq files are located

    Example:
        CONSENSO_D DENV_1 "~/Documents/fast_files"
          
""")
    exit_program() 


def use_gui():
    print("***********************************************")
    print("*          Usando CONSENSO_D en GUI           *")
    print("***********************************************")
    print("Lanzando interfaz gráfica...")

    # Obtener la ruta del directorio donde está el script actual
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Construir la ruta relativa al archivo de la GUI
    gui_script_path = os.path.join(current_dir, 'GUI CONSENSO_D.py')
    
    # Ejecutar el script de la GUI
    os.system(f'python3 "{gui_script_path}"')


def guided_console():
    viruses = ["DENV_1", "DENV_2", "DENV_3", "DENV_4", "SARS_COV_2", "RABV"]

    print("\n¿Qué virus desea ensamblar?")
    print("(Which virus do you want to assemble?)")
    for i, virus in enumerate(viruses, 1):
        print(f"{i}. {virus}")
    
    virus_choice = input("Elija el número del virus (choose the number of the virus): ")
    try:
        virus_choice = int(virus_choice)
        if virus_choice < 1 or virus_choice > len(viruses):
            raise ValueError

        virus = viruses[virus_choice - 1]
    except ValueError:
        print("*************************************************")
        print("* Elección inválida, por favor intente de nuevo! *")
        print("*       Invalid choice, please try again!        *")
        print("*************************************************")
        return

    path = input("\nArrastre el path donde tiene los archivos fastq (Drag the path where the fastq files are located): ")

    command = f"./CONSENSO {virus} \"{path}\""
    print("\nEjecutando el siguiente comando:")
    print(f"(Running the following command:)\n{command}")

    os.system(command)


def GLUE_lineage():
    viruses = ["DENV_1", "DENV_2", "DENV_3", "DENV_4"]

    print("\n¿Qué serotipo desea analizar?")
    print("(Which serotype do you want to analyze?)")
    print("(DENV_1, DENV_2, DENV_3, DENV_4)")
    
    virus_choice = input("Ingrese el serotipo (enter the serotype): ").upper().strip()
    if virus_choice not in viruses:
        print("*************************************************")
        print("* Serotipo inválido, por favor intente de nuevo! *")
        print("*    Invalid serotype, please try again!        *")
        print("*************************************************")
        return

    path = input("\nArrastre el archivo fasta de sus secuencias: ")
    
    if not os.path.exists(path):
        print("El archivo especificado no existe. (The specified file does not exist.)")
        return

    output_dir = os.path.dirname(path)  # Usa el mismo directorio del archivo de entrada
    command = f"bash '/home/ics2/CONSENSO_D/Scripts/GLUE_lineage.sh' \"{path}\" \"{output_dir}\" \"{virus_choice}\""
    print("\nEjecutando el siguiente comando:")
    print(f"(Running the following command:)\n{command}")

    os.system(command)

def exit_program():
    print("***********************************************")
    print("*          Saliendo de CONSENSO_D...          *")
    print("        ¡Gracias por usar CONSENSO_D          *")
    print("***********************************************")
    exit()

def invalid_choice():
    print("***********************************************")
    print("*      Opción inváliad, Intentar de nuevo      *")
    print("***********************************************")

def main():
    welcome_message()

    while True:
        main_menu()
        choice = input("Ingrese su elección (1-5) / Enter your choice (1-5): ")

        if choice == '1':
            use_console()
        elif choice == '2':
            guided_console()
        elif choice == '3':
            use_gui()
        elif choice == '4':
            GLUE_lineage()
        elif choice == '5':
            exit_program()
        else:
            invalid_choice()

if __name__ == "__main__":
    main()