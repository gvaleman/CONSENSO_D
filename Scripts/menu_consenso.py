import time
import os

def welcome_message():
    print("""
=============================================
=          Welcome to CONSENSO_D            =
=       your Viral Genome Assembler         =
=                                           =
=    developed by / desarollado por:        =
=   Sustainable Sciences Institute Team     =
=============================================""")
    time.sleep(2)
    print()  # Adding a new line for better readability

def main_menu():
    print("""
=================================================
=                 MENU PRINCIPAL                =
=                   MAIN MENU                   =
===============================¡¡¡¡¡=============
   |  1. Usar CONSENSO_D en consola           |
   |      (Use CONSENSO_D in console)         |
   |                                          |
   |  2. Usar CONSENSO_D de manera guiada     |
   |      (Use CONSENSO_D in guided console)  |
   |                                          |
   |  3. Usar CONSENSO_D en GUI               |
   |       (Use CONSENSO_D in GUI)            |
   |                                          |
   |  4. Salir  / Exit                        |
   |                                          |
=================================================""")

def use_console():
    print("""
***********************************************
*            Using CONSENSO_D in Console      *
***********************************************

ENGLISH:
You can now use the console.
    Parameters:
        CONSENSO_D: Name of the script (constant)
        VIRUS: Name of the virus (DENV_1, DENV_2, DENV_3, DENV_4, RABV)
        path: Directory location where the raw fastq files are located

    Example:
        CONSENSO_D DENV_1 "~Documents/fast_files"

ESPAÑOL:
A continuación puedes usar la consola.
    Parámetros:
        CONSENSO_D: Nombre del script (es una constante)
        VIRUS: Nombre del virus (DENV_1, DENV_2, DENV_3, DENV_4, RABV)
        path: Ubicación del directorio donde están los archivos fastq crudos

    Ejemplo:
        CONSENSO_D DENV_1 "~Documentos/fast_files"
""")
    exit_program()  # Exiting to let the user use the console

def use_gui():
    print("***********************************************")
    print("*            Using CONSENSO_D in GUI           *")
    print("***********************************************")
    print("Launching graphical user interface...")
    os.system('python3 /home/ics2/CONSENSO_D/Scripts/GUI\\ CONSENSO_D.py')

def guided_console():
    viruses = ["DENV_1", "DENV_2", "DENV_3", "DENV_4", "RABV"]

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
        print("***********************************************")
        print("*      Elección inválida, por favor intente de nuevo!      *")
        print("*      Invalid choice, please try again!      *")
        print("***********************************************")
        return

    path = input("\nArrastre el path donde tiene los archivos fastq (Drag the path where the fastq files are located): ")

    command = f"./CONSENSO {virus} \"{path}\""
    print("\nEjecutando el siguiente comando:")
    print(f"(Running the following command:)\n{command}")

    os.system(command)

def exit_program():
    print("***********************************************")
    print("*          Exiting CONSENSO_D...              *")
    print("*     Thank you for using our tool!           *")
    print("***********************************************")
    exit()

def invalid_choice():
    print("***********************************************")
    print("*      Invalid choice, please try again!      *")
    print("***********************************************")

def main():
    welcome_message()

    while True:
        main_menu()
        choice = input("Ingrese su elección (1-4) / Enter your choice (1-4): ")

        if choice == '1':
            use_console()
        elif choice == '2':
            guided_console()
        elif choice == '3':
            use_gui()
        elif choice == '4':
            exit_program()
        else:
            invalid_choice()

if __name__ == "__main__":
    main()