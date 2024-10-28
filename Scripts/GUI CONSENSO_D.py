import os
import tkinter as tk
from tkinter import filedialog, messagebox, Toplevel, Text, Scrollbar
from PIL import Image, ImageTk
import subprocess
from threading import Thread
from queue import Queue, Empty

def ejecutar_consenso_d():
    sequence_type_display = sequence_type_var.get()
    # Convertir el valor mostrado al valor interno
    sequence_type = opciones_secuenciacion_map[sequence_type_display]
    virus = virus_var.get()
    ruta = ruta_var.get()

    if not sequence_type or not virus or not ruta:
        messagebox.showerror("Error", "Seleccione un tipo de secuenciación, tipo de virus y/o ruta.")
        return

    ventana_progreso = Toplevel(ventana)
    ventana_progreso.title("Procesando...")
    tk.Label(ventana_progreso, text="Procesando...").pack(pady=10)

    area_texto = Text(ventana_progreso, wrap='word', height=20, width=80)
    area_texto.pack(side='left', fill='both', expand=True)

    barra_desplazamiento = Scrollbar(ventana_progreso, command=area_texto.yview)
    barra_desplazamiento.pack(side='right', fill='y')
    area_texto['yscrollcommand'] = barra_desplazamiento.set

    def enqueue_output(out, queue):
        for line in iter(out.readline, b''):
            queue.put(line)
        out.close()

    def ejecutar_script_y_actualizar_salida():
        try:
            # Obtener el directorio actual del script
            current_dir = os.path.dirname(os.path.abspath(__file__))

            # Construir la ruta relativa para el archivo ASSEMBLER.sh
            script_path = os.path.join(current_dir, 'ASSEMBLER.sh')

            proceso = subprocess.Popen([script_path, sequence_type, virus, ruta],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       text=True)

            queue = Queue()
            Thread(target=enqueue_output, args=(proceso.stdout, queue)).start()
            Thread(target=enqueue_output, args=(proceso.stderr, queue)).start()

            while True:
                try:
                    line = queue.get_nowait()
                except Empty:
                    if proceso.poll() is not None:
                        break
                else:
                    area_texto.insert(tk.END, line)
                    area_texto.see(tk.END)

            proceso.wait()
            if proceso.returncode == 0:
                messagebox.showinfo("Éxito", "¡CONSENSO_D ejecutado con éxito!")
            else:
                messagebox.showerror("Error", "Ocurrió un error durante la ejecución.")
        except Exception as e:
            messagebox.showerror("Error", f"Ocurrió un error: {e}")
        finally:
            ventana_progreso.destroy()

    hilo = Thread(target=ejecutar_script_y_actualizar_salida)
    hilo.start()

def buscar_ruta():
    directorio = filedialog.askdirectory()
    ruta_var.set(directorio)

def salir_app():
    ventana.quit()

def acerca_de_consenso_d():
    messagebox.showinfo("Acerca de CONSENSO_D", "CONSENSO_D Viral Genome Assembler\nVersión 0.2")

# Crear ventana principal
ventana = tk.Tk()
ventana.title("CONSENSO_D: Viral Genome Assembler")

# Configurar el icono de la ventana usando iconphoto
try:
    current_dir = os.path.dirname(os.path.abspath(__file__))
    ruta_icono = os.path.join(current_dir, '..', 'Resources', 'icon.png')
    icono_imagen = Image.open(ruta_icono)
    icono = ImageTk.PhotoImage(icono_imagen)
    ventana.iconphoto(True, icono)
except Exception as e:
    print(f"Error al cargar el icono: {e}")
    messagebox.showerror("Error", f"No se pudo cargar el icono: {e}")

# Crear barra de menú 
menubar = tk.Menu(ventana)
ventana.config(menu=menubar)

menu_archivo = tk.Menu(menubar, tearoff=0)
menu_archivo.add_command(label="Salir", command=salir_app)
menubar.add_cascade(label="Archivo", menu=menu_archivo)

menu_ayuda = tk.Menu(menubar, tearoff=0)
menu_ayuda.add_command(label="Acerca de CONSENSO_D", command=acerca_de_consenso_d)
menubar.add_cascade(label="Ayuda", menu=menu_ayuda)

# Crear y configurar el marco principal
marco_principal = tk.Frame(ventana, padx=10, pady=10)
marco_principal.pack(padx=10, pady=10)

# Añadir la imagen a la interfaz
try:
    ruta_imagen = os.path.join(current_dir, '..', 'Resources', 'consenso logo.png')
    imagen = Image.open(ruta_imagen)
    mostrar = ImageTk.PhotoImage(imagen)
    etiqueta_imagen = tk.Label(marco_principal, image=mostrar)
    etiqueta_imagen.image = mostrar
    etiqueta_imagen.grid(row=0, column=0, columnspan=3, pady=10)
except Exception as e:
    messagebox.showerror("Error", f"No se pudo cargar la imagen: {e}")

# Diccionario para mapear la representación amigable al valor interno
opciones_secuenciacion_map = {"Nanopore": "NANO", "Illumina": "ILLUMINA"}

# Etiqueta y campo desplegable para seleccionar el tipo de secuenciación
tk.Label(marco_principal, text="Seleccionar Tipo de Secuenciación:").grid(row=1, column=0, sticky="w")
sequence_type_var = tk.StringVar()
opciones_secuenciacion_display = list(opciones_secuenciacion_map.keys())
tk.OptionMenu(marco_principal, sequence_type_var, *opciones_secuenciacion_display).grid(row=1, column=1, sticky="ew")

# Etiqueta y campo desplegable para seleccionar el tipo de virus
tk.Label(marco_principal, text="Seleccionar Tipo de Virus:").grid(row=2, column=0, sticky="w")
virus_var = tk.StringVar()
opciones_virus = ["DENV_1", "DENV_2", "DENV_3", "DENV_4", "SARS_COV_2", "RABV"]
tk.OptionMenu(marco_principal, virus_var, *opciones_virus).grid(row=2, column=1, sticky="ew")

# Etiqueta y campo para seleccionar el directorio
tk.Label(marco_principal, text="Seleccionar Ruta del Directorio:").grid(row=3, column=0, sticky="w")
ruta_var = tk.StringVar()
tk.Entry(marco_principal, textvariable=ruta_var, width=50).grid(row=3, column=1, sticky="ew")
tk.Button(marco_principal, text="Buscar", command=buscar_ruta).grid(row=3, column=2, sticky="ew")

# Botón para ejecutar CONSENSO_D con estilo personalizado
boton_ejecutar = tk.Button(
    marco_principal, 
    text="Ejecutar CONSENSO_D", 
    command=ejecutar_consenso_d,
    bg="#1e97f3", 
    fg="white", 
    activebackground="#0056b3", 
    activeforeground="white",
    font=("Helvetica", 12, "bold"),  # Fuente en negrita y tamaño más grande
    padx=8,  # Margen horizontal
    pady=4   # Margen vertical
)
boton_ejecutar.grid(row=4, column=0, columnspan=3, pady=10)

# Ejecutar el bucle principal de la interfaz gráfica
ventana.mainloop()