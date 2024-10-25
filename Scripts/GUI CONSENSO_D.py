import tkinter as tk
from tkinter import filedialog, messagebox, Toplevel, Label, Text, Scrollbar
import subprocess
from threading import Thread
from queue import Queue, Empty
from PIL import Image, ImageTk

def ejecutar_consenso_d():
    virus = virus_var.get()
    ruta = ruta_var.get()

    if not virus or not ruta:
        messagebox.showerror("Error", "Seleccione un tipo de virus y/o ruta.")
        return

    # Crear ventana de progreso
    ventana_progreso = Toplevel(root)
    ventana_progreso.title("Procesando...")
    Label(ventana_progreso, text="Procesando...").pack(pady=10)

    # Crear un área de texto con barra de desplazamiento
    area_texto = Text(ventana_progreso, wrap='word', height=20, width=80)
    area_texto.pack(side='left', fill='both', expand=True)

    barra_desplazamiento = Scrollbar(ventana_progreso, command=area_texto.yview)
    barra_desplazamiento.pack(side='right', fill='y')
    area_texto['yscrollcommand'] = barra_desplazamiento.set

    def enqueue_output(out, queue):
        for line in iter(out.readline, b''):
            queue.put(line)
        out.close()


#Funcion Externa
    def ejecutar_script_y_actualizar_salida():
        try:
            # Ejecutar el script CONSENSO_D con los parámetros seleccionados
            proceso = subprocess.Popen(['/home/ics2/CONSENSO_D/Scripts/CONSENSO', virus, ruta],
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

    # Ejecutar el script en un hilo separado para no bloquear la interfaz gráfica
    hilo = Thread(target=ejecutar_script_y_actualizar_salida)
    hilo.start()

def buscar_ruta():
    directorio = filedialog.askdirectory()
    ruta_var.set(directorio)

def salir_app():
    root.quit()

def acerca_de_consenso_d():
    messagebox.showinfo("Acerca de CONSENSO_D", "CONSENSO_D Viral Genome Assembler\nVersión 0.2")

# Crear ventana principal
ventana = tk.Tk()
ventana.title("CONSENSO_D Viral Genome Assembler")

# Crear barra de menú 
menubar = tk.Menu(ventana)
ventana.config(menu=menubar)

# Crear menú File
menu_archivo = tk.Menu(menubar, tearoff=0)
menu_archivo.add_command(label="Salir", command=salir_app)
menubar.add_cascade(label="Archivo", menu=menu_archivo)

# Crear menú Help
menu_ayuda = tk.Menu(menubar, tearoff=0)
menu_ayuda.add_command(label="Acerca de CONSENSO_D", command=acerca_de_consenso_d)
menubar.add_cascade(label="Ayuda", menu=menu_ayuda)

# Subir imagen
ruta_imagen = '/home/ics2/CONSENSO_D/Resources/consenso logo.png'
imagen = Image.open(ruta_imagen)
mostrar = ImageTk.PhotoImage(imagen)

# Crear y configurar el marco principal
marco_principal = tk.Frame(ventana, padx=10, pady=10)
marco_principal.pack(padx=10, pady=10)

# Añadir la imagen a la interfaz
etiqueta_imagen = tk.Label(marco_principal, image=mostrar)
etiqueta_imagen.image = mostrar  # Guardar una referencia para evitar recolección de basura
etiqueta_imagen.grid(row=0, column=0, columnspan=3, pady=10)

# Etiqueta y campo desplegable para seleccionar el tipo de virus
tk.Label(marco_principal, text="Seleccionar Tipo de Virus:").grid(row=1, column=0, sticky="w")
virus_var = tk.StringVar()
opciones_virus = ["DENV_1", "DENV_2", "DENV_3", "DENV_4", "RABV"]
tk.OptionMenu(marco_principal, virus_var, *opciones_virus).grid(row=1, column=1, sticky="ew")

# Etiqueta y campo para seleccionar el directorio
tk.Label(marco_principal, text="Seleccionar Ruta del Directorio:").grid(row=2, column=0, sticky="w")
ruta_var = tk.StringVar()
tk.Entry(marco_principal, textvariable=ruta_var, width=50).grid(row=2, column=1, sticky="ew")
tk.Button(marco_principal, text="Buscar", command=buscar_ruta).grid(row=2, column=2, sticky="ew")

# Botón para ejecutar CONSENSO_D
tk.Button(marco_principal, text="Ejecutar CONSENSO_D", command=ejecutar_consenso_d).grid(row=3, column=0, columnspan=3, pady=10)

# Ejecutar el bucle principal de la interfaz gráfica
ventana.mainloop()