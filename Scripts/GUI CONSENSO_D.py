import tkinter as tk
from tkinter import filedialog, messagebox, Toplevel, Label, Text, Scrollbar
import subprocess
from threading import Thread
from queue import Queue, Empty
from PIL import Image, ImageTk

def run_consenso_d():
    virus = virus_var.get()
    path = path_var.get()

    if not virus or not path:
        messagebox.showerror("Error", "Please select a virus type and a directory path.")
        return

    # Crear ventana de progreso
    progress_window = Toplevel(root)
    progress_window.title("Processing")
    Label(progress_window, text="Processing...").pack(pady=10)

    # Crear un área de texto con barra de desplazamiento
    text_area = Text(progress_window, wrap='word', height=20, width=80)
    text_area.pack(side='left', fill='both', expand=True)

    scrollbar = Scrollbar(progress_window, command=text_area.yview)
    scrollbar.pack(side='right', fill='y')
    text_area['yscrollcommand'] = scrollbar.set

    def enqueue_output(out, queue):
        for line in iter(out.readline, b''):
            queue.put(line)
        out.close()

    def run_script_and_update_output():
        try:
            # Ejecutar el script CONSENSO_D con los parámetros seleccionados
            process = subprocess.Popen(['/home/ics2/CONSENSO_D/Scripts/CONSENSO', virus, path],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       text=True)

            queue = Queue()
            Thread(target=enqueue_output, args=(process.stdout, queue)).start()
            Thread(target=enqueue_output, args=(process.stderr, queue)).start()

            while True:
                try:
                    line = queue.get_nowait()
                except Empty:
                    if process.poll() is not None:
                        break
                else:
                    text_area.insert(tk.END, line)
                    text_area.see(tk.END)

            process.wait()
            if process.returncode == 0:
                messagebox.showinfo("Success", "CONSENSO_D executed successfully!")
            else:
                messagebox.showerror("Error", "An error occurred during execution.")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {e}")
        finally:
            progress_window.destroy()

    # Ejecutar el script en un hilo separado para no bloquear la interfaz gráfica
    thread = Thread(target=run_script_and_update_output)
    thread.start()

def browse_path():
    directory = filedialog.askdirectory()
    path_var.set(directory)

def exit_app():
    root.quit()

def about_consenso_d():
    messagebox.showinfo("About CONSENSO_D", "CONSENSO_D Viral Genome Assembler\nVersion 1.0")

# Crear ventana principal
root = tk.Tk()
root.title("CONSENSO_D Viral Genome Assembler")

# Crear barra de menú
menubar = tk.Menu(root)
root.config(menu=menubar)

# Crear menú File
file_menu = tk.Menu(menubar, tearoff=0)
file_menu.add_command(label="Exit", command=exit_app)
menubar.add_cascade(label="File", menu=file_menu)

# Crear menú Help
help_menu = tk.Menu(menubar, tearoff=0)
help_menu.add_command(label="About CONSENSO_D", command=about_consenso_d)
menubar.add_cascade(label="Help", menu=help_menu)

# Subir imagen
image_path = '/home/ics2/CONSENSO_D/Resources/consenso logo.png'
image = Image.open(image_path)
display = ImageTk.PhotoImage(image)

# Crear y configurar el marco principal
main_frame = tk.Frame(root, padx=10, pady=10)
main_frame.pack(padx=10, pady=10)

# Añadir la imagen a la interfaz
image_label = tk.Label(main_frame, image=display)
image_label.image = display  # Guardar una referencia para evitar recolección de basura
image_label.grid(row=0, column=0, columnspan=3, pady=10)

# Etiqueta y campo desplegable para seleccionar el tipo de virus
tk.Label(main_frame, text="Select Virus Type:").grid(row=1, column=0, sticky="w")
virus_var = tk.StringVar()
virus_options = ["DENV_1", "DENV_2", "DENV_3", "DENV_4", "RABV"]
tk.OptionMenu(main_frame, virus_var, *virus_options).grid(row=1, column=1, sticky="ew")

# Etiqueta y campo para seleccionar el directorio
tk.Label(main_frame, text="Select Directory Path:").grid(row=2, column=0, sticky="w")
path_var = tk.StringVar()
tk.Entry(main_frame, textvariable=path_var, width=50).grid(row=2, column=1, sticky="ew")
tk.Button(main_frame, text="Browse", command=browse_path).grid(row=2, column=2, sticky="ew")

# Botón para ejecutar CONSENSO_D
tk.Button(main_frame, text="Run CONSENSO_D", command=run_consenso_d).grid(row=3, column=0, columnspan=3, pady=10)

# Ejecutar el bucle principal de la interfaz gráfica
root.mainloop()