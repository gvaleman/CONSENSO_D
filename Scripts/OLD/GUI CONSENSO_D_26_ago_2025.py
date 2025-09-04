import os
import tkinter as tk
from tkinter import filedialog, messagebox, Toplevel, Text, Scrollbar, ttk
from PIL import Image, ImageTk
import subprocess
from threading import Thread
from queue import Queue, Empty

class ConsensoGUI:
    def __init__(self):
        self.ventana = tk.Tk()
        self.ventana.title("CONSENSO_D: Viral Genome Assembler")
        
        # Variables
        self.sequence_type_var = tk.StringVar()
        self.virus_var = tk.StringVar()
        self.ruta_var = tk.StringVar()
        self.variant_serotipo_var = tk.StringVar()
        self.variant_fasta_var = tk.StringVar()
        
        # Mapeo de opciones
        self.opciones_secuenciacion_map = {"Nanopore": "NANO", "Illumina": "ILLUMINA"}
        
        self.setup_gui()
        
    def setup_gui(self):
        self.configurar_icono()
        self.crear_menu()
        self.crear_layout_principal()
        
    def configurar_icono(self):
        try:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            ruta_icono = os.path.join(current_dir, '..', 'Resources', 'icon.png')
            icono_imagen = Image.open(ruta_icono)
            icono = ImageTk.PhotoImage(icono_imagen)
            self.ventana.iconphoto(True, icono)
        except Exception as e:
            print(f"Error al cargar el icono: {e}")
    
    def crear_menu(self):
        menubar = tk.Menu(self.ventana)
        self.ventana.config(menu=menubar)
        
        menu_archivo = tk.Menu(menubar, tearoff=0)
        menu_archivo.add_command(label="Salir", command=self.salir_app)
        menubar.add_cascade(label="Archivo", menu=menu_archivo)
        
        menu_ayuda = tk.Menu(menubar, tearoff=0)
        menu_ayuda.add_command(label="Acerca de CONSENSO_D", command=self.acerca_de_consenso_d)
        menubar.add_cascade(label="Ayuda", menu=menu_ayuda)
    
    def crear_layout_principal(self):
        # Frame principal con padding
        main_frame = tk.Frame(self.ventana, padx=20, pady=20)
        main_frame.pack(fill='both', expand=True)
        
        # Logo en la parte superior
        self.agregar_logo(main_frame)
        
        # Frame horizontal para contenido
        content_frame = tk.Frame(main_frame)
        content_frame.pack(fill='both', expand=True, pady=20)
        
        # Panel izquierdo - Navegación
        self.crear_panel_navegacion(content_frame)
        
        # Panel derecho - Contenido
        self.crear_panel_contenido(content_frame)
    
    def agregar_logo(self, parent):
        try:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            ruta_imagen = os.path.join(current_dir, '..', 'Resources', 'consenso logo.png')
            imagen = Image.open(ruta_imagen)
            # Redimensionar imagen si es necesario
            imagen = imagen.resize((300, 100), Image.Resampling.LANCZOS)
            mostrar = ImageTk.PhotoImage(imagen)
            etiqueta_imagen = tk.Label(parent, image=mostrar)
            etiqueta_imagen.image = mostrar
            etiqueta_imagen.pack(pady=10)
        except Exception as e:
            print(f"Error al cargar la imagen: {e}")
    
    def crear_panel_navegacion(self, parent):
        # Frame izquierdo para navegación
        nav_frame = tk.Frame(parent, bg="#f0f0f0", relief="raised", bd=2)
        nav_frame.pack(side='left', fill='y', padx=(0, 20), pady=10)
        nav_frame.pack_propagate(False)
        nav_frame.config(width=250)
        
        # Título del panel
        titulo = tk.Label(nav_frame, text="FUNCIONALIDADES", 
                         font=("Helvetica", 14, "bold"), 
                         bg="#f0f0f0", fg="#2c3e50")
        titulo.pack(pady=20)
        
        # Separador
        separator = tk.Frame(nav_frame, height=2, bg="#3498db")
        separator.pack(fill='x', padx=20, pady=10)
        
        # Botones de navegación
        self.btn_ensamblaje = tk.Button(nav_frame, 
                                       text="🧬 Ensamblaje de Genomas",
                                       command=lambda: self.mostrar_panel("ensamblaje"),
                                       bg="#3498db", fg="white",
                                       font=("Helvetica", 11, "bold"),
                                       relief="flat", padx=20, pady=10,
                                       activebackground="#2980b9")
        self.btn_ensamblaje.pack(fill='x', padx=20, pady=5)
        
        self.btn_variantes = tk.Button(nav_frame,
                                      text="🧪 Llamado de Variantes", 
                                      command=lambda: self.mostrar_panel("variantes"),
                                      bg="#e74c3c", fg="white",
                                      font=("Helvetica", 11, "bold"),
                                      relief="flat", padx=20, pady=10,
                                      activebackground="#c0392b")
        self.btn_variantes.pack(fill='x', padx=20, pady=5)
        
        self.btn_lineajes = tk.Button(nav_frame,
                                     text="🌳 Viral-Branch (Linajes)",
                                     command=lambda: self.mostrar_panel("lineajes"),
                                     bg="#27ae60", fg="white",
                                     font=("Helvetica", 11, "bold"),
                                     relief="flat", padx=20, pady=10,
                                     activebackground="#229954")
        self.btn_lineajes.pack(fill='x', padx=20, pady=5)
        
        # Información del proyecto
        info_label = tk.Label(nav_frame, 
                             text="Instituto de Ciencias\nSostenibles", 
                             font=("Helvetica", 9, "italic"),
                             bg="#f0f0f0", fg="#7f8c8d")
        info_label.pack(side='bottom', pady=20)
    
    def crear_panel_contenido(self, parent):
        # Frame derecho para contenido
        self.content_frame = tk.Frame(parent, bg="white", relief="sunken", bd=2)
        self.content_frame.pack(side='right', fill='both', expand=True, pady=10)
        
        # Crear todos los paneles (inicialmente ocultos)
        self.crear_panel_ensamblaje()
        self.crear_panel_variantes() 
        self.crear_panel_lineajes()
        
        # Mostrar panel de ensamblaje por defecto
        self.mostrar_panel("ensamblaje")
    
    def crear_panel_ensamblaje(self):
        self.panel_ensamblaje = tk.Frame(self.content_frame, bg="white")
        
        # Título
        titulo = tk.Label(self.panel_ensamblaje, 
                         text="🧬 ENSAMBLAJE DE GENOMAS VIRALES", 
                         font=("Helvetica", 16, "bold"),
                         bg="white", fg="#2c3e50")
        titulo.pack(pady=20)
        
        # Frame para formulario
        form_frame = tk.Frame(self.panel_ensamblaje, bg="white")
        form_frame.pack(pady=20, padx=40, fill='x')
        
        # Tipo de secuenciación
        tk.Label(form_frame, text="Tipo de Secuenciación:", 
                font=("Helvetica", 11, "bold"), bg="white").grid(row=0, column=0, sticky="w", pady=10)
        opciones_secuenciacion_display = list(self.opciones_secuenciacion_map.keys())
        tk.OptionMenu(form_frame, self.sequence_type_var, *opciones_secuenciacion_display).grid(row=0, column=1, sticky="ew", padx=10)
        
        # Tipo de virus
        tk.Label(form_frame, text="Tipo de Virus:", 
                font=("Helvetica", 11, "bold"), bg="white").grid(row=1, column=0, sticky="w", pady=10)
        opciones_virus = ["DENV_1", "DENV_2", "DENV_3", "DENV_4", "SARS_COV_2", "RABV", "N_RABV"]
        tk.OptionMenu(form_frame, self.virus_var, *opciones_virus).grid(row=1, column=1, sticky="ew", padx=10)
        
        # Directorio
        tk.Label(form_frame, text="Directorio de Datos:", 
                font=("Helvetica", 11, "bold"), bg="white").grid(row=2, column=0, sticky="w", pady=10)
        tk.Entry(form_frame, textvariable=self.ruta_var, width=40).grid(row=2, column=1, sticky="ew", padx=10)
        tk.Button(form_frame, text="Buscar", command=self.buscar_ruta).grid(row=2, column=2, padx=10)
        
        # Configurar grid
        form_frame.columnconfigure(1, weight=1)
        
        # Botón ejecutar
        boton_ejecutar = tk.Button(self.panel_ensamblaje, 
                                  text="🚀 Ejecutar Ensamblaje", 
                                  command=self.ejecutar_consenso_d,
                                  bg="#3498db", fg="white",
                                  font=("Helvetica", 12, "bold"),
                                  padx=20, pady=10)
        boton_ejecutar.pack(pady=20)
    
    def crear_panel_variantes(self):
        self.panel_variantes = tk.Frame(self.content_frame, bg="white")
        
        # Título
        titulo = tk.Label(self.panel_variantes, 
                         text="🧪 LLAMADO DE VARIANTES", 
                         font=("Helvetica", 16, "bold"),
                         bg="white", fg="#e74c3c")
        titulo.pack(pady=20)
        
        # Descripción
        desc = tk.Label(self.panel_variantes, 
                       text="Generar archivos VCF a partir de secuencias consenso multiFASTA",
                       font=("Helvetica", 10, "italic"),
                       bg="white", fg="#7f8c8d")
        desc.pack(pady=10)
        
        # Frame para formulario
        form_frame = tk.Frame(self.panel_variantes, bg="white")
        form_frame.pack(pady=20, padx=40, fill='x')
        
        # Serotipo/Virus
        tk.Label(form_frame, text="Serotipo/Virus:", 
                font=("Helvetica", 11, "bold"), bg="white").grid(row=0, column=0, sticky="w", pady=10)
        opciones_variant = ["DENV_1", "DENV_2", "DENV_3", "DENV_4", "RABV"]
        tk.OptionMenu(form_frame, self.variant_serotipo_var, *opciones_variant).grid(row=0, column=1, sticky="ew", padx=10)
        
        # Archivo FASTA
        tk.Label(form_frame, text="Archivo multiFASTA:", 
                font=("Helvetica", 11, "bold"), bg="white").grid(row=1, column=0, sticky="w", pady=10)
        tk.Entry(form_frame, textvariable=self.variant_fasta_var, width=40).grid(row=1, column=1, sticky="ew", padx=10)
        tk.Button(form_frame, text="Buscar", command=self.buscar_fasta_variantes).grid(row=1, column=2, padx=10)
        
        # Configurar grid
        form_frame.columnconfigure(1, weight=1)
        
        # Botón ejecutar
        boton_ejecutar = tk.Button(self.panel_variantes, 
                                  text="🧪 Ejecutar Variant Calling", 
                                  command=self.ejecutar_variant_calling,
                                  bg="#e74c3c", fg="white",
                                  font=("Helvetica", 12, "bold"),
                                  padx=20, pady=10)
        boton_ejecutar.pack(pady=20)
    
    def crear_panel_lineajes(self):
        self.panel_lineajes = tk.Frame(self.content_frame, bg="white")
        
        # Título
        titulo = tk.Label(self.panel_lineajes, 
                         text="🌳 VIRAL-BRANCH (Determinación de Linajes)", 
                         font=("Helvetica", 16, "bold"),
                         bg="white", fg="#27ae60")
        titulo.pack(pady=20)
        
        # Descripción
        desc = tk.Label(self.panel_lineajes, 
                       text="Sistema de Machine Learning para determinación de linajes de dengue",
                       font=("Helvetica", 10, "italic"),
                       bg="white", fg="#7f8c8d")
        desc.pack(pady=10)
        
        # Mensaje de desarrollo
        dev_msg = tk.Label(self.panel_lineajes, 
                          text="🚧 Módulo en desarrollo\n\nEsta funcionalidad estará disponible próximamente",
                          font=("Helvetica", 12),
                          bg="white", fg="#f39c12")
        dev_msg.pack(pady=50)
    
    def mostrar_panel(self, panel_nombre):
        # Ocultar todos los paneles
        for widget in self.content_frame.winfo_children():
            widget.pack_forget()
        
        # Resetear colores de botones
        self.btn_ensamblaje.config(bg="#3498db", activebackground="#2980b9")
        self.btn_variantes.config(bg="#e74c3c", activebackground="#c0392b")  
        self.btn_lineajes.config(bg="#27ae60", activebackground="#229954")
        
        # Mostrar panel seleccionado y resaltar botón
        if panel_nombre == "ensamblaje":
            self.panel_ensamblaje.pack(fill='both', expand=True)
            self.btn_ensamblaje.config(bg="#2980b9", activebackground="#2980b9")
        elif panel_nombre == "variantes":
            self.panel_variantes.pack(fill='both', expand=True)
            self.btn_variantes.config(bg="#c0392b", activebackground="#c0392b")
        elif panel_nombre == "lineajes":
            self.panel_lineajes.pack(fill='both', expand=True)
            self.btn_lineajes.config(bg="#229954", activebackground="#229954")
    
    def buscar_ruta(self):
        directorio = filedialog.askdirectory()
        self.ruta_var.set(directorio)
    
    def buscar_fasta_variantes(self):
        archivo = filedialog.askopenfilename(
            title="Seleccionar archivo multiFASTA",
            filetypes=[("FASTA files", "*.fasta *.fa *.fas"), ("All files", "*.*")]
        )
        self.variant_fasta_var.set(archivo)
    
    def ejecutar_variant_calling(self):
        serotipo = self.variant_serotipo_var.get()
        archivo_fasta = self.variant_fasta_var.get()
        
        if not serotipo or not archivo_fasta:
            messagebox.showerror("Error", "Seleccione un serotipo y archivo FASTA.")
            return
        
        self.ejecutar_script_con_progreso("Variant_caller.sh", [serotipo, archivo_fasta], 
                                         "Variant Calling")
    
    def ejecutar_consenso_d(self):
        sequence_type_display = self.sequence_type_var.get()
        sequence_type = self.opciones_secuenciacion_map.get(sequence_type_display, "")
        virus = self.virus_var.get()
        ruta = self.ruta_var.get()

        if not sequence_type or not virus or not ruta:
            messagebox.showerror("Error", "Complete todos los campos.")
            return

        self.ejecutar_script_con_progreso("ASSEMBLER.sh", [sequence_type, virus, ruta], 
                                         "Ensamblaje")
    
    def ejecutar_script_con_progreso(self, script_name, args, titulo):
        ventana_progreso = Toplevel(self.ventana)
        ventana_progreso.title(f"Procesando {titulo}...")
        ventana_progreso.geometry("800x500")
        
        tk.Label(ventana_progreso, text=f"Ejecutando {titulo}...", 
                font=("Helvetica", 12, "bold")).pack(pady=10)

        area_texto = Text(ventana_progreso, wrap='word', height=25, width=100)
        area_texto.pack(side='left', fill='both', expand=True, padx=10)

        barra_desplazamiento = Scrollbar(ventana_progreso, command=area_texto.yview)
        barra_desplazamiento.pack(side='right', fill='y', padx=(0, 10))
        area_texto['yscrollcommand'] = barra_desplazamiento.set

        def ejecutar_script_thread():
            try:
                current_dir = os.path.dirname(os.path.abspath(__file__))
                script_path = os.path.join(current_dir, script_name)
                
                proceso = subprocess.Popen([script_path] + args,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT,
                                         text=True,
                                         universal_newlines=True)

                for line in iter(proceso.stdout.readline, ''):
                    area_texto.insert(tk.END, line)
                    area_texto.see(tk.END)
                    ventana_progreso.update()

                proceso.wait()
                
                if proceso.returncode == 0:
                    messagebox.showinfo("Éxito", f"¡{titulo} ejecutado con éxito!")
                else:
                    messagebox.showerror("Error", f"Error durante la ejecución de {titulo}.")
            except Exception as e:
                messagebox.showerror("Error", f"Error: {e}")
            finally:
                ventana_progreso.destroy()

        Thread(target=ejecutar_script_thread).start()
    
    def salir_app(self):
        self.ventana.quit()
    
    def acerca_de_consenso_d(self):
        messagebox.showinfo("Acerca de CONSENSO_D", 
                           "CONSENSO_D Viral Genome Assembler\n" +
                           "Versión 1.0\n\n" +
                           "Funcionalidades:\n" +
                           "• Ensamblaje de genomas virales\n" +
                           "• Llamado de variantes\n" +
                           "• Determinación de linajes (Viral-Branch)\n\n" +
                           "Instituto de Ciencias Sostenibles")
    
    def run(self):
        self.ventana.mainloop()

# Ejecutar la aplicación
if __name__ == "__main__":
    app = ConsensoGUI()
    app.run()
