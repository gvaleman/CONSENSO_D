import os
import tkinter as tk
from tkinter import filedialog, messagebox, Toplevel, Text, Scrollbar, ttk
from PIL import Image, ImageTk
import subprocess
from threading import Thread
from queue import Queue, Empty
import sys
import re
import shutil
import csv
import io

# Importar la función de organización de FASTQ
sys.path.append(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Scripts")
)
try:
    from organizar_fastq import organizar_fastq
except ImportError:
    # Definir la función aquí como respaldo
    def organizar_fastq(directorio=None):
        """
        Organiza archivos FASTQ en carpetas según el nombre de la muestra.
        """
        # Usar directorio actual si no se especifica
        if directorio is None:
            directorio = os.getcwd()

        # Cambiar al directorio especificado
        os.chdir(directorio)

        # Patrón para identificar archivos R1
        patron = re.compile(r"(.+?)_S\d+_L001_R1_001\.fastq\.gz")

        # Contador para estadísticas
        muestras_procesadas = 0
        archivos_movidos = 0

        print(f"Organizando archivos FASTQ en: {directorio}")

        # Recorrer todos los archivos en el directorio
        for archivo in os.listdir("."):
            # Buscar archivos R1
            coincidencia = patron.match(archivo)
            if coincidencia:
                nombre_muestra = coincidencia.group(1)
                archivo_r1 = archivo
                archivo_r2 = archivo.replace("_R1_", "_R2_")

                # Verificar si existe el archivo R2 correspondiente
                if os.path.isfile(archivo_r2):
                    print(f"Procesando muestra: {nombre_muestra}")

                    # Crear directorio si no existe
                    if not os.path.exists(nombre_muestra):
                        os.makedirs(nombre_muestra)

                    # Mover archivos R1 y R2 al directorio
                    shutil.move(archivo_r1, os.path.join(nombre_muestra, archivo_r1))
                    shutil.move(archivo_r2, os.path.join(nombre_muestra, archivo_r2))

                    print(f"  → Archivos movidos a {nombre_muestra}/")
                    muestras_procesadas += 1
                    archivos_movidos += 2
                else:
                    print(f"ADVERTENCIA: No se encontró archivo R2 para {archivo_r1}")

        print(
            f"Organización completada. {muestras_procesadas} muestras procesadas, {archivos_movidos} archivos movidos."
        )
        return muestras_procesadas, archivos_movidos


# Definir la clase SampleSheet aquí como respaldo
try:
    from sample_sheet import SampleSheet
except ImportError:

    class SampleSheet:
        def __init__(self, sample_sheet_file):
            self.sample_sheet_file = sample_sheet_file
            self.samples = self._parse_sample_sheet()

        def _parse_sample_sheet(self):
            samples = {}
            with open(self.sample_sheet_file, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    sample_id = row.get("Sample_ID")
                    sample_name = row.get("Sample_Name")
                    if sample_id and sample_name:
                        samples[sample_id] = sample_name
            return samples

        def rename_consensus_files(self, consensus_file, output_dir):
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            renamed_files = []
            if os.path.isfile(consensus_file):
                # Asumir que el nombre del archivo de consenso es el Sample_ID
                base_name = os.path.basename(consensus_file)
                sample_id = os.path.splitext(base_name)[0]

                if sample_id in self.samples:
                    new_name = f"{self.samples[sample_id]}.fasta"
                    output_path = os.path.join(output_dir, new_name)
                    shutil.copy(consensus_file, output_path)
                    renamed_files.append(output_path)
                    print(f"Renombrado '{consensus_file}' a '{new_name}'")
            else:
                print("El archivo de consenso no existe.")

            return renamed_files


class ConsensoGUI:
    def __init__(self):
        self.ventana = tk.Tk()
        self.ventana.title("CONSENSO_D: Viral Genome Assembler")

        # Variables existentes
        self.sequence_type_var = tk.StringVar()
        self.denovo_sequence_type_var = tk.StringVar()
        self.virus_var = tk.StringVar()
        self.ruta_var = tk.StringVar()
        self.variant_serotipo_var = tk.StringVar()
        self.variant_fasta_var = tk.StringVar()
        self.carpeta_fastq_var = tk.StringVar()

        # Nuevas variables
        self.sample_sheet_var = tk.StringVar()
        self.consensus_file_var = tk.StringVar()
        self.output_dir_var = tk.StringVar()
        self.total_samples = 0
        self.processed_samples = 0
        self.progress_var = tk.DoubleVar()
        self.progress_label_var = tk.StringVar()

        # Nuevas variables para primers
        self.primer_option_var = tk.StringVar(value="None")
        self.primer_protocol_var = tk.StringVar()
        self.custom_primer_file_var = tk.StringVar()

        # Nuevas variables para el conversor de CSV a BED
        self.csv_converter_input_var = tk.StringVar()
        self.csv_converter_virus_var = tk.StringVar()

        # Mapeo de opciones
        self.opciones_secuenciacion_map = {"Nanopore": "NANO", "Illumina": "ILLUMINA"}
        self.opciones_denovo_map = {"Illumina": "illumina", "Nanopore": "nanopore"}

        self.setup_gui()

    def setup_gui(self):
        self.configurar_icono()
        self.crear_menu()
        self.crear_layout_principal()

    def configurar_icono(self):
        try:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            ruta_icono = os.path.join(current_dir, "..", "Resources", "icon.png")
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
        menu_ayuda.add_command(
            label="Acerca de CONSENSO_D", command=self.acerca_de_consenso_d
        )
        menubar.add_cascade(label="Ayuda", menu=menu_ayuda)

    def crear_layout_principal(self):
        # Frame principal con padding
        main_frame = tk.Frame(self.ventana, padx=20, pady=20)
        main_frame.pack(fill="both", expand=True)

        # Logo en la parte superior
        self.agregar_logo(main_frame)

        # Frame horizontal para contenido
        content_frame = tk.Frame(main_frame)
        content_frame.pack(fill="both", expand=True, pady=20)

        # Panel izquierdo - Navegación
        self.crear_panel_navegacion(content_frame)

        # Panel derecho - Contenido
        self.crear_panel_contenido(content_frame)

    def agregar_logo(self, parent):
        try:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            ruta_imagen = os.path.join(
                current_dir, "..", "Resources", "consenso logo.png"
            )
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
        nav_frame.pack(side="left", fill="y", padx=(0, 20), pady=10)
        nav_frame.pack_propagate(False)
        nav_frame.config(width=250)

        # Título del panel
        titulo = tk.Label(
            nav_frame,
            text="FUNCIONALIDADES",
            font=("Helvetica", 14, "bold"),
            bg="#f0f0f0",
            fg="#2c3e50",
        )
        titulo.pack(pady=20)

        # Separador
        separator = tk.Frame(nav_frame, height=2, bg="#3498db")
        separator.pack(fill="x", padx=20, pady=10)

        # Botones de navegación
        self.btn_ensamblaje = tk.Button(
            nav_frame,
            text="🧬 Ensamblaje por Referencia",
            command=lambda: self.mostrar_panel("ensamblaje"),
            bg="#3498db",
            fg="white",
            font=("Helvetica", 11, "bold"),
            relief="flat",
            padx=20,
            pady=10,
            activebackground="#2980b9",
        )
        self.btn_ensamblaje.pack(fill="x", padx=20, pady=5)

        self.btn_de_novo = tk.Button(
            nav_frame,
            text="🔬 Ensamblaje De Novo",
            command=lambda: self.mostrar_panel("de_novo"),
            bg="#16a085",
            fg="white",
            font=("Helvetica", 11, "bold"),
            relief="flat",
            padx=20,
            pady=10,
            activebackground="#117a65",
        )
        self.btn_de_novo.pack(fill="x", padx=20, pady=5)

        self.btn_variantes = tk.Button(
            nav_frame,
            text="🧪 Llamado de Variantes",
            command=lambda: self.mostrar_panel("variantes"),
            bg="#f39c12",
            fg="white",
            font=("Helvetica", 11, "bold"),
            relief="flat",
            padx=20,
            pady=10,
            activebackground="#d35400",
        )
        self.btn_variantes.pack(fill="x", padx=20, pady=5)

        self.btn_lineajes = tk.Button(
            nav_frame,
            text="🌳 Viral-Branch (Linajes)",
            command=lambda: self.mostrar_panel("lineajes"),
            bg="#8e44ad",
            fg="white",
            font=("Helvetica", 11, "bold"),
            relief="flat",
            padx=20,
            pady=10,
            activebackground="#7d3c98",
        )
        self.btn_lineajes.pack(fill="x", padx=20, pady=5)

        # Nuevo botón para Otras Funciones
        self.btn_otras = tk.Button(
            nav_frame,
            text="🔧 Otras Funciones",
            command=lambda: self.mostrar_panel("otras"),
            bg="#7f8c8d",
            fg="white",
            font=("Helvetica", 11, "bold"),
            relief="flat",
            padx=20,
            pady=10,
            activebackground="#566573",
        )
        self.btn_otras.pack(fill="x", padx=20, pady=5)

        # Información del proyecto
        info_label = tk.Label(
            nav_frame,
            text="Instituto de Ciencias\nSostenibles",
            font=("Helvetica", 9, "italic"),
            bg="#f0f0f0",
            fg="#7f8c8d",
        )
        info_label.pack(side="bottom", pady=20)

    def crear_panel_contenido(self, parent):
        # Frame derecho para contenido
        self.content_frame = tk.Frame(parent, bg="white", relief="sunken", bd=2)
        self.content_frame.pack(side="right", fill="both", expand=True, pady=10)

        # Crear todos los paneles (inicialmente ocultos)
        self.crear_panel_ensamblaje()
        self.crear_panel_de_novo()
        self.crear_panel_variantes()
        self.crear_panel_lineajes()
        self.crear_panel_otras_funciones()

        # Mostrar panel de ensamblaje por defecto
        self.mostrar_panel("ensamblaje")

    def crear_panel_de_novo(self):
        self.panel_de_novo = tk.Frame(self.content_frame, bg="white")

        # Título
        titulo = tk.Label(
            self.panel_de_novo,
            text="🔬 ENSAMBLAJE DE NOVO",
            font=("Helvetica", 16, "bold"),
            bg="white",
            fg="#16a085",
        )
        titulo.pack(pady=20)

        # Frame para formulario
        form_frame = tk.Frame(self.panel_de_novo, bg="white")
        form_frame.pack(pady=20, padx=40, fill="x")

        # Directorio
        tk.Label(
            form_frame,
            text="Directorio de Datos:",
            font=("Helvetica", 11, "bold"),
            bg="white",
        ).grid(row=1, column=0, sticky="w", pady=10)
        tk.Entry(form_frame, textvariable=self.ruta_var, width=40).grid(
            row=1, column=1, sticky="ew", padx=10
        )
        tk.Button(form_frame, text="Buscar", command=self.buscar_directorio).grid(
            row=1, column=2, padx=10
        )

        # Configurar el grid para expandir la entrada
        form_frame.columnconfigure(1, weight=1)

        # Botón para ejecutar el ensamblaje De Novo
        boton_de_novo = tk.Button(
            self.panel_de_novo,
            text="🧬 Ejecutar Ensamblaje De Novo",
            command=self.ejecutar_de_novo_ensamblaje,
            bg="#16a085",
            fg="white",
            font=("Helvetica", 12, "bold"),
            padx=20,
            pady=10,
        )
        boton_de_novo.pack(pady=30)

        # Opciones de tecnología de secuenciación
        tk.Label(
            form_frame,
            text="Tecnología de Secuenciación:",
            font=("Helvetica", 11, "bold"),
            bg="white",
        ).grid(row=2, column=0, sticky="w", pady=10)
        opciones_denovo_tech_display = ["illumina", "nanopore"]
        self.denovo_sequence_type_var.set(
            opciones_denovo_tech_display[0]
        )  # Set default
        tk.OptionMenu(
            form_frame, self.denovo_sequence_type_var, *opciones_denovo_tech_display
        ).grid(row=2, column=1, sticky="ew", padx=10)

    def crear_panel_ensamblaje(self):
        self.panel_ensamblaje = tk.Frame(self.content_frame, bg="white")

        # Título
        titulo = tk.Label(
            self.panel_ensamblaje,
            text="🧬 ENSAMBLAJE POR REFERENCIA",
            font=("Helvetica", 16, "bold"),
            bg="white",
            fg="#3498db",
        )
        titulo.pack(pady=20)

        # Frame para formulario
        form_frame = tk.Frame(self.panel_ensamblaje, bg="white")
        form_frame.pack(pady=20, padx=40, fill="x")

        # Tipo de secuenciación
        tk.Label(
            form_frame,
            text="Tipo de Secuenciación:",
            font=("Helvetica", 11, "bold"),
            bg="white",
        ).grid(row=0, column=0, sticky="w", pady=10)
        opciones_secuenciacion_display = list(self.opciones_secuenciacion_map.keys())
        tk.OptionMenu(
            form_frame, self.sequence_type_var, *opciones_secuenciacion_display
        ).grid(row=0, column=1, sticky="ew", padx=10)

        # Tipo de virus
        tk.Label(
            form_frame,
            text="Tipo de Virus:",
            font=("Helvetica", 11, "bold"),
            bg="white",
        ).grid(row=1, column=0, sticky="w", pady=10)
        opciones_virus = [
            "DENV_1",
            "DENV_2",
            "DENV_3",
            "DENV_4",
            "DENV-AUTO",
            "SARS_COV_2",
            "RABV",
            "N_RABV",
        ]
        tk.OptionMenu(form_frame, self.virus_var, *opciones_virus).grid(
            row=1, column=1, sticky="ew", padx=10
        )

        # Selección de Primers
        tk.Label(
            form_frame, text="Primers:", font=("Helvetica", 11, "bold"), bg="white"
        ).grid(row=2, column=0, sticky="w", pady=10)

        primer_options_frame = tk.Frame(form_frame, bg="white")
        primer_options_frame.grid(row=2, column=1, sticky="ew", padx=10)

        tk.Radiobutton(
            primer_options_frame,
            text="None",
            variable=self.primer_option_var,
            value="None",
            bg="white",
            command=self.toggle_primer_options,
        ).pack(side="left", padx=5)
        tk.Radiobutton(
            primer_options_frame,
            text="Pre-defined",
            variable=self.primer_option_var,
            value="Pre-defined",
            bg="white",
            command=self.toggle_primer_options,
        ).pack(side="left", padx=5)
        tk.Radiobutton(
            primer_options_frame,
            text="Custom",
            variable=self.primer_option_var,
            value="Custom",
            bg="white",
            command=self.toggle_primer_options,
        ).pack(side="left", padx=5)

        # Opciones para primers pre-definidos
        self.predefined_primer_frame = tk.Frame(form_frame, bg="white")
        self.predefined_primer_frame.grid(row=3, column=1, sticky="ew", padx=10, pady=5)

        tk.Label(
            self.predefined_primer_frame,
            text="Protocolo:",
            font=("Helvetica", 10),
            bg="white",
        ).pack(side="left")
        current_dir = os.path.dirname(os.path.abspath(__file__))
        primer_dir = os.path.join(current_dir, "..", "Primers")
        protocolos = [
            d
            for d in os.listdir(primer_dir)
            if os.path.isdir(os.path.join(primer_dir, d))
        ]
        self.primer_protocol_menu = tk.OptionMenu(
            self.predefined_primer_frame, self.primer_protocol_var, *protocolos
        )
        self.primer_protocol_menu.pack(side="left", padx=5)

        # Opciones para primers custom
        self.custom_primer_frame = tk.Frame(form_frame, bg="white")
        self.custom_primer_frame.grid(row=3, column=1, sticky="ew", padx=10, pady=5)

        tk.Entry(
            self.custom_primer_frame, textvariable=self.custom_primer_file_var, width=30
        ).pack(side="left", expand=True, fill="x")
        tk.Button(
            self.custom_primer_frame,
            text="Buscar",
            command=self.buscar_custom_primer_file,
        ).pack(side="left", padx=5)

        # Directorio
        tk.Label(
            form_frame,
            text="Directorio de Datos:",
            font=("Helvetica", 11, "bold"),
            bg="white",
        ).grid(row=4, column=0, sticky="w", pady=10)
        tk.Entry(form_frame, textvariable=self.ruta_var, width=40).grid(
            row=4, column=1, sticky="ew", padx=10
        )
        tk.Button(form_frame, text="Buscar", command=self.buscar_directorio).grid(
            row=4, column=2, padx=10
        )

        # Configurar el grid para expandir la entrada
        form_frame.columnconfigure(1, weight=1)

        # Botón para ejecutar el ensamblaje
        boton_ensamblar = tk.Button(
            self.panel_ensamblaje,
            text="🚀 Ejecutar Ensamblaje de Referencia",
            command=self.ejecutar_ensamblaje,
            bg="#3498db",
            fg="white",
            font=("Helvetica", 12, "bold"),
            padx=20,
            pady=10,
        )
        boton_ensamblar.pack(pady=30)

        self.toggle_primer_options()  # Inicializar visibilidad

    def toggle_primer_options(self):
        option = self.primer_option_var.get()
        if option == "Pre-defined":
            self.predefined_primer_frame.grid()
            self.custom_primer_frame.grid_remove()
        elif option == "Custom":
            self.predefined_primer_frame.grid_remove()
            self.custom_primer_frame.grid()
        else:  # None
            self.predefined_primer_frame.grid_remove()
            self.custom_primer_frame.grid_remove()

    def buscar_custom_primer_file(self):
        archivo = filedialog.askopenfilename(
            title="Seleccionar archivo de primers FASTA",
            filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")],
        )
        self.custom_primer_file_var.set(archivo)

    def crear_panel_variantes(self):
        self.panel_variantes = tk.Frame(self.content_frame, bg="white")
        # Contenido de Llamado de Variantes
        titulo = tk.Label(
            self.panel_variantes,
            text="🧪 LLAMADO DE VARIANTES",
            font=("Helvetica", 16, "bold"),
            bg="white",
            fg="#f39c12",
        )
        titulo.pack(pady=20)

        form_frame = tk.Frame(self.panel_variantes, bg="white")
        form_frame.pack(pady=20, padx=40, fill="x")

        # Archivo multiFASTA
        tk.Label(
            form_frame,
            text="Archivo multiFASTA:",
            font=("Helvetica", 11, "bold"),
            bg="white",
        ).grid(row=0, column=0, sticky="w", pady=10)
        tk.Entry(form_frame, textvariable=self.variant_fasta_var, width=40).grid(
            row=0, column=1, sticky="ew", padx=10
        )
        tk.Button(form_frame, text="Buscar", command=self.buscar_fasta_variantes).grid(
            row=0, column=2, padx=10
        )

        # Serotipo para variantes
        tk.Label(
            form_frame,
            text="Tipo de Virus:",
            font=("Helvetica", 11, "bold"),
            bg="white",
        ).grid(row=1, column=0, sticky="w", pady=10)
        opciones_virus_variantes = [
            "DENV_1",
            "DENV_2",
            "DENV_3",
            "DENV_4",
            "SARS_COV_2",
            "RABV",
            "N_RABV",
        ]
        tk.OptionMenu(
            form_frame, self.variant_serotipo_var, *opciones_virus_variantes
        ).grid(row=1, column=1, sticky="ew", padx=10)

        form_frame.columnconfigure(1, weight=1)

        boton_variantes = tk.Button(
            self.panel_variantes,
            text="🔬 Ejecutar Llamado de Variantes",
            command=self.ejecutar_llamado_variantes,
            bg="#f39c12",
            fg="white",
            font=("Helvetica", 12, "bold"),
            padx=20,
            pady=10,
        )
        boton_variantes.pack(pady=30)

    def crear_panel_lineajes(self):
        self.panel_lineajes = tk.Frame(self.content_frame, bg="white")
        # Contenido de Viral-Branch (Linajes)
        titulo = tk.Label(
            self.panel_lineajes,
            text="🌳 VIRAL-BRANCH (LINAJES)",
            font=("Helvetica", 16, "bold"),
            bg="white",
            fg="#8e44ad",
        )
        titulo.pack(pady=20)

        form_frame = tk.Frame(self.panel_lineajes, bg="white")
        form_frame.pack(pady=20, padx=40, fill="x")

        # Directorio de FASTQ
        tk.Label(
            form_frame,
            text="Directorio de Archivos FASTQ:",
            font=("Helvetica", 11, "bold"),
            bg="white",
        ).grid(row=0, column=0, sticky="w", pady=10)
        tk.Entry(form_frame, textvariable=self.carpeta_fastq_var, width=40).grid(
            row=0, column=1, sticky="ew", padx=10
        )
        tk.Button(form_frame, text="Buscar", command=self.buscar_directorio_fastq).grid(
            row=0, column=2, padx=10
        )

        form_frame.columnconfigure(1, weight=1)

        boton_lineajes = tk.Button(
            self.panel_lineajes,
            text="📊 Ejecutar Viral-Branch",
            command=self.ejecutar_viral_branch,
            bg="#8e44ad",
            fg="white",
            font=("Helvetica", 12, "bold"),
            padx=20,
            pady=10,
        )
        boton_lineajes.pack(pady=30)

    def crear_panel_otras_funciones(self):
        self.panel_otras = tk.Frame(self.content_frame, bg="white")

        # Título de la sección
        titulo = tk.Label(
            self.panel_otras,
            text="🔧 OTRAS FUNCIONES",
            font=("Helvetica", 16, "bold"),
            bg="white",
            fg="#7f8c8d",
        )
        titulo.pack(pady=20)

        # Usar un Notebook para la interfaz con pestañas
        self.notebook_otras = ttk.Notebook(self.panel_otras)
        self.notebook_otras.pack(pady=10, expand=True, fill="both")

        # Pestaña para 'Organizar FASTQ'
        tab_organizar = tk.Frame(self.notebook_otras, bg="white")
        self.notebook_otras.add(tab_organizar, text="Organizar FASTQ")
        self.crear_panel_organizar_fastq(tab_organizar)

        # Pestaña para 'Procesar Sample Sheet'
        tab_sample_sheet = tk.Frame(self.notebook_otras, bg="white")
        self.notebook_otras.add(tab_sample_sheet, text="Procesar Sample Sheet")
        self.crear_panel_sample_sheet(tab_sample_sheet)

        # Pestaña para 'CSV to BED Converter'
        tab_csv_to_bed = tk.Frame(self.notebook_otras, bg="white")
        self.notebook_otras.add(tab_csv_to_bed, text="CSV to BED Converter")
        self.crear_panel_csv_to_bed(tab_csv_to_bed)

    def crear_panel_organizar_fastq(self, parent):
        panel_organizar = parent

        # Título de la sección
        seccion_title = tk.Label(
            panel_organizar,
            text="Organizar Archivos FASTQ",
            font=("Helvetica", 12, "bold"),
            bg="white",
            fg="#9b59b6",
        )
        seccion_title.pack(pady=(20, 10))

        # Descripción
        desc = tk.Label(
            panel_organizar,
            text="Organiza archivos .fastq.gz de Illumina en carpetas de muestra",
            font=("Helvetica", 10, "italic"),
            bg="white",
            fg="#7f8c8d",
        )
        desc.pack(pady=(0, 15))

        # Frame para formulario
        form_frame = tk.Frame(panel_organizar, bg="white")
        form_frame.pack(fill="x", padx=20)

        # Directorio
        tk.Label(
            form_frame,
            text="Directorio a Organizar:",
            font=("Helvetica", 10, "bold"),
            bg="white",
        ).grid(row=0, column=0, sticky="w", pady=10)
        tk.Entry(form_frame, textvariable=self.carpeta_fastq_var, width=40).grid(
            row=0, column=1, sticky="ew", padx=10
        )
        tk.Button(form_frame, text="Buscar", command=self.buscar_directorio_fastq).grid(
            row=0, column=2, padx=10
        )

        form_frame.columnconfigure(1, weight=1)

        # Botón ejecutar
        boton_organizar = tk.Button(
            panel_organizar,
            text="📂 Organizar FASTQ",
            command=self.ejecutar_organizar_fastq,
            bg="#9b59b6",
            fg="white",
            font=("Helvetica", 11, "bold"),
            padx=15,
            pady=8,
        )
        boton_organizar.pack(pady=20)

    # Añadir esta función para crear el panel de sample sheet
    def crear_panel_sample_sheet(self, parent):
        self.panel_sample_sheet = parent

        # Título de la sección
        seccion_title = tk.Label(
            self.panel_sample_sheet,
            text="Procesar Sample Sheet de Illumina",
            font=("Helvetica", 12, "bold"),
            bg="white",
            fg="#3498db",
        )
        seccion_title.pack(pady=(20, 10))

        # Descripción
        desc = tk.Label(
            self.panel_sample_sheet,
            text="Renombra archivos de consenso usando información del sample sheet",
            font=("Helvetica", 10, "italic"),
            bg="white",
            fg="#7f8c8d",
        )
        desc.pack(pady=(0, 15))

        # Frame para formulario
        form_frame = tk.Frame(self.panel_sample_sheet, bg="white")
        form_frame.pack(fill="x", padx=20)

        # Sample Sheet
        tk.Label(
            form_frame,
            text="Sample Sheet (CSV):",
            font=("Helvetica", 10, "bold"),
            bg="white",
        ).grid(row=0, column=0, sticky="w", pady=10)
        tk.Entry(form_frame, textvariable=self.sample_sheet_var, width=40).grid(
            row=0, column=1, sticky="ew", padx=10
        )
        tk.Button(form_frame, text="Buscar", command=self.buscar_sample_sheet).grid(
            row=0, column=2, padx=10
        )

        # Archivo de consenso
        tk.Label(
            form_frame,
            text="Archivo de consenso:",
            font=("Helvetica", 10, "bold"),
            bg="white",
        ).grid(row=1, column=0, sticky="w", pady=10)
        tk.Entry(form_frame, textvariable=self.consensus_file_var, width=40).grid(
            row=1, column=1, sticky="ew", padx=10
        )
        tk.Button(form_frame, text="Buscar", command=self.buscar_archivo_consenso).grid(
            row=1, column=2, padx=10
        )

        # Directorio de salida
        tk.Label(
            form_frame,
            text="Directorio de salida:",
            font=("Helvetica", 10, "bold"),
            bg="white",
        ).grid(row=2, column=0, sticky="w", pady=10)
        tk.Entry(form_frame, textvariable=self.output_dir_var, width=40).grid(
            row=2, column=1, sticky="ew", padx=10
        )
        tk.Button(form_frame, text="Buscar", command=self.buscar_dir_salida).grid(
            row=2, column=2, padx=10
        )

        # Configurar grid
        form_frame.columnconfigure(1, weight=1)

        # Botón ejecutar
        boton_procesar = tk.Button(
            self.panel_sample_sheet,
            text="🧬 Procesar Sample Sheet",
            command=self.procesar_sample_sheet,
            bg="#3498db",
            fg="white",
            font=("Helvetica", 11, "bold"),
            padx=15,
            pady=8,
        )
        boton_procesar.pack(pady=20)

    def crear_panel_csv_to_bed(self, parent):
        panel_converter = parent

        # Título de la sección
        seccion_title = tk.Label(
            panel_converter,
            text="Convertidor de Primers CSV a BED",
            font=("Helvetica", 12, "bold"),
            bg="white",
            fg="#16a085",
        )
        seccion_title.pack(pady=(20, 10))

        # Descripción
        desc = tk.Label(
            panel_converter,
            text="Convierte un archivo CSV de primers (nombre,secuencia) a formato BED alineando contra un genoma de referencia.",
            font=("Helvetica", 10, "italic"),
            bg="white",
            fg="#7f8c8d",
            wraplength=500,
        )
        desc.pack(pady=(0, 15))

        # Frame para formulario
        form_frame = tk.Frame(panel_converter, bg="white")
        form_frame.pack(fill="x", padx=20)

        # Archivo CSV de entrada
        tk.Label(
            form_frame,
            text="Archivo CSV de Primers:",
            font=("Helvetica", 10, "bold"),
            bg="white",
        ).grid(row=0, column=0, sticky="w", pady=10)
        tk.Entry(form_frame, textvariable=self.csv_converter_input_var, width=40).grid(
            row=0, column=1, sticky="ew", padx=10
        )
        tk.Button(form_frame, text="Buscar", command=self.buscar_csv_file).grid(
            row=0, column=2, padx=10
        )

        # Tipo de virus para referencia
        tk.Label(
            form_frame,
            text="Genoma de Referencia (Virus):",
            font=("Helvetica", 10, "bold"),
            bg="white",
        ).grid(row=1, column=0, sticky="w", pady=10)
        opciones_virus = [
            "DENV_1",
            "DENV_2",
            "DENV_3",
            "DENV_4",
            "SARS_COV_2",
            "RABV",
            "N_RABV",
        ]
        tk.OptionMenu(form_frame, self.csv_converter_virus_var, *opciones_virus).grid(
            row=1, column=1, sticky="ew", padx=10
        )

        form_frame.columnconfigure(1, weight=1)

        # Frame para botones de acción
        action_frame = tk.Frame(panel_converter, bg="white")
        action_frame.pack(pady=20)

        # Botón ejecutar
        boton_convertir = tk.Button(
            action_frame,
            text="🧬 Convertir a BED",
            command=self.ejecutar_csv_to_bed,
            bg="#16a085",
            fg="white",
            font=("Helvetica", 11, "bold"),
            padx=15,
            pady=8,
        )
        boton_convertir.pack(side="left", padx=10)

        # Botón para descargar plantilla
        boton_plantilla = tk.Button(
            action_frame,
            text="📄 Descargar Ejemplo",
            command=self.descargar_plantilla_csv,
            bg="#7f8c8d",
            fg="white",
            font=("Helvetica", 11, "bold"),
            padx=15,
            pady=8,
        )
        boton_plantilla.pack(side="left", padx=10)

    def add_to_log(self, text_widget, message):
        """Helper function to add messages to a text widget log."""
        text_widget.config(state="normal")
        text_widget.insert(tk.END, message)
        text_widget.see(tk.END)  # Auto-scroll
        text_widget.config(state="disabled")

    def descargar_plantilla_csv(self):
        try:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            template_path = os.path.join(
                current_dir, "..", "Resources", "primer_template.csv"
            )

            if not os.path.exists(template_path):
                messagebox.showerror("Error", "No se encontró el archivo de plantilla.")
                return

            save_path = filedialog.asksaveasfilename(
                title="Guardar plantilla CSV como...",
                defaultextension=".csv",
                initialfile="primer_template.csv",
                filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            )

            if save_path:
                shutil.copy(template_path, save_path)
                messagebox.showinfo("Éxito", f"Plantilla guardada en: {save_path}")
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo guardar la plantilla: {e}")

    def buscar_csv_file(self):
        archivo = filedialog.askopenfilename(
            title="Seleccionar archivo CSV de primers",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        self.csv_converter_input_var.set(archivo)

    def ejecutar_csv_to_bed(self):
        csv_file = self.csv_converter_input_var.get()
        virus_type = self.csv_converter_virus_var.get()

        if not csv_file or not virus_type:
            messagebox.showerror(
                "Error",
                "Por favor, seleccione un archivo CSV y un tipo de virus de referencia.",
            )
            return

        # Crear ventana de progreso
        ventana_progreso = Toplevel(self.ventana)
        ventana_progreso.title("Convirtiendo CSV a BED...")
        ventana_progreso.geometry("700x400")

        tk.Label(
            ventana_progreso,
            text="Procesando archivo de primers...",
            font=("Helvetica", 12, "bold"),
        ).pack(pady=10)

        area_texto = Text(
            ventana_progreso, wrap="word", height=20, width=80, state="disabled"
        )
        area_texto.pack(side="left", fill="both", expand=True, padx=10)

        barra_desplazamiento = Scrollbar(ventana_progreso, command=area_texto.yview)
        barra_desplazamiento.pack(side="right", fill="y", padx=(0, 10))
        area_texto["yscrollcommand"] = barra_desplazamiento.set

        def ejecutar_thread():
            try:
                current_dir = os.path.dirname(os.path.abspath(__file__))
                script_path = os.path.join(current_dir, "csv_to_bed.py")

                if not os.path.exists(script_path):
                    self.ventana.after(
                        0,
                        lambda: messagebox.showerror(
                            "Error",
                            f"El script de conversión no se encuentra: {script_path}",
                        ),
                    )
                    ventana_progreso.destroy()
                    return

                command = ["python3", script_path, csv_file, virus_type]
                proceso = subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    encoding="utf-8",
                )

                def leer_salida():
                    for line in iter(proceso.stdout.readline, ""):
                        self.ventana.after(0, self.add_to_log, area_texto, line)
                    for line in iter(proceso.stderr.readline, ""):
                        self.ventana.after(0, self.add_to_log, area_texto, line)

                Thread(target=leer_salida, daemon=True).start()
                proceso.wait()

                def update_gui():
                    if proceso.returncode == 0:
                        output_bed = os.path.splitext(csv_file)[0] + ".bed"
                        messagebox.showinfo(
                            "Éxito",
                            f"Conversión completada.\nArchivo guardado en: {output_bed}",
                        )
                    else:
                        messagebox.showerror(
                            "Error",
                            f"La conversión falló. Revise la ventana de logs para más detalles.",
                        )
                    ventana_progreso.destroy()

                self.ventana.after(0, update_gui)

            except Exception as e:

                def show_error():
                    messagebox.showerror("Error", f"Ocurrió un error inesperado: {e}")
                    ventana_progreso.destroy()

                self.ventana.after(0, show_error)

        Thread(target=ejecutar_thread).start()

    def mostrar_panel(self, panel_name):
        # Oculta todos los paneles de contenido
        self.panel_ensamblaje.pack_forget()
        self.panel_de_novo.pack_forget()
        self.panel_variantes.pack_forget()
        self.panel_lineajes.pack_forget()
        self.panel_otras.pack_forget()

        # Muestra el panel seleccionado
        if panel_name == "ensamblaje":
            self.panel_ensamblaje.pack(fill="both", expand=True)
        elif panel_name == "de_novo":
            self.panel_de_novo.pack(fill="both", expand=True)
        elif panel_name == "variantes":
            self.panel_variantes.pack(fill="both", expand=True)
        elif panel_name == "lineajes":
            self.panel_lineajes.pack(fill="both", expand=True)
        elif panel_name == "otras":
            self.panel_otras.pack(fill="both", expand=True)

    def buscar_directorio(self):
        directorio = filedialog.askdirectory(title="Seleccionar directorio de datos")
        self.ruta_var.set(directorio)

    def buscar_fasta_variantes(self):
        archivo = filedialog.askopenfilename(
            title="Seleccionar archivo multiFASTA",
            filetypes=[("FASTA files", "*.fasta *.fa *.fas"), ("All files", "*.*")],
        )
        self.variant_fasta_var.set(archivo)

    def buscar_directorio_fastq(self):
        directorio = filedialog.askdirectory(
            title="Seleccionar directorio de archivos FASTQ"
        )
        self.carpeta_fastq_var.set(directorio)

    def ejecutar_ensamblaje(self):
        sequence_type_display = self.sequence_type_var.get()
        sequence_type_value = self.opciones_secuenciacion_map.get(sequence_type_display)
        virus = self.virus_var.get()
        ruta = self.ruta_var.get()
        primer_option = self.primer_option_var.get()

        # Validaciones iniciales
        if not sequence_type_value or not virus or not ruta:
            messagebox.showerror(
                "Error", "Seleccione el tipo de secuenciación, virus y directorio."
            )
            return

        # Lógica para determinar el script y los argumentos
        is_auto_denv = virus == "DENV-AUTO"

        if is_auto_denv:
            script_name = "ASSEMBLER_DENV_AUTO.sh"
        else:
            # Mantener los nombres de script anteriores para otros virus
            script_name = "ASSEMBLER.sh"  # O 'ASSEMBLERV1.sh' si se prefiere

        current_dir = os.path.dirname(os.path.abspath(__file__))
        script_path = os.path.join(current_dir, script_name)

        if not os.path.exists(script_path):
            # Fallback para el nombre alternativo del script original
            if not is_auto_denv and script_name == "ASSEMBLER.sh":
                script_path = os.path.join(current_dir, "ASSEMBLERV1.sh")

        if not os.path.exists(script_path):
            messagebox.showerror(
                "Error",
                f"El script '{script_name}' no se encuentra en el directorio de scripts.",
            )
            return

        primer_file_path = "none"
        if primer_option == "Pre-defined":
            protocol = self.primer_protocol_var.get()
            if not protocol:
                messagebox.showerror("Error", "Seleccione un protocolo de primers.")
                return

            # Para DENV-AUTO, no podemos saber el serotipo de antemano.
            # El usuario debe seleccionar un archivo de primers custom que contenga todos los necesarios
            # o el script de ensamblaje debe manejarlo. Por ahora, se deshabilita para DENV-AUTO.
            if is_auto_denv:
                messagebox.showwarning(
                    "Advertencia",
                    "La opción 'Pre-defined' no está disponible para DENV-AUTO. Por favor, use 'Custom' y seleccione un archivo BED o FASTA que contenga los primers para todos los serotipos.",
                )
                return

            virus_filename = virus.replace("_", "") + ".fasta"
            primer_file_path = os.path.join(
                current_dir, "..", "Primers", protocol, virus_filename
            )

            if not os.path.exists(primer_file_path):
                messagebox.showerror(
                    "Error", f"El archivo de primers no existe: {primer_file_path}"
                )
                return

        elif primer_option == "Custom":
            primer_file_path = self.custom_primer_file_var.get()
            if not primer_file_path or not os.path.exists(primer_file_path):
                messagebox.showerror(
                    "Error", "Seleccione un archivo de primers válido."
                )
                return

        # Construir el comando
        command = ["bash", script_path, sequence_type_value]
        if is_auto_denv:
            # El script AUTO no necesita el argumento 'virus'
            command.extend([ruta, primer_file_path])
        else:
            # El script original sí lo necesita
            command.extend([virus, ruta, primer_file_path])

        # Contar el número total de muestras
        try:
            self.total_samples = len(
                [d for d in os.listdir(ruta) if os.path.isdir(os.path.join(ruta, d))]
            )
        except FileNotFoundError:
            self.total_samples = 0
            messagebox.showerror("Error", f"El directorio no existe: {ruta}")
            return

        self.processed_samples = 0

        # Crear ventana de progreso
        ventana_progreso = Toplevel(self.ventana)
        ventana_progreso.title("Proceso de Ensamblaje...")
        ventana_progreso.geometry("700x450")

        tk.Label(
            ventana_progreso,
            text="Ejecutando script de ensamblaje...",
            font=("Helvetica", 12, "bold"),
        ).pack(pady=10)

        self.progress_label_var.set(f"0/{self.total_samples} muestras ensambladas")
        tk.Label(
            ventana_progreso,
            textvariable=self.progress_label_var,
            font=("Helvetica", 10),
        ).pack(pady=5)

        barra_progreso = ttk.Progressbar(
            ventana_progreso,
            orient="horizontal",
            mode="determinate",
            length=600,
            variable=self.progress_var,
        )
        barra_progreso.pack(pady=10)
        self.progress_var.set(0)

        area_texto = Text(
            ventana_progreso, wrap="word", height=20, width=80, state="disabled"
        )
        area_texto.pack(side="left", fill="both", expand=True, padx=10)

        barra_desplazamiento = Scrollbar(ventana_progreso, command=area_texto.yview)
        barra_desplazamiento.pack(side="right", fill="y", padx=(0, 10))
        area_texto["yscrollcommand"] = barra_desplazamiento.set

        def ejecutar_script_thread():
            try:

                def update_gui(
                    line=None, progress=None, progress_text=None, final_message=None
                ):
                    if line:
                        self.add_to_log(area_texto, line)
                    if progress is not None:
                        self.progress_var.set(progress)
                    if progress_text:
                        self.progress_label_var.set(progress_text)
                    if final_message:
                        success, msg = final_message
                        if success:
                            messagebox.showinfo("Éxito", msg)
                        else:
                            messagebox.showerror("Error", msg)
                        ventana_progreso.destroy()

                proceso = subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                    encoding="utf-8",
                )

                # Patrones para detectar el progreso
                sample_start_pattern_old = re.compile(r"Procesando archivo\s+:\s+\S+")
                sample_start_pattern_new = re.compile(r"Procesando muestra:\s+\S+")

                for line in iter(proceso.stdout.readline, ""):
                    self.ventana.after(0, update_gui, line)

                    if sample_start_pattern_old.search(
                        line
                    ) or sample_start_pattern_new.search(line):
                        self.processed_samples += 1
                        progress = (
                            (self.processed_samples / self.total_samples) * 100
                            if self.total_samples > 0
                            else 0
                        )
                        progress_text = f"{self.processed_samples}/{self.total_samples} muestras ensambladas"
                        self.ventana.after(0, update_gui, None, progress, progress_text)

                proceso.wait()

                if proceso.returncode == 0:
                    self.ventana.after(
                        0,
                        update_gui,
                        None,
                        None,
                        None,
                        (True, "¡Ensamblaje ejecutado con éxito!"),
                    )
                else:
                    self.ventana.after(
                        0,
                        update_gui,
                        None,
                        None,
                        None,
                        (
                            False,
                            f"Error durante la ejecución del ensamblaje. Código de salida: {proceso.returncode}",
                        ),
                    )
            except Exception as e:

                def show_error():
                    messagebox.showerror("Error", f"Error: {e}")
                    ventana_progreso.destroy()

                self.ventana.after(0, show_error)

        Thread(target=ejecutar_script_thread).start()

    def ejecutar_de_novo_ensamblaje(self):
        ruta = self.ruta_var.get()
        sequence_type_value = self.denovo_sequence_type_var.get()

        if not ruta or not sequence_type_value:
            messagebox.showerror(
                "Error", "Seleccione el directorio de datos y el tipo de secuenciación."
            )
            return

        try:
            self.total_samples = len(
                [d for d in os.listdir(ruta) if os.path.isdir(os.path.join(ruta, d))]
            )
        except FileNotFoundError:
            self.total_samples = 0
            messagebox.showerror("Error", f"El directorio no existe: {ruta}")
            return

        self.processed_samples = 0

        ventana_progreso = Toplevel(self.ventana)
        ventana_progreso.title("Proceso de Ensamblaje De Novo...")
        ventana_progreso.geometry("700x450")

        tk.Label(
            ventana_progreso,
            text="Ejecutando script de ensamblaje De Novo...",
            font=("Helvetica", 12, "bold"),
        ).pack(pady=10)

        self.progress_label_var.set(f"0/{self.total_samples} muestras ensambladas")
        tk.Label(
            ventana_progreso,
            textvariable=self.progress_label_var,
            font=("Helvetica", 10),
        ).pack(pady=5)

        barra_progreso = ttk.Progressbar(
            ventana_progreso,
            orient="horizontal",
            mode="determinate",
            length=600,
            variable=self.progress_var,
        )
        barra_progreso.pack(pady=10)

        self.progress_var.set(0)

        area_texto = Text(
            ventana_progreso, wrap="word", height=20, width=80, state="disabled"
        )
        area_texto.pack(side="left", fill="both", expand=True, padx=10)

        barra_desplazamiento = Scrollbar(ventana_progreso, command=area_texto.yview)
        barra_desplazamiento.pack(side="right", fill="y", padx=(0, 10))
        area_texto["yscrollcommand"] = barra_desplazamiento.set

        def ejecutar_script_thread():
            try:
                current_dir = os.path.dirname(os.path.abspath(__file__))
                script_path = os.path.join(current_dir, "De_Novo_Assembly.sh")
                if not os.path.exists(script_path):
                    self.ventana.after(
                        0,
                        lambda: messagebox.showerror(
                            "Error", f"El script no se encuentra: De_Novo_Assembly.sh"
                        ),
                    )
                    return

                def update_gui(
                    line=None, progress=None, progress_text=None, final_message=None
                ):
                    if line:
                        self.add_to_log(area_texto, line)
                    if progress is not None:
                        self.progress_var.set(progress)
                    if progress_text:
                        self.progress_label_var.set(progress_text)
                    if final_message:
                        success, msg = final_message
                        if success:
                            messagebox.showinfo("Éxito", msg)
                        else:
                            messagebox.showerror("Error", msg)
                        ventana_progreso.destroy()

                command = ["bash", script_path, ruta, sequence_type_value]
                proceso = subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                    encoding="utf-8",
                )

                sample_start_pattern = re.compile(r"Processing sample: (.*)")

                for line in iter(proceso.stdout.readline, ""):
                    self.ventana.after(0, update_gui, line)

                    match = sample_start_pattern.search(line)
                    if match:
                        self.processed_samples += 1
                        progress = 0
                        if self.total_samples > 0:
                            progress = (
                                self.processed_samples / self.total_samples
                            ) * 100
                        progress_text = f"{self.processed_samples}/{self.total_samples} muestras ensambladas"
                        self.ventana.after(0, update_gui, None, progress, progress_text)

                proceso.wait()

                if proceso.returncode == 0:
                    self.ventana.after(
                        0,
                        update_gui,
                        None,
                        None,
                        None,
                        (True, "¡Ensamblaje De Novo ejecutado con éxito!"),
                    )
                else:
                    self.ventana.after(
                        0,
                        update_gui,
                        None,
                        None,
                        None,
                        (
                            False,
                            f"Error durante la ejecución del ensamblaje De Novo. Código de salida: {proceso.returncode}",
                        ),
                    )
            except FileNotFoundError:

                def show_error():
                    messagebox.showerror(
                        "Error",
                        "Asegúrese de que el script de ensamblaje De Novo exista y sea ejecutable.",
                    )
                    ventana_progreso.destroy()

                self.ventana.after(0, show_error)
            except Exception as e:

                def show_error():
                    messagebox.showerror("Error", f"Error: {e}")
                    ventana_progreso.destroy()

                self.ventana.after(0, show_error)

        Thread(target=ejecutar_script_thread).start()

    def ejecutar_llamado_variantes(self):
        # Implementación de la función de llamado de variantes
        messagebox.showinfo("Llamado de Variantes", "Función en desarrollo...")

    def ejecutar_viral_branch(self):
        # Implementación de la función de Viral-Branch
        messagebox.showinfo("Viral-Branch", "Función en desarrollo...")

    def ejecutar_organizar_fastq(self):
        directorio = self.carpeta_fastq_var.get()
        if not directorio:
            messagebox.showerror("Error", "Seleccione un directorio a organizar.")
            return

        ventana_progreso = Toplevel(self.ventana)
        ventana_progreso.title("Organizando FASTQ...")
        ventana_progreso.geometry("700x400")

        tk.Label(
            ventana_progreso,
            text="Organizando archivos FASTQ...",
            font=("Helvetica", 12, "bold"),
        ).pack(pady=10)

        area_texto = Text(
            ventana_progreso, wrap="word", height=20, width=80, state="disabled"
        )
        area_texto.pack(side="left", fill="both", expand=True, padx=10)

        barra_desplazamiento = Scrollbar(ventana_progreso, command=area_texto.yview)
        barra_desplazamiento.pack(side="right", fill="y", padx=(0, 10))
        area_texto["yscrollcommand"] = barra_desplazamiento.set

        def ejecutar_thread():
            try:
                # Redirigir stdout para capturar la salida
                stdout_original = sys.stdout
                salida = io.StringIO()
                sys.stdout = salida

                organizar_fastq(directorio)

                # Restaurar stdout
                sys.stdout = stdout_original
                output = salida.getvalue()

                # Schedule GUI updates in the main thread
                def update_gui():
                    self.add_to_log(area_texto, output)
                    messagebox.showinfo("Éxito", "Organización completada.")
                    ventana_progreso.destroy()

                self.ventana.after(0, update_gui)

            except Exception as e:
                # Schedule error message in the main thread
                def show_error():
                    messagebox.showerror("Error", f"Error durante la organización: {e}")
                    ventana_progreso.destroy()

                self.ventana.after(0, show_error)

        Thread(target=ejecutar_thread).start()

    # Añadir estas funciones para manejar el sample sheet
    def buscar_sample_sheet(self):
        archivo = filedialog.askopenfilename(
            title="Seleccionar Sample Sheet",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        self.sample_sheet_var.set(archivo)

    def buscar_archivo_consenso(self):
        archivo = filedialog.askopenfilename(
            title="Seleccionar archivo de consenso",
            filetypes=[("FASTA files", "*.fasta *.fa *.fas"), ("All files", "*.*")],
        )
        self.consensus_file_var.set(archivo)

    def buscar_dir_salida(self):
        directorio = filedialog.askdirectory(title="Seleccionar directorio de salida")
        self.output_dir_var.set(directorio)

    def procesar_sample_sheet(self):
        sample_sheet = self.sample_sheet_var.get()
        consensus_file = self.consensus_file_var.get()
        output_dir = self.output_dir_var.get()

        if not sample_sheet or not consensus_file:
            messagebox.showerror(
                "Error", "Seleccione un sample sheet y un archivo de consenso."
            )
            return

        # Crear ventana de progreso
        ventana_progreso = Toplevel(self.ventana)
        ventana_progreso.title("Procesando Sample Sheet...")
        ventana_progreso.geometry("700x400")

        tk.Label(
            ventana_progreso,
            text="Procesando Sample Sheet y renombrando archivos...",
            font=("Helvetica", 12, "bold"),
        ).pack(pady=10)

        area_texto = Text(
            ventana_progreso, wrap="word", height=20, width=80, state="disabled"
        )
        area_texto.pack(side="left", fill="both", expand=True, padx=10)

        barra_desplazamiento = Scrollbar(ventana_progreso, command=area_texto.yview)
        barra_desplazamiento.pack(side="right", fill="y", padx=(0, 10))
        area_texto["yscrollcommand"] = barra_desplazamiento.set

        def ejecutar_thread():
            try:
                # Redirigir stdout a un StringIO para capturar la salida
                stdout_original = sys.stdout
                salida = io.StringIO()
                sys.stdout = salida

                # Importar el módulo de sample sheet
                sys.path.append(
                    os.path.join(
                        os.path.dirname(os.path.abspath(__file__)), "..", "Scripts"
                    )
                )
                from sample_sheet import SampleSheet

                # Procesar sample sheet
                sheet = SampleSheet(sample_sheet)
                archivos = sheet.rename_consensus_files(consensus_file, output_dir)

                # Restaurar stdout
                sys.stdout = stdout_original
                output = salida.getvalue()

                def update_gui():
                    self.add_to_log(area_texto, output)
                    if archivos:
                        messagebox.showinfo(
                            "Éxito",
                            f"Procesamiento completado.\n{len(archivos)} archivos generados.",
                        )
                    else:
                        messagebox.showwarning(
                            "Advertencia", "No se generaron archivos."
                        )
                    ventana_progreso.destroy()

                self.ventana.after(0, update_gui)

            except ImportError:

                def show_error():
                    messagebox.showerror(
                        "Error",
                        "No se pudo encontrar el módulo 'sample_sheet.py'. Asegúrese de que esté en la ruta correcta.",
                    )
                    ventana_progreso.destroy()

                self.ventana.after(0, show_error)
            except Exception as e:

                def show_error():
                    messagebox.showerror(
                        "Error", f"Error durante el procesamiento: {e}"
                    )
                    ventana_progreso.destroy()

                self.ventana.after(0, show_error)

        Thread(target=ejecutar_thread).start()

    def salir_app(self):
        self.ventana.quit()

    def acerca_de_consenso_d(self):
        messagebox.showinfo(
            "Acerca de CONSENSO_D",
            "CONSENSO_D Viral Genome Assembler\n"
            + "Versión 1.0\n\n"
            + "Funcionalidades:\n"
            + "• Ensamblaje de genomas virales\n"
            + "• Llamado de variantes\n"
            + "• Determinación de linajes (Viral-Branch)\n"
            + "• Organización de archivos FASTQ\n"
            + "• Procesamiento de Sample Sheet\n\n"
            + "Instituto de Ciencias Sostenibles",
        )

    def run(self):
        self.ventana.mainloop()


# Ejecutar la aplicación
if __name__ == "__main__":
    app = ConsensoGUI()
    app.run()
