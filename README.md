# CONSENSO_D
CONSENSO_D es un conjunto de scripts diseñado por el equipo de secuenciación y genómica del Instituto de Ciencias Sostenibles de Nicaragua con la finalidad de mapear, ensamblar y generar secuencias consenso de los 4 serotipos de Dengue y del virus de la Rabia a partir de archivos fastq de secuenciadores NANOPORE.

CONSENSO_D usa multiples Scripts python y Bash que permiten realizar diferentes tareas:

  - Interfaz grafica intuitiva para manejar todos los procesos
  - Ensamblaje de Genomas del los 4 serotipos del dengue, Virus de la Rabia y el SARS CoV 2
  - Permite realizar ensamblaje de plataformas Nanopore e Illumina
  - Realiza control de calidad de las secuencias mediante fastWC y Multi QC
  - Realiza determinacion de linajes del virus del dengue bajo la nomenclatura Dengue-Lineages, implementando un modelo Random Forest de Aprendizaje automatico para clasificar los linajes (https://github.com/gvaleman/Viral-Branch)
  - Determinación de mutaciones para los serotipos del dengue

Configurando CONSENSO

  CLONAR EL REPOSITORIO

En la terminal ejecutar el siguientes comando
```
cd
git clone https://github.com/gvaleman/CONSENSO_D.git
```
para la configuración, es necesario tener instalado Miniconda3. Es necesario activar un ambiente conda de la siguiente manera
```
cd # sale a la carpeta home
cd CONSENSO_D # entrar a la carpeta de CONSENSO
conda env create -f environment.yml
```
Despues, solo hay que esperar a que se realice la instalación del entorno para iniciar a utilizar CONSENSO.

# USANDO CONSENSO
INSTRUCCIONES:
- Acceder a la carpeta CONSENSO_D y entrar a la carpeta bin.
- Abrir una terminal y ejecutar el Script CONSENSO.py. esto se puede hacer de la siguiente manera abriendo una terminal
  ```
  cd # sale a la carpeta home
  cd CONSENSO_D # entrar a la carpeta de CONSENSO
  cd bin #entrar a la carpeta bin
  python3 CONSENSO.py -- gui #Ejecutar el Script
  ```
- Esto abrira una interfaz gráfica de usuario.
- Ensamblaje de genomas: Seleccionar "Ensamblaje de genomas a la Izquierda
- En la pantalla principal seleccionar la tecnologia de secuenciación, serotipo o virus, esquema de primers (si aplica), y dar click en Ejecutar
- El Ensamblaje de genomas se ejecutara, mostrando el progreso. al finalizr un archivo llamado "all_consnsus.fasta" se guardara en ell directorio indicado

Ejemplo:
```
./CONSENSO DENV_1 '/home/ics2/CONSENSO_D/DENV1_fastq_pass'
```
  donde:
  
   ./CONSENSO : Nombre del script
   
   DENV_1 : nombre del serotipo
   
   '/home/ics2/CONSENSO_D/DENV1_fastq_pass' : Dirección exacta de la carpeta "DENV1_fastq_pass" que contiene los archivos fastq de dengue 1


*---------------------------------* 


*Acerca del script contarX*

contarX es un script que cuenta el número de gaps (X) en una o más secuencias consenso.
Obtener el número de gaps es importante porque brinda una idea de la covertura de la secuencia de referencia obtenida con respecto a una secuencia completa de referencia. Esto brinda una idea de la calidad de la secuencia y de manera indirecta brinda nociones sobre la calidad de la secuenciación, calidad de la muestra, calidad de los primers y/o reactivos.

contarX toma cada una de la secuencias de nucleótidos y las trata como una cadena de texto, evaluando caracter por caracter. Todos los caracteres "X" los cuenta por cada secuencia y finalmente muestra en la terminal el nombre de la secuencia y el número de X

INSTRUCCIONES
- Abrir una terminal, haciendo click sobre el icono de la terminal o con la combinacion de teclas CTRL+ALT+T
- escrinir el nombre del script: ./contarX
- Escribir la dirección del archivo fasta o arrastrarlo hacia la terminal.
- Presionar la tecla enter

Ejemplo
```
./contarX '/home/ics2/CONSENSO_D/ALL_CONSENSUS.fasta'
```
      
 Donde
 
   ./contarX: nombre del script
   
   '/home/ics2/CONSENSO_D/ALL_CONSENSUS.fasta': dirección exacta del archivo fasta 
   
   
   *---------------------------------* 


*Acerca del script RENOMBRAR*

RENOMBRAR es un script que permite renombrar una o más secuencias consensos de manera autómatica en un archivo multifasta.
Se requeren dos elementos:
 1- Una archivo fasta con las secuencias sin renombrar
 2- Un archivo txt con el nombre antiguo y nuevo de las secuencias. Este archivo txt se puede crear a partor de una hoja de excel, copiando y pegando el contenido de la hoja de excel a un editor de texto. es importante procurar que no quede ninguna linea vacia al final del texto. A continuación se muestra un ejemplo de la estructura final del archivo txt
 ```
barcode01 nuevo_nombre_1
barcode02	nuevo_nombre_2
barcode03	nuevo_nombre_3
barcode04	nuevo_nombre_4
```

INSTRUCCIONES
- Abrir una terminal, haciendo click sobre el icono de la terminal o con la combinacion de teclas CTRL+ALT+T
- Escrinir el nombre del script: ./RENOMBRAR
- Escribir la dirección del de texto con los viejos y nuevos nombres o arrastrarlo hacia la terminal.
- seguido de un espacio, escribir la dirección del archivo fasta sin renombrar o arrastrarlo hacia la terminal.
- Escribir el simbolo ">" y la ubicación exacta y el nombre donde es guardará el nuevo archivo fasta

Ejemplo
```
./RENOMBRAR '/home/ics2/mi_carpeta/nombres.txt'  '/home/ics2/mi_carpeta/viejos_nombres.fasta' > /home/ics2/mi_carpeta/nuevo.fasta'
```

  donde
  
 ./RENOMBRAR: Nombre del script
 
 '/home/ics2/mi_carpeta/nombres.txt': Ubicación exacta del archivo de texto que contiene los nombres antiguos y nuevos
 
 '/home/ics2/mi_carpeta/viejos_nombres.fasta': Ubicación exacta del archivo fasta que contiene las secuencias con los viejos nombres
 
  > /home/ics2/mi_carpeta/nuevo.fasta': Ubicación excta y nommbre del archivo que contendrá los archivos fasta renombrados. En este caso, el archivo se guardará en la la carpeta home, ics2, mi_carpeta, con el nombre "nuevo.fasta". Es importante guardarlo siempre con la extención ".fasta"
 
