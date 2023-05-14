# CONSENSO_D
CONSENSO_D es un conjunto de 3 scripts diseñado por el equipo de secuenciación y genómica del instituto de ciencias sostenibles con la finalidad de mapear, ensamblar y generar una secuencia consenso de los 4 serotipos de Dengue a partir de archivos fastq de secuenciadores NANOPORE.
CONSENSO automatiza el paso a paso del proceso para la generación de una secuencia consenso.

Configurando los scripts

CLONAR EL REPOSITORIO
En la terminal ejecutar el siguientes comando
```
git clone https://github.com/gvaleman/CONSENSO_D.git
```
para la configuración de los scripts, escribir en la terminal
```
cd CONSENSO_D
chmod +x CONFIGURAR
./CONFIGURAR
```
*Acerca del script CONENSO*

CONSENSO es un script diseñado por el equipo de secuenciación y genómica del ICS con la finalidad de mapear, ensamblar y generar una secuencia consenso de los 4 serotipos de Dengue a partir de archivos fastq de secuenciadores NANOPORE.
CONSENSO automatiza el paso a paso del proceso para la generación de una secuencia consenso.

EL proceso paso a paso para la obtención de una secuencia consenso incluye por cada secuencia
   1- Descompreción: LOs archivos fastq productos de la secuenciación nanopore están comprimidos en un archivo gz. El objetivo de almacenar los archivos en este formato es no utilizar mucho espacio en el disco de la PC. Los programas no utilizan estos formatos por lo que es necesario descomprimirlos
   2- Concatenación: Cada muestra secuenciada genera varios archivos fastq comprimido. Cada archivo contiene entre 1000 y 4000 reads según la configuración de la secuenciación. Todos los archivos fastq se deben concatenar para dar paso a un solo archivo que contenga todas las lecturas.
   3- Mapeo contra referencia: Las lecturas obtenidas de la secuenciación (que ahora están en el archivo fastq concatenado) se deben alinear contra una secuencia de referencia. En este paso se realiza la alineación y se genera un archivo que contiene la información del alineamiento de todos los segmentos (reads) a una secuencia de referencia. Este archivo se genera en formato SAM
   4- Convertir de SAM a BAM: El archivo SAM es convertido a un archivo de formato BAM que es una versión binaria comprimida del archivo SAM. que es útil porque es mucho más compacto que el archivo SAM y se puede indexar para un acceso rápido a regiones específicas, lo que lo hace ideal para el análisis de grandes cantidades de datos de secuenciación.
   5- Ordenar y filtrar las lecturas alineadas: Se realiza esta filtración y ordenamiento para mejorar la calidad de los datos y aumentar la precisión de la secuencia consenso resultante.
   6-Indexación del mapeo: Una vez teniendo la información de la alineación de los fragmentos filtrado y ordenado, se genera un archivo de indexación. Es decir, un archivo que contiene la información y las coordenadas de los resultados y cambios de las bases nitrogenadas de los segmentos alineados con respecto a la secuencia de referencia. 
   7- Obtner un archivo VCF:el archivo BCF se genera a partir del archivo indexado. Este archivo contiene información específican (cambio y posición del cambio) de las bases nitrogenadas del ADN/ARN secuenciado. 
   8- Secuencia consenso: En este paso, se toma la información del archivo VCF que contiene los cambios del material secuenciado e inserta estos cambios a una secuencia de referencia. 
   9- Renombrar secuencia de referencia: UNa vez que se insertaron los cambios a la secuencia se debe referencia se debe renombrar esta secuencia con el nombre de la muestra o barcode analizado.
   10- fastq a fasta: La secuencia definitiva es almacenada en formato fastq, por lo que se debe convertir a formato fasta

Estos 10 pasos son realizados muestra por muestra y de manera automática por el scipt CONSENSO. Adicionalmente, una vez generadas todas las secuencias consenso, el script genera un solo archivo "ALL_CONSENSUS.FASTA" que contiene a todas las secuencias

INSTRUCCIONES:
- seleccionar las carpetas de los barcode de serotipo que desea ensamblar y guardarlas en un directorio de trabajo
- Abrir una terminal, haciendo click sobre el icono de la terminal o con la combinacion de teclas CTRL+ALT+T
- En la terminal escribir ./CONSENSO: Esto hace el llamado al scrip ejecutable
- Escribir el serotipo que se desea ensamblar: El nombre del virus seguido de un guión bajo, seguido del serotipo, sin espacios (ejemplo: DENV_1, DENV_2, DENV_3 o DENV_4 )
- Escribir o arrastrar el directorio de trabajo que contiene las carpetas con los archivos fastq
- presionar ENTER y esperar que el proceso de ensamblaje termine

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
 
