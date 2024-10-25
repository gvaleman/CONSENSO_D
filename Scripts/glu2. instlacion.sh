#!/bin/bash

# Variables
ENV_NAME="CONSENSO_D"
GLUE_HOME="$HOME/CONSENSO_D/DENV_GLU/gluetools"
BLAST_DIR="$HOME/CONSENSO_D/blast-2.2.31-pl526he19e7b1_5"
MYSQL_USER="gluetools"
MYSQL_PASSWORD="glue12345"
MYSQL_DATABASE="GLUE_TOOLS"
SQL_FILE_PATH="$HOME/CONSENSO_D/dengue_glue.sql" # Asegúrate de ajustar el path a donde está tu archivo dengue_glue.sql

# Inicializar conda
eval "$(conda shell.bash hook)"

# Verificar si el entorno conda ya existe
if conda env list | grep -q $ENV_NAME; then
    echo "El entorno conda '$ENV_NAME' ya existe. Activando..."
else
    echo "El entorno conda '$ENV_NAME' no existe. Creando..."
    conda env create -f environment.yml
fi

# Activar el entorno conda
conda activate $ENV_NAME

# Obtener las rutas de MAFFT y RAxML
MAFFT_DIR=$(which mafft)
RAXML_DIR=$(which raxmlHPC)

# Instalar MySQL/MariaDB si no está instalado
if ! command -v mariadb &> /dev/null; then
    echo "MariaDB no está instalado. Instalando..."
    conda install -c conda-forge mariadb
else
    echo "MariaDB ya está instalado."
fi

# Verificar si el servicio de MariaDB ya está iniciado
if ! pgrep mysqld &> /dev/null; then
    echo "Iniciando el servicio de MariaDB..."
    mariadb-install-db --user="$USER" --basedir="$(conda info --base)/envs/$ENV_NAME" --datadir="$(conda info --base)/envs/$ENV_NAME/var/lib/mysql"
    mysqld_safe --datadir="$(conda info --base)/envs/$ENV_NAME/var/lib/mysql" &
    # Esperar a que MariaDB se inicie
    sleep 10
else
    echo "El servicio de MariaDB ya está en ejecución."
fi

# Configurar MariaDB si no lo está
if ! mariadb -u $MYSQL_USER -p$MYSQL_PASSWORD -e "USE $MYSQL_DATABASE;" &> /dev/null; then
    echo "Configurando MariaDB..."
    mariadb -u root <<EOF
CREATE USER IF NOT EXISTS '$MYSQL_USER'@'localhost' IDENTIFIED BY '$MYSQL_PASSWORD';
CREATE DATABASE IF NOT EXISTS $MYSQL_DATABASE CHARACTER SET UTF8;
GRANT ALL PRIVILEGES ON $MYSQL_DATABASE.* TO '$MYSQL_USER'@'localhost';
FLUSH PRIVILEGES;
EOF
else
    echo "MariaDB ya está configurado con el usuario y base de datos especificados."
fi

# Verificar conexión a MariaDB
mariadb -u $MYSQL_USER -p$MYSQL_PASSWORD -e "USE $MYSQL_DATABASE; SHOW TABLES;"

# Inicializar la base de datos con el archivo SQL si no se ha hecho
if ! mariadb -u $MYSQL_USER -p$MYSQL_PASSWORD -e "USE $MYSQL_DATABASE; SHOW TABLES;" | grep -q 'some_table_name'; then
    echo "Inicializando la base de datos con dengue_glue.sql..."
    mariadb -u $MYSQL_USER -p$MYSQL_PASSWORD $MYSQL_DATABASE < $SQL_FILE_PATH
else
    echo "La base de datos ya está inicializada."
fi

# Configuración de GLUE
export GLUE_HOME=$GLUE_HOME
export PATH=${PATH}:${GLUE_HOME}/bin

# Guardar las variables de entorno en ~/.bashrc para futuras sesiones si no están ya presentes
if ! grep -q 'export GLUE_HOME=' ~/.bashrc; then
    echo 'export GLUE_HOME=$HOME/CONSENSO_D/DENV_GLU/gluetools' >> ~/.bashrc
fi

if ! grep -q 'export PATH=' ~/.bashrc | grep '${GLUE_HOME}/bin'; then
    echo 'export PATH=${PATH}:${GLUE_HOME}/bin' >> ~/.bashrc
fi

# Hacer el script ejecutable si no lo es
if [ ! -x "${GLUE_HOME}/bin/gluetools.sh" ]; then
    chmod u+x ${GLUE_HOME}/bin/gluetools.sh
fi

# Configuración del archivo XML si no existe
if [ ! -f "${GLUE_HOME}/conf/gluetools-config.xml" ]; then
    cat <<EOL > ${GLUE_HOME}/conf/gluetools-config.xml
<gluetools>
    <database>
        <username>$MYSQL_USER</username>
        <password>$MYSQL_PASSWORD</password>
        <vendor>MySQL</vendor>
        <jdbcUrl>jdbc:mysql://localhost:3306/$MYSQL_DATABASE?characterEncoding=UTF-8</jdbcUrl>
    </database>
    <properties>
        <property>
            <name>gluetools.core.programs.blast.blastn.executable</name>
            <value>${BLAST_DIR}/bin/blastn</value>
        </property>
        <property>
            <name>gluetools.core.programs.blast.tblastn.executable</name>
            <value>${BLAST_DIR}/bin/tblastn</value>
        </property>
        <property>
            <name>gluetools.core.programs.blast.makeblastdb.executable</name>
            <value>${BLAST_DIR}/bin/makeblastdb</value>
        </property>
        <property>
            <name>gluetools.core.programs.blast.temp.dir</name>
            <value>${GLUE_HOME}/tmp/blastfiles</value>
        </property>
        <property>
            <name>gluetools.core.programs.blast.db.dir</name>
            <value>${GLUE_HOME}/tmp/blastdbs</value>
        </property>
        <property>
            <name>gluetools.core.programs.blast.search.threads</name>
            <value>4</value>
        </property>
        <!-- MAFFT-specific config -->
        <property>
            <name>gluetools.core.programs.mafft.executable</name>
            <value>${MAFFT_DIR}</value>
        </property>
        <property>
            <name>gluetools.core.programs.mafft.cpus</name>
            <value>4</value>
        </property>
        <property>
            <name>gluetools.core.programs.mafft.temp.dir</name>
            <value>${GLUE_HOME}/tmp/mafftfiles</value>
        </property>
        <!-- RAxML-specific config -->
        <property>
            <name>gluetools.core.programs.raxml.raxmlhpc.executable</name>
            <value>${RAXML_DIR}</value>
        </property>
        <property>
            <name>gluetools.core.programs.raxml.raxmlhpc.cpus</name>
            <value>4</value>
        </property>
        <property>
            <name>gluetools.core.programs.raxml.temp.dir</name>
            <value>${GLUE_HOME}/tmp/raxmlfiles</value>
        </property>
    </properties>
</gluetools>
EOL
fi

# Crear directorios temporales si no existen
mkdir -p ${GLUE_HOME}/tmp/blastfiles
mkdir -p ${GLUE_HOME}/tmp/blastdbs
mkdir -p ${GLUE_HOME}/tmp/mafftfiles
mkdir -p ${GLUE_HOME}/tmp/raxmlfiles

# Mensaje final
echo "Instalación de GLUE completada con éxito. Ejecuta '${GLUE_HOME}/bin/gluetools.sh' para iniciar GLUE."

# Ejecutar GLUE para verificar que funciona
${GLUE_HOME}/bin/gluetools.sh
