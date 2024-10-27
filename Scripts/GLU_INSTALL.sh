#!/bin/bash

# Variables
ENV_NAME="CONSENSO_D"
GLUE_HOME="$HOME/CONSENSO_D/DENV_GLU/gluetools"
BLAST_DIR="$HOME/CONSENSO_D/blast-2.2.31-pl526he19e7b1_5"
MYSQL_USER="gluetools"
MYSQL_PASSWORD="glue12345"
MYSQL_DATABASE="GLUE_TOOLS"
SQL_FILE_PATH="$HOME/CONSENSO_D/dengue_glue.sql"

# Verificar Java
if ! command -v java &> /dev/null; then
    echo "Error: Java no está instalado. Por favor, instale Java 1.8.0 o superior."
    exit 1
fi

java_version=$(java -version 2>&1 | head -n 1 | cut -d'"' -f2 | cut -d'.' -f2)
if [ "$java_version" -lt "8" ]; then
    echo "Error: Se requiere Java 1.8.0 o superior."
    exit 1
fi

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

# Verificar BLAST+
if [ ! -f "${BLAST_DIR}/bin/blastn" ] || [ ! -f "${BLAST_DIR}/bin/tblastn" ] || [ ! -f "${BLAST_DIR}/bin/makeblastdb" ]; then
    echo "Error: No se encontraron los ejecutables de BLAST+ en ${BLAST_DIR}/bin"
    echo "Por favor, instale BLAST+ 2.2.31 desde la página de NCBI"
    exit 1
fi

# Verificar MAFFT y RAxML
MAFFT_DIR=$(which mafft)
RAXML_DIR=$(which raxmlHPC)

if [ -z "$MAFFT_DIR" ]; then
    echo "Error: MAFFT no está instalado"
    exit 1
fi

if [ -z "$RAXML_DIR" ]; then
    echo "Error: RAxML no está instalado"
    exit 1
fi

# Verificar/Instalar MariaDB
if ! command -v mariadb &> /dev/null; then
    echo "MariaDB no está instalado. Instalando..."
    conda install -c conda-forge mariadb
fi

# Iniciar MariaDB si no está corriendo
if ! pgrep mysqld &> /dev/null; then
    echo "Iniciando el servicio de MariaDB..."
    # Verificar si el directorio de datos existe
    DB_DIR="$(conda info --base)/envs/$ENV_NAME/var/lib/mysql"
    if [ ! -d "$DB_DIR" ]; then
        mariadb-install-db --user="$USER" --basedir="$(conda info --base)/envs/$ENV_NAME" --datadir="$DB_DIR"
    fi
    mysqld_safe --datadir="$DB_DIR" &
    sleep 10
fi

# Verificar/Instalar MariaDB en el entorno conda
if ! conda list | grep -q "mariadb"; then
    echo "Instalando MariaDB en el entorno conda $ENV_NAME..."
    conda install -c conda-forge mariadb-connector-c mariadb -y
    if [ $? -ne 0 ]; then
        echo "Error: No se pudo instalar MariaDB"
        exit 1
    fi
fi

# Asegurarse de que los binarios de MariaDB estén en el PATH
export PATH="$(conda info --base)/envs/$ENV_NAME/bin:$PATH"

# Definir el directorio de datos
DB_DIR="$(conda info --base)/envs/$ENV_NAME/var/lib/mysql"
mkdir -p "$DB_DIR"

# Verificar si el servicio de MariaDB está corriendo
if ! pgrep mysqld &> /dev/null; then
    echo "Iniciando el servicio de MariaDB..."
    
    # Inicializar la base de datos si no está inicializada
    if [ ! -f "$DB_DIR/ibdata1" ]; then
        echo "Inicializando la base de datos MariaDB..."
        mysql_install_db --datadir="$DB_DIR" --auth-root-authentication-method=normal
        if [ $? -ne 0 ]; then
            echo "Error: No se pudo inicializar la base de datos"
            exit 1
        fi
    fi

    # Iniciar el servidor MariaDB
    mysqld_safe --datadir="$DB_DIR" --skip-grant-tables &
    
    # Esperar a que el servidor esté disponible
    echo "Esperando a que el servidor MariaDB esté disponible..."
    sleep 10
    
    # Asegurar la instalación
    mysql -u root <<EOF
FLUSH PRIVILEGES;
ALTER USER 'root'@'localhost' IDENTIFIED BY 'root';
CREATE USER IF NOT EXISTS '$MYSQL_USER'@'localhost' IDENTIFIED BY '$MYSQL_PASSWORD';
CREATE DATABASE IF NOT EXISTS $MYSQL_DATABASE CHARACTER SET UTF8;
GRANT ALL PRIVILEGES ON $MYSQL_DATABASE.* TO '$MYSQL_USER'@'localhost';
FLUSH PRIVILEGES;
EOF

    # Reiniciar el servidor con la configuración normal
    pkill mysqld
    sleep 5
    mysqld_safe --datadir="$DB_DIR" &
    sleep 10
fi

# Verificar la conexión
echo "Verificando la conexión a MariaDB..."
if mysql -u $MYSQL_USER -p$MYSQL_PASSWORD -e "USE $MYSQL_DATABASE;" &> /dev/null; then
    echo "Conexión a MariaDB establecida correctamente"
else
    echo "Error: No se pudo conectar a MariaDB con las credenciales proporcionadas"
    exit 1
fi

# Configurar GLUE_HOME
echo "Configurando GLUE_HOME..."
if [ ! -d "$GLUE_HOME" ]; then
    echo "Error: No se encontró el directorio GLUE_HOME en $GLUE_HOME"
    exit 1
fi

# Actualizar .bashrc
{
    echo "export GLUE_HOME=$GLUE_HOME"
    echo 'export PATH=${PATH}:${GLUE_HOME}/bin'
} >> ~/.bashrc

# Hacer ejecutable gluetools.sh
chmod u+x "${GLUE_HOME}/bin/gluetools.sh"

# Crear y configurar gluetools-config.xml
CONFIG_FILE="${GLUE_HOME}/conf/gluetools-config.xml"
mkdir -p "${GLUE_HOME}/conf"

cat > "$CONFIG_FILE" <<EOL
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

# Crear directorios temporales
for DIR in blastfiles blastdbs mafftfiles raxmlfiles; do
    mkdir -p "${GLUE_HOME}/tmp/$DIR"
done

echo "Instalación completada. Verificando GLUE..."
source ~/.bashrc
${GLUE_HOME}/bin/gluetools.sh
