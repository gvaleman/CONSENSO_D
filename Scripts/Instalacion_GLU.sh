#!/bin/bash

# Variables
ENV_NAME="CONSENSO_D"
GLUE_HOME="$HOME/CONSENSO_D/DENV_GLU/gluetools"
BLAST_DIR="$HOME/CONSENSO_D/blast-2.2.31-pl526he19e7b1_5"
MYSQL_USER="gluetools"
MYSQL_PASSWORD="glue12345"
MYSQL_DATABASE="GLUE_TOOLS"

# Inicializar conda
eval "$(conda shell.bash hook)"

# Crear el entorno conda si no existe
if ! conda env list | grep -q $ENV_NAME; then
    conda env create -f environment.yml
fi

# Activar el entorno conda
conda activate $ENV_NAME

# Obtener las rutas de MAFFT y RAxML
MAFFT_DIR=$(which mafft)
RAXML_DIR=$(which raxmlHPC)

# Instalar MySQL/MariaDB si no está instalado
if ! command -v mariadb &> /dev/null; then
    conda install -c conda-forge mariadb
fi

# Iniciar el servicio de MariaDB
mariadb-install-db --user="$USER" --basedir="$(conda info --base)/envs/$ENV_NAME" --datadir="$(conda info --base)/envs/$ENV_NAME/var/lib/mysql"
mysqld_safe --datadir="$(conda info --base)/envs/$ENV_NAME/var/lib/mysql" &

# Esperar a que MariaDB se inicie
sleep 10

# Configurar MariaDB
mariadb -u root <<EOF
CREATE USER IF NOT EXISTS '$MYSQL_USER'@'localhost' IDENTIFIED BY '$MYSQL_PASSWORD';
CREATE DATABASE IF NOT EXISTS $MYSQL_DATABASE CHARACTER SET UTF8;
GRANT ALL PRIVILEGES ON $MYSQL_DATABASE.* TO '$MYSQL_USER'@'localhost';
FLUSH PRIVILEGES;
EOF

# Verificar conexión
mariadb -u $MYSQL_USER -p$MYSQL_PASSWORD -e "USE $MYSQL_DATABASE; SHOW TABLES;"

# Configuración de GLUE
export GLUE_HOME=$GLUE_HOME
export PATH=${PATH}:${GLUE_HOME}/bin

# Hacer el script ejecutable
chmod u+x ${GLUE_HOME}/bin/gluetools.sh

# Configuración del archivo XML
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

# Crear directorios temporales
mkdir -p ${GLUE_HOME}/tmp/blastfiles
mkdir -p ${GLUE_HOME}/tmp/blastdbs
mkdir -p ${GLUE_HOME}/tmp/mafftfiles
mkdir -p ${GLUE_HOME}/tmp/raxmlfiles

# Mensaje final
echo "Instalación de GLUE completada con éxito. Ejecuta '${GLUE_HOME}/bin/gluetools.sh' para iniciar GLUE."

# Ejecutar GLUE para verificar que funciona
${GLUE_HOME}/bin/gluetools.sh
