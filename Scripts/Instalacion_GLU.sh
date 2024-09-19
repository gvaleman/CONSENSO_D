#!/bin/bash

# Variables
ENV_NAME="CONSENSO_D"
GLUE_HOME="$HOME/CONSENSO_D/DENV_GLU/gluetools"
BLAST_DIR="$HOME/CONSENSO_D/blast-2.2.31-pl526he19e7b1_5"
MYSQL_USER="gluetools"
MYSQL_PASSWORD="glue12345"
MYSQL_DATABASE="GLUE_TOOLS"

# Activar entorno Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $ENV_NAME

# Instalar MySQL si no está instalado
if ! mysql --version; then
    sudo apt-get update
    sudo apt-get install -y mysql-server
    sudo sed -i 's/^#socket/socket/' /etc/mysql/mysql.conf.d/mysqld.cnf
    sudo sed -i 's/^#pid-file/pid-file/' /etc/mysql/mysql.conf.d/mysqld.cnf
    sudo sed -i 's/^#port/port/' /etc/mysql/mysql.conf.d/mysqld.cnf
    sudo service mysql restart
    sleep 5
fi

# Configurar MySQL
sudo mysql <<EOF
CREATE USER IF NOT EXISTS '$MYSQL_USER'@'localhost' IDENTIFIED BY '$MYSQL_PASSWORD';
CREATE DATABASE IF NOT EXISTS $MYSQL_DATABASE CHARACTER SET UTF8;
GRANT ALL PRIVILEGES ON $MYSQL_DATABASE.* TO '$MYSQL_USER'@'localhost';
FLUSH PRIVILEGES;
EOF

# Verificar conexión
mysql -u $MYSQL_USER -p$MYSQL_PASSWORD -e "USE $MYSQL_DATABASE; SHOW TABLES;" --socket=/var/run/mysqld/mysqld.sock

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
    </properties>
</gluetools>
EOL

# Crear directorios temporales
mkdir -p ${GLUE_HOME}/tmp/blastfiles
mkdir -p ${GLUE_HOME}/tmp/blastdbs

# Mensaje final
echo "Instalación de GLUE completada con éxito. Ejecuta '${GLUE_HOME}/bin/gluetools.sh' para iniciar GLUE."

# Ejecutar GLUE para verificar que funciona
${GLUE_HOME}/bin/gluetools.sh
