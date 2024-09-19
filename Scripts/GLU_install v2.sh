#!/bin/bash

# Variables de entorno
CONDA_ENV="CONSENSO_D"
GLUE_DIR="/home/ics2/CONSENSO_D/DENV_GLU"
GLUE_INSTALL_DIR="$GLUE_DIR/gluetools"
GLUE_ZIP="$GLUE_DIR/glueInstallDir-1.1.113.zip"
GLUE_JAR="$GLUE_DIR/gluetools-core-1.1.113.jar"
MYSQL_ROOT_PASSWORD="root_password"  # Cambia esto a la contraseña de root de MySQL
GLUE_DB_PASSWORD="glue12345"
MYSQL_DATA_DIR="$HOME/mysql_data"
MYSQL_LOG_DIR="$HOME/mysql_logs"
MYSQL_SOCKET_DIR="$HOME/mysql_socket"
MYSQL_CNF="$HOME/mysql_config/my.cnf"

# Crear directorios necesarios
mkdir -p "$MYSQL_DATA_DIR" "$MYSQL_LOG_DIR" "$MYSQL_SOCKET_DIR" "$(dirname "$MYSQL_CNF")"

# Asegurarse de que el script se detenga en caso de error
set -e

# Activar el entorno de Conda
echo "Activando el entorno de Conda..."
source ~/anaconda3/etc/profile.d/conda.sh || { echo "Error: No se pudo cargar Conda."; exit 1; }
conda activate "$CONDA_ENV" || { echo "Error: No se pudo activar el entorno Conda."; exit 1; }

# Paso 1: Descomprimir el archivo ZIP de GLUE
echo "Descomprimiendo gluetools..."
unzip -o "$GLUE_ZIP" -d "$GLUE_DIR"

# Paso 2: Mover el archivo JAR del motor de GLUE a la carpeta lib
echo "Moviendo el archivo JAR del motor de GLUE a la carpeta lib..."
mkdir -p "$GLUE_INSTALL_DIR/lib"
mv "$GLUE_JAR" "$GLUE_INSTALL_DIR/lib/" || { echo "Error: No se pudo mover el archivo JAR."; exit 1; }

# Paso 3: Configurar variables de entorno
echo "Configurando variables de entorno..."
echo "export GLUE_HOME=$GLUE_INSTALL_DIR" >> ~/.bashrc
echo "export PATH=\$PATH:\$GLUE_HOME/bin" >> ~/.bashrc
source ~/.bashrc || { echo "Error: No se pudo recargar las variables de entorno."; exit 1; }

# Paso 4: Asegurarse de que el script GLUE es ejecutable
echo "Haciendo el script GLUE ejecutable..."
chmod u+x "$GLUE_INSTALL_DIR/bin/gluetools.sh" || { echo "Error: No se pudo hacer ejecutable el script GLUE."; exit 1; }

# Paso 5: Instalar dependencias en el entorno de Conda
echo "Instalando dependencias en el entorno de Conda..."
conda install -y -c conda-forge openjdk=8 || { echo "Error: No se pudo instalar OpenJDK."; exit 1; }
conda install -y -c anaconda mysql || { echo "Error: No se pudo instalar MySQL."; exit 1; }

# Paso 6: Configurar MySQL
echo "Configurando MySQL..."
chmod 750 "$MYSQL_DATA_DIR" "$MYSQL_LOG_DIR" "$MYSQL_SOCKET_DIR"

# Crear archivo de configuración de MySQL
echo "[mysqld]" > "$MYSQL_CNF"
echo "datadir=$MYSQL_DATA_DIR" >> "$MYSQL_CNF"
echo "socket=$MYSQL_SOCKET_DIR/mysql.sock" >> "$MYSQL_CNF"
echo "log-error=$MYSQL_LOG_DIR/error.log" >> "$MYSQL_CNF"
echo "pid-file=$MYSQL_DATA_DIR/mysql.pid" >> "$MYSQL_CNF"
echo "lc-messages-dir=$(dirname $(which mysqld))/share" >> "$MYSQL_CNF"
echo "[client]" >> "$MYSQL_CNF"
echo "socket=$MYSQL_SOCKET_DIR/mysql.sock" >> "$MYSQL_CNF"
echo "user=root" >> "$MYSQL_CNF"

# Asegurarse de que el directorio de datos tenga los permisos correctos
sudo chown -R $(whoami):$(whoami) "$MYSQL_DATA_DIR" "$MYSQL_LOG_DIR" "$MYSQL_SOCKET_DIR"

# Inicializar MySQL con el archivo de configuración personalizado
mysqld --defaults-file="$MYSQL_CNF" --initialize-insecure || { echo "Error: No se pudo inicializar MySQL."; exit 1; }

# Iniciar MySQL utilizando el archivo de configuración personalizado
mysqld_safe --defaults-file="$MYSQL_CNF" &

# Esperar a que MySQL inicie
sleep 10

# Verificar si MySQL está corriendo
if ! mysqladmin ping --socket="$MYSQL_SOCKET_DIR/mysql.sock" --silent; then
    echo "Error: MySQL no está corriendo."
    echo "Revisar el archivo de log para más detalles: $MYSQL_LOG_DIR/error.log"
    cat "$MYSQL_LOG_DIR/error.log"
    exit 1
fi

mysql --socket="$MYSQL_SOCKET_DIR/mysql.sock" -e "
CREATE USER 'gluetools'@'localhost' IDENTIFIED BY '$GLUE_DB_PASSWORD';
CREATE DATABASE GLUE_TOOLS CHARACTER SET UTF8;
GRANT ALL PRIVILEGES ON GLUE_TOOLS.* TO 'gluetools'@'localhost';
" || { echo "Error: No se pudo configurar la base de datos MySQL."; exit 1; }

# Paso 7: Configurar el archivo XML de GLUE
echo "Configurando el archivo XML de GLUE..."
GLUE_CONFIG_FILE="$GLUE_INSTALL_DIR/conf/gluetools-config.xml"

# Editar el archivo XML para ajustar las rutas
sed -i 's|/home/fred|/home/ics2|g' "$GLUE_CONFIG_FILE"

# Asegurarse de que las rutas de los directorios temporales existen
mkdir -p /home/ics2/gluetools/tmp/blastfiles
mkdir -p /home/ics2/gluetools/tmp/blastdbs
mkdir -p /home/ics2/gluetools/tmp/raxmlfiles
mkdir -p /home/ics2/gluetools/tmp/mafftfiles
mkdir -p /home/ics2/gluetools/tmp/samfiles

# Paso 8: Probar el funcionamiento de GLUE
echo "Probando el funcionamiento de GLUE..."
$GLUE_INSTALL_DIR/bin/gluetools.sh << EOF
quit
EOF

echo "Instalación de GLUE finalizada."

