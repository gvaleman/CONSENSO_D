#!/bin/bash

# Variables
JAVA_VERSION="1.8.0"
MYSQL_VERSION="5.6"
GLUE_USER="gluetools"
GLUE_PASSWORD="glue12345"
GLUE_DB="GLUE_TOOLS"
GLUE_INSTALL_DIR="/home/ics2/gluetools"
LOCAL_GLUE_INSTALL_ZIP="/home/ics2/CONSENSO_D/DENV_GLU/glueInstallDir-1.1.113.zip"
LOCAL_GLUE_ENGINE_JAR="/home/ics2/CONSENSO_D/DENV_GLU/gluetools-core-1.1.113.jar"

# Function to check command existence
command_exists() {
    command -v "$1" &>/dev/null
}

# Install Java
if ! command_exists java || [[ "$(java -version 2>&1)" != *"$JAVA_VERSION"* ]]; then
    echo "Installing Java $JAVA_VERSION..."
    sudo apt update
    sudo apt install -y openjdk-8-jdk
    export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64
    export PATH=$JAVA_HOME/bin:$PATH
fi

# Install MySQL
if ! command_exists mysql || [[ "$(mysql --version)" != *"Distrib $MYSQL_VERSION"* ]]; then
    echo "Installing MySQL $MYSQL_VERSION..."
    sudo apt update
    sudo apt install -y mysql-server
    sudo service mysql start
fi

# Ensure MySQL is running
echo "Checking MySQL status..."
sudo service mysql status
if [[ $? -ne 0 ]]; then
    echo "Starting MySQL service..."
    sudo service mysql start
    if [[ $? -ne 0 ]]; then
        echo "Failed to start MySQL service. Please check MySQL logs."
        exit 1
    fi
fi

# Create MySQL user and database if they don't exist
mysql -u root -e "CREATE USER IF NOT EXISTS '$GLUE_USER'@'localhost' IDENTIFIED BY '$GLUE_PASSWORD';"
mysql -u root -e "CREATE DATABASE IF NOT EXISTS $GLUE_DB CHARACTER SET UTF8;"
mysql -u root -e "GRANT ALL PRIVILEGES ON $GLUE_DB.* TO '$GLUE_USER'@'localhost';"

# Install GLUE
if [ -f "$LOCAL_GLUE_INSTALL_ZIP" ]; then
    echo "Installing GLUE..."
    unzip -o "$LOCAL_GLUE_INSTALL_ZIP" -d "$(dirname "$GLUE_INSTALL_DIR")"
    export GLUE_HOME="$GLUE_INSTALL_DIR"
    export PATH="$PATH:$GLUE_HOME/bin"
    echo "export GLUE_HOME=$GLUE_HOME" >> ~/.bashrc
    echo "export PATH=\$PATH:\$GLUE_HOME/bin" >> ~/.bashrc
else
    echo "GLUE installation zip not found at $LOCAL_GLUE_INSTALL_ZIP."
    exit 1
fi

# Place GLUE engine jar
if [ -f "$LOCAL_GLUE_ENGINE_JAR" ]; then
    mkdir -p "$GLUE_HOME/lib"
    cp "$LOCAL_GLUE_ENGINE_JAR" "$GLUE_HOME/lib/gluetools-core-1.1.113.jar"
else
    echo "GLUE engine jar not found at $LOCAL_GLUE_ENGINE_JAR."
    exit 1
fi

# Configure GLUE
GLUE_CONF="$GLUE_HOME/conf/gluetools-config.xml"
if [ -f "$GLUE_CONF" ]; then
    echo "Configuring GLUE..."
    sed -i "s/<username>.*<\/username>/<username>$GLUE_USER<\/username>/" "$GLUE_CONF"
    sed -i "s/<password>.*<\/password>/<password>$GLUE_PASSWORD<\/password>/" "$GLUE_CONF"
    sed -i "s|<jdbcUrl>.*<\/jdbcUrl>|<jdbcUrl>jdbc:mysql://localhost:3306/$GLUE_DB?useSSL=false&amp;characterEncoding=UTF-8<\/jdbcUrl>|" "$GLUE_CONF"
    sed -i "/<vendor>.*<\/vendor>/d" "$GLUE_CONF"  # Remove deprecated vendor element
else
    echo "Failed to find GLUE configuration file at $GLUE_CONF."
    exit 1
fi

# Make sure the GLUE bash script is executable
if [ -f "$GLUE_HOME/bin/gluetools.sh" ]; then
    chmod u+x "$GLUE_HOME/bin/gluetools.sh"
else
    echo "Failed to find GLUE tools script at $GLUE_HOME/bin/gluetools.sh."
    exit 1
fi

# Testing GLUE command line
echo "Testing GLUE command line..."
if command -v gluetools.sh &>/dev/null; then
    gluetools.sh <<EOF
quit
EOF
else
    echo "Failed to run GLUE command line. Please check the installation."
    exit 1
fi

echo "GLUE installation and configuration complete."
