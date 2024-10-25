#!/bin/bash

# Variables
SCRIPTS_DIR="$HOME/CONSENSO_D/Scripts"
REPO_URL="https://github.com/giffordlabcvr/Dengue-GLUE.git"
REPO_DIR="$HOME/CONSENSO_D/Dengue-GLUE"

# Clonar el repositorio si no existe
if [ ! -d "${REPO_DIR}" ]; then
    echo "Clonando el repositorio Dengue-GLUE en ${REPO_DIR}..."
    git clone $REPO_URL $REPO_DIR
else
    echo "El repositorio Dengue-GLUE ya está clonado en ${REPO_DIR}."
fi

# Hacer los scripts ejecutables
if [ -f "$SCRIPTS_DIR/Main.py" ]; then
    chmod +x $SCRIPTS_DIR/Main.py
    echo "Permisos de ejecución asignados a Main.py"
else
    echo "Main.py no encontrado en $SCRIPTS_DIR"
fi

if [ -f "$SCRIPTS_DIR/CONSENSO" ]; then
    chmod +x $SCRIPTS_DIR/CONSENSO
    echo "Permisos de ejecución asignados a CONSENSO"
else
    echo "CONSENSO no encontrado en $SCRIPTS_DIR"
fi

# Mensaje final
echo "Instalación simplificada completada con éxito."
