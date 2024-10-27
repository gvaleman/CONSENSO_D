#!/bin/bash

# Verificar si expect está instalado
if ! command -v expect &> /dev/null; then
    echo "Error: 'expect' no está instalado. Por favor, instálalo con:"
    echo "sudo apt-get install expect"
    exit 1
fi

# Verificar argumentos
if [ "$#" -ne 3 ]; then
    echo "Uso: $0 <archivo_fasta> <directorio_salida> <serotipo>"
    echo "Ejemplo: $0 /ruta/secuencias.fasta /ruta/resultados DENV_3"
    exit 1
fi

# Asignar argumentos a variables
FASTA_FILE="$1"
OUTPUT_DIR="$2"
SEROTYPE="$3"

# Validar serotipo
if ! [[ "$SEROTYPE" =~ ^DENV_[1-4]$ ]]; then
    echo "Error: El serotipo debe ser uno de: DENV_1, DENV_2, DENV_3, DENV_4"
    exit 1
fi

# Crear directorio de salida si no existe
mkdir -p "$OUTPUT_DIR"

# Extraer el número del serotipo para el comando
SEROTYPE_NUM=$(echo $SEROTYPE | cut -d'_' -f2)

# Crear script expect
cat << 'EOF' > /tmp/glue_commands.exp
#!/usr/bin/expect -f

# Deshabilitar timeout automático para permitir la finalización completa
set timeout -1

# Obtener argumentos
set fasta_file [lindex $argv 0]
set output_dir [lindex $argv 1]
set serotype_num [lindex $argv 2]

# Función para manejar la paginación
proc handle_pagination {} {
    expect {
        -re {Rows.*Q:quit} {
            send "Q\r"
            exp_continue
        }
        "Mode path: /project/dengue" {
            return
        }
        timeout {
            puts "Timeout durante la paginación"
            exit 1
        }
    }
}

# Iniciar gluetools
spawn gluetools.sh

# Esperar al prompt de GLUE y enviar comandos
expect {
    timeout { puts "Timeout esperando GLUE>"; exit 1 }
    "GLUE>" { send "project dengue\r" }
}

expect {
    timeout { puts "Timeout en project dengue"; exit 1 }
    "Mode path: /project/dengue" { send "console set cmd-output-file-format tab\r" }
}

expect {
    timeout { puts "Timeout en set format"; exit 1 }
    "Mode path: /project/dengue" { send "console set next-cmd-output-file $output_dir/resultados_genotipo_DENV_$serotype_num.txt\r" }
}

expect {
    timeout { puts "Timeout en set output file"; exit 1 }
    "Mode path: /project/dengue" { send "module denv${serotype_num}MaxLikelihoodGenotyper genotype file -f '$fasta_file'\r" }
}

# Manejar la paginación y esperar a que termine
handle_pagination

# Salir limpiamente
send "exit\r"
expect {
    timeout { puts "Timeout en exit"; exit 1 }
    eof { exit 0 }
}
EOF

# Dar permisos de ejecución al script expect
chmod +x /tmp/glue_commands.exp

echo "----------------------------------------------------------------------------"
echo "--   Iniciando análisis de genotipificación para DENV-$SEROTYPE_NUM...    --"
echo "---------------------------------------------------------------------------"
echo " "

# Ejecutar el script expect
/tmp/glue_commands.exp "$FASTA_FILE" "$OUTPUT_DIR" "$SEROTYPE_NUM" 2>&1 | tee /tmp/glue_output.log

# Capturar el código de salida
EXPECT_EXIT_CODE=$?

# Verificar si la ejecución fue exitosa y mostrar mensaje personalizado
if [ $EXPECT_EXIT_CODE -eq 0 ] && [ -f "${OUTPUT_DIR}/resultados_genotipo_${SEROTYPE}.txt" ]; then
    echo "✓ Genotipificación completada exitosamente"
    echo "✓ Resultados guardados en: ${OUTPUT_DIR}/resultados_genotipo_${SEROTYPE}.txt"
    
    # Mostrar un resumen de los resultados
    echo -e "\nResumen de resultados:"
    echo "------------------------"
    echo "Total de secuencias analizadas: $(grep -c ">" "${OUTPUT_DIR}/resultados_genotipo_${SEROTYPE}.txt")"
    echo -e "\nPuedes encontrar los resultados completos en:"
    echo "${OUTPUT_DIR}/resultados_genotipo_${SEROTYPE}.txt"
else
    echo "✗ Error durante la ejecución de GLUE"
    echo "Consulta el log de errores en: /tmp/glue_output.log"
    exit 1
fi

# Limpiar archivos temporales
rm /tmp/glue_commands.exp /tmp/glue_output.log
