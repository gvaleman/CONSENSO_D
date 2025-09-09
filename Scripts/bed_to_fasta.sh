#!/bin/bash

# Script para convertir archivos BED a FASTA de primers
# Uso: ./bed_to_fasta.sh <reference_fasta> <bed_file> <output_fasta>

if [ $# -ne 3 ]; then
    echo "Uso: $0 <reference_fasta> <bed_file> <output_fasta>"
    echo ""
    echo "Ejemplos:"
    echo "  $0 reference.fasta primers.bed primers.fasta"
    echo "  $0 /path/to/ref.fasta /path/to/primers.bed output_primers.fasta"
    exit 1
fi

REFERENCE_FASTA="$1"
BED_FILE="$2"
OUTPUT_FASTA="$3"

# Verificar que los archivos de entrada existan
if [ ! -f "$REFERENCE_FASTA" ]; then
    echo "Error: Archivo de referencia no encontrado: $REFERENCE_FASTA"
    exit 1
fi

if [ ! -f "$BED_FILE" ]; then
    echo "Error: Archivo BED no encontrado: $BED_FILE"
    exit 1
fi

echo "Convirtiendo BED a FASTA..."
echo "Referencia: $REFERENCE_FASTA"
echo "BED file: $BED_FILE"
echo "Output: $OUTPUT_FASTA"

# Usar bedtools para extraer las secuencias
if command -v bedtools &> /dev/null; then
    echo "Usando bedtools para extraer secuencias..."
    bedtools getfasta -fi "$REFERENCE_FASTA" -bed "$BED_FILE" -fo "$OUTPUT_FASTA" -name
    
    if [ $? -eq 0 ]; then
        echo "✅ Conversión completada exitosamente"
        echo "📊 Primers extraídos: $(grep -c "^>" "$OUTPUT_FASTA")"
        echo "📁 Archivo guardado en: $OUTPUT_FASTA"
    else
        echo "❌ Error durante la extracción con bedtools"
        exit 1
    fi
else
    echo "⚠️  bedtools no está instalado. Usando método alternativo con Python..."
    
    # Método alternativo usando Python
    python3 << EOF
import sys

def read_fasta(filename):
    """Lee un archivo FASTA y retorna un diccionario con las secuencias"""
    sequences = {}
    current_seq = ""
    current_name = ""
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = current_seq
                current_name = line[1:].split()[0]  # Solo el primer campo después de >
                current_seq = ""
            else:
                current_seq += line
    
    if current_name:
        sequences[current_name] = current_seq
    
    return sequences

def process_bed_to_fasta(reference_fasta, bed_file, output_fasta):
    """Convierte un archivo BED a FASTA usando una referencia"""
    
    # Leer archivo de referencia
    print("Leyendo archivo de referencia...")
    ref_sequences = read_fasta(reference_fasta)
    
    # Procesar archivo BED
    primer_count = 0
    with open(bed_file, 'r') as bed, open(output_fasta, 'w') as fasta:
        for line in bed:
            line = line.strip()
            if line and not line.startswith('#'):
                fields = line.split('\t')
                if len(fields) >= 3:
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    
                    # Nombre del primer (usar campo 4 si existe, sino crear uno)
                    if len(fields) > 3 and fields[3]:
                        primer_name = fields[3]
                    else:
                        primer_name = f"primer_{primer_count + 1}_{chrom}_{start}_{end}"
                    
                    # Extraer secuencia
                    if chrom in ref_sequences:
                        ref_seq = ref_sequences[chrom]
                        if start < len(ref_seq) and end <= len(ref_seq):
                            primer_seq = ref_seq[start:end]
                            fasta.write(f">{primer_name}\n{primer_seq}\n")
                            primer_count += 1
                        else:
                            print(f"⚠️  Coordenadas fuera de rango para {chrom}: {start}-{end}")
                    else:
                        print(f"⚠️  Cromosoma no encontrado en referencia: {chrom}")
    
    return primer_count

# Ejecutar la conversión
try:
    count = process_bed_to_fasta("$REFERENCE_FASTA", "$BED_FILE", "$OUTPUT_FASTA")
    print(f"✅ Conversión completada exitosamente")
    print(f"📊 Primers extraídos: {count}")
    print(f"📁 Archivo guardado en: $OUTPUT_FASTA")
except Exception as e:
    print(f"❌ Error durante la conversión: {e}")
    sys.exit(1)
EOF

fi

echo ""
echo "Verificando el archivo de salida..."
if [ -f "$OUTPUT_FASTA" ]; then
    echo "📋 Resumen del archivo generado:"
    echo "   - Total de secuencias: $(grep -c "^>" "$OUTPUT_FASTA")"
    echo "   - Tamaño del archivo: $(du -h "$OUTPUT_FASTA" | cut -f1)"
    echo "   - Primeras 3 secuencias:"
    head -6 "$OUTPUT_FASTA" | grep -E "^>|^[ATCGN]" | head -6
else
    echo "❌ Error: No se pudo generar el archivo de salida"
    exit 1
fi
