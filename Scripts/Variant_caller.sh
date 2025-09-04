#!/bin/bash

# Configuración
serotype=$1
input_fasta=$2
max_parallel=19  # Ajustar según tu hardware

# Validaciones iniciales
if [ $# -ne 2 ]; then
    echo "CONSENSO_D - Variant Calling Pipeline"
    echo "Uso: $0 <VIRUS/SEROTIPO> <ARCHIVO_FASTA>"
    echo ""
    echo "Virus/Serotipos soportados:"
    echo "  DENV_1, DENV_2, DENV_3, DENV_4 - Virus del Dengue"
    echo "  RABV                           - Virus de la Rabia"
    echo ""
    echo "Funcionalidad:"
    echo "  - Procesa archivos multiFASTA con genomas consenso"
    echo "  - Genera archivos VCF individuales por secuencia"
    echo "  - Paralelización automática para archivos grandes"
    exit 1
fi

# Obtener rutas
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
consenso_dir="$(dirname "$script_dir")"  # Directorio padre de Scripts
input_dir=$(dirname "$input_fasta")

# Definir referencia (buscar en el directorio CONSENSO_D)
case "$serotype" in
    "DENV_1") reference_fasta="$consenso_dir/Ref_DENV/Reference_DV_1.fasta" ;;
    "DENV_2") reference_fasta="$consenso_dir/Ref_DENV/Reference_DV_2.fasta" ;;
    "DENV_3") reference_fasta="$consenso_dir/Ref_DENV/Reference_DV_3.fasta" ;;
    "DENV_4") reference_fasta="$consenso_dir/Ref_DENV/Reference_DV_4.fasta" ;;
    "RABV") reference_fasta="$consenso_dir/Ref_DENV/RABV_Reference.fasta" ;;
    *) echo "Error: Serotipo no válido: $serotype"; exit 1 ;;
esac

# Validar archivos
for file in "$input_fasta" "$reference_fasta"; do
    if [ ! -f "$file" ]; then
        echo "Error: Archivo no encontrado: $file"
        exit 1
    fi
done

echo "================================================="
echo "       CONSENSO_D - VARIANT CALLER"
echo "           Virus/Serotipo: $serotype"
echo "           Archivo: $(basename $input_fasta)"
echo "           Procesamiento paralelo: $max_parallel"
echo "================================================="

# Crear directorios
result_dir="$input_dir/variant_calling_analysis"
mkdir -p "$result_dir"/{vcfs,temp,logs}

# Indexar referencia si no existe
if [ ! -f "$reference_fasta.fai" ]; then
    echo "Indexando referencia..."
    samtools faidx "$reference_fasta"
fi

# Función de procesamiento
process_sequence() {
    local seq_id=$1
    local seq_file=$2
    local ref=$3
    local out_dir=$4
    
    echo "Procesando: $seq_id" >> "$out_dir/logs/progress.log"
    
    # Alineamiento
    minimap2 -a "$ref" "$seq_file" > "$out_dir/temp/aln_$seq_id.sam" 2>/dev/null
    
    # Conversión y ordenamiento
    samtools view -S -b "$out_dir/temp/aln_$seq_id.sam" | \
    samtools sort -o "$out_dir/temp/sorted_$seq_id.bam" 2>/dev/null
    
    # Indexar BAM
    samtools index "$out_dir/temp/sorted_$seq_id.bam" 2>/dev/null
    
    # Variant calling
    bcftools mpileup -Ou -f "$ref" "$out_dir/temp/sorted_$seq_id.bam" 2>/dev/null | \
    bcftools call -mv -Oz -o "$out_dir/vcfs/${seq_id}.vcf.gz" 2>/dev/null
    
    # Limpiar archivos temporales
    rm -f "$out_dir/temp/aln_$seq_id.sam" "$out_dir/temp/sorted_$seq_id.bam"*
    
    echo "Completado: $seq_id" >> "$out_dir/logs/progress.log"
}

# Procesar multiFASTA
seq_count=0
current_seq=""
seq_id=""

while IFS= read -r line; do
    if [[ $line == ">"* ]]; then
        # Procesar secuencia anterior si existe
        if [[ -n "$current_seq" && -n "$seq_id" ]]; then
            echo "$current_seq" > "$result_dir/temp/sequence_$seq_id.fasta"
            
            # Ejecutar en paralelo
            process_sequence "$seq_id" "$result_dir/temp/sequence_$seq_id.fasta" \
                           "$reference_fasta" "$result_dir" &
            
            # Control de paralelización
            if (( seq_count % max_parallel == 0 )); then
                wait
            fi
        fi
        
        # Preparar nueva secuencia
        seq_count=$((seq_count + 1))
        seq_id=$(echo "$line" | sed 's/>//' | cut -d' ' -f1 | tr '/' '_')
        current_seq="$line"
        
        echo "Preparando secuencia $seq_count: $seq_id"
        
    else
        current_seq="$current_seq"$'\n'"$line"
    fi
done < "$input_fasta"

# Procesar última secuencia
if [[ -n "$current_seq" && -n "$seq_id" ]]; then
    echo "$current_seq" > "$result_dir/temp/sequence_$seq_id.fasta"
    process_sequence "$seq_id" "$result_dir/temp/sequence_$seq_id.fasta" \
                   "$reference_fasta" "$result_dir" &
fi

# Esperar que terminen todos los procesos
wait

# Limpiar archivos temporales de secuencias
rm -f "$result_dir/temp/sequence_"*.fasta

echo "================================================="
echo "CONSENSO_D - PROCESO COMPLETADO"
echo "Virus/Serotipo: $serotype"
echo "Total de secuencias procesadas: $seq_count"
echo "Archivos VCF generados en: $result_dir/vcfs/"
echo "Log de progreso en: $result_dir/logs/progress.log"
echo "================================================="
