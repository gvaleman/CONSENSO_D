#!/bin/bash

#Este Script genera secuencias de Referencias del dengue con filtro de calidad 50x y ensamblaje multi-serotipo

# Salir a la carpeta home para asegurar que las rutas relativas sean correctas
cd

# Obtener los argumentos
sequence_type=$1
serotype=$2
dir=$3
primer_bed_file=$4

# Obtener la ruta del directorio actual (la carpeta home del usuario)
current_dir=$PWD

# Redirigir la salida a un archivo temporal en la carpeta home del usuario
LOG_FILE_TEMP="$HOME/proceso_log_temp.txt"
exec > >(tee -a "$LOG_FILE_TEMP") 2>&1

# Calcular número óptimo de hilos (80% de los disponibles, mínimo 1)
available_threads=$(nproc)
optimal_threads=$(echo "scale=0; ($available_threads * 0.8)/1" | bc)
threads=${optimal_threads:-1}
threads=${threads%%.*}

echo "Usando $threads hilos para el procesamiento"

# Detectar la ubicación de los scripts
if [ -f "$current_dir/CONSENSO_D/Scripts/DENV_Universal_Classifier.py" ]; then
    classifier_script="$current_dir/CONSENSO_D/Scripts/DENV_Universal_Classifier.py"
    primer_extractor_script="$current_dir/CONSENSO_D/Scripts/extract_primers.py"
elif [ -f "$current_dir/Scripts/DENV_Universal_Classifier.py" ]; then
    classifier_script="$current_dir/Scripts/DENV_Universal_Classifier.py"
    primer_extractor_script="$current_dir/Scripts/extract_primers.py"
elif [ -f "./DENV_Universal_Classifier.py" ]; then
    classifier_script="./DENV_Universal_Classifier.py"
    primer_extractor_script="./extract_primers.py"
else
    classifier_script=""
    primer_extractor_script=""
fi

# Definir todas las referencias y modelos de DENV
declare -A DENV_REFS
declare -A DENV_MODELS
declare -A DENV_ENCODERS

DENV_REFS["DENV_1"]="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_1.fasta"
DENV_REFS["DENV_2"]="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_2.fasta"
DENV_REFS["DENV_3"]="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_3.fasta"
DENV_REFS["DENV_4"]="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_4.fasta"

DENV_MODELS["DENV_1"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV1/dengue_lineage_classifier.joblib.gz"
DENV_MODELS["DENV_2"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV2/dengue_lineage_classifier.joblib.gz"
DENV_MODELS["DENV_3"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV3/dengue_lineage_classifier.joblib.gz"
DENV_MODELS["DENV_4"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV4/dengue_lineage_classifier.joblib.gz"

DENV_ENCODERS["DENV_1"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV1/label_encoder.joblib"
DENV_ENCODERS["DENV_2"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV2/label_encoder.joblib"
DENV_ENCODERS["DENV_3"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV3/label_encoder.joblib"
DENV_ENCODERS["DENV_4"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV4/label_encoder.joblib"

# Variables adicionales para otros virus
case "$serotype" in
    "SARS_COV_2")
        fasta_file="$current_dir/CONSENSO_D/Ref_DENV/REF_NC_045512_SARS_COV_2.fasta"
        ;;
    "RABV")
        fasta_file="$current_dir/CONSENSO_D/Ref_DENV/RABV_Reference.fasta"
        ;;
    "N_RABV")
        fasta_file="$current_dir/CONSENSO_D/Ref_DENV/Lisa_RABV_N.fasta"
        ;;
    "ALL_DENV")
        echo "Modo multi-serotipo activado: se ensamblarán todos los serotipos de DENV"
        ;;
    "DENV_1"|"DENV_2"|"DENV_3"|"DENV_4")
        fasta_file="${DENV_REFS[$serotype]}"
        ;;
    *)
        echo "Error: serotipo / Virus no válido"
        echo "
Intentar con las opciones:
DENV_1, DENV_2, DENV_3, DENV_4
ALL_DENV (para ensamblar todos los serotipos)
SARS_COV_2
RABV
N_RABV"
        exit 1
        ;;
esac

# Colores
CYAN="\e[36m"
YELLOW="\e[33m"
GREEN="\e[32m"
RED="\e[31m"
WHITE="\e[97m"
GRAY="\e[90m"
BOLD="\e[1m"
RESET="\e[0m"

echo -e "${CYAN}===========================================${RESET}"
echo -e "${BOLD}${WHITE}  CONSENSO ${GRAY}v0.06-beta${RESET}"
echo -e "${CYAN}===========================================${RESET}"
echo ""
echo -e "  ${BOLD}Ensamblando${RESET}"
echo -e "  Serotipo: ${YELLOW}${BOLD}${serotype}${RESET}"
echo -e "  ${GREEN}Filtro de calidad: ≥50x cobertura${RESET}"
echo ""
echo -e "${CYAN}-------------------------------------------${RESET}"
echo -e "  ${GRAY}Instituto De Ciencias Sostenibles${RESET}"
echo -e "${CYAN}===========================================${RESET}"

if [ ! -d "$dir" ]; then
    echo "No se puede acceder al directorio: $dir"
    exit 1
fi

cd "$dir"

# Crear directorios para resultados
mkdir -p QC_Reports
mkdir -p Lineage_Classification
mkdir -p Multi_Serotype_Results

# Preparar archivos de adaptadores
if [ -f "$current_dir/CONSENSO_D/Adapters/adapters.csv" ]; then
    awk -F',' 'NR>1 && $1!="" && $2!="" {print ">"$1"\n"$2}' "$current_dir/CONSENSO_D/Adapters/adapters.csv" > adapters.fasta
    ADAPTERS_FASTA="adapters.fasta"
else
    ADAPTERS_FASTA=""
fi

# Preparar archivo de primers
PRIMERS_FASTA=""
primer_file_path=$4

if [ -n "$primer_file_path" ] && [ "$primer_file_path" != "none" ] && [ -f "$primer_file_path" ]; then
    echo "Usando archivo de primers proporcionado: $primer_file_path"
    PRIMERS_FASTA="$primer_file_path"
fi

# Función para calcular estadísticas de cobertura y calidad
calculate_coverage_stats() {
    local coverage_file=$1
    local output_prefix=$2
    
    awk '{
        total_pos++;
        sum_cov += $3;
        if ($3 >= 50) cov_50x++;
        if ($3 >= 100) cov_100x++;
        if ($3 > max_cov) max_cov = $3;
    } END {
        avg_cov = (total_pos > 0) ? sum_cov/total_pos : 0;
        pct_50x = (total_pos > 0) ? (cov_50x/total_pos)*100 : 0;
        pct_100x = (total_pos > 0) ? (cov_100x/total_pos)*100 : 0;
        print avg_cov","pct_50x","pct_100x","max_cov","total_pos
    }' "$coverage_file" > "${output_prefix}_stats.txt"
}

# Función para contar bases no ambiguas
count_non_ambiguous_bases() {
    local fasta_file=$1
    
    # Contar bases no ambiguas (A, T, G, C)
    grep -v "^>" "$fasta_file" | tr -d '\n' | grep -o "[ATGC]" | wc -l
}

# Función para ensamblar contra una referencia específica
assemble_against_reference() {
    local ref_name=$1
    local ref_file=$2
    local input_file=$3
    local sample_name=$4
    local is_paired=$5
    local R2_file=$6
    
    echo "  → Ensamblando contra $ref_name..."
    
    # Crear subdirectorio para este serotipo
    local output_dir="${sample_name}_${ref_name}"
    mkdir -p "$output_dir"
    cd "$output_dir"
    
    if [ "$is_paired" == "true" ]; then
        # Mapeo paired-end
        bwa mem -t $threads "$ref_file" "$input_file" "$R2_file" > "Output_${ref_name}.aln.sam"
    else
        # Mapeo single-end o Nanopore
        minimap2 -t $threads -a "$ref_file" "$input_file" > "Output_${ref_name}.aln.sam"
    fi
    
    # Proceso común
    samtools view -@ $threads -S -b "Output_${ref_name}.aln.sam" > "Output_${ref_name}.aln.bam"
    samtools sort -@ $threads "Output_${ref_name}.aln.bam" -o "Output_${ref_name}.sorted.bam"
    samtools index "Output_${ref_name}.sorted.bam"
    
    # Verificar calidad
    samtools flagstat "Output_${ref_name}.sorted.bam" > "${sample_name}_${ref_name}.flagstat"
    samtools stats "Output_${ref_name}.sorted.bam" > "${sample_name}_${ref_name}.stats"
    
    # Cobertura
    samtools depth "Output_${ref_name}.sorted.bam" > "${sample_name}_${ref_name}.coverage"
    
    # Calcular estadísticas
    calculate_coverage_stats "${sample_name}_${ref_name}.coverage" "${sample_name}_${ref_name}"
    
    # Leer estadísticas
    read avg_cov pct_50x pct_100x max_cov total_pos < <(cat "${sample_name}_${ref_name}_stats.txt" | tr ',' ' ')
    
    echo "     Cobertura promedio: ${avg_cov}x"
    echo "     Posiciones ≥50x: ${pct_50x}%"
    
    # Generar consenso solo si cumple el filtro de 50x
    local pass_filter="NO"
    if (( $(echo "$pct_50x >= 80" | bc -l) )); then
        pass_filter="YES"
        echo "     ✓ PASA filtro de calidad (≥80% posiciones con 50x)"
        
        # Generar VCF y consenso
        bcftools mpileup -Ou -f "$ref_file" "Output_${ref_name}.sorted.bam" | \
            bcftools call -c -Oz -o "${sample_name}_${ref_name}_calls.vcf.gz"
        
        bcftools norm -f "$ref_file" "${sample_name}_${ref_name}_calls.vcf.gz" \
            -Oz -o "${sample_name}_${ref_name}_normalized.vcf.gz"
        
        bcftools view -i 'QUAL>30' "${sample_name}_${ref_name}_normalized.vcf.gz" | \
            vcfutils.pl vcf2fq | \
            seqtk seq -aQ64 -q30 -n N > "SAMPLE_${ref_name}_cns.fasta"
        
        echo ">${sample_name}_${ref_name}" > "${sample_name}_${ref_name}.fasta"
        tail -n +2 "SAMPLE_${ref_name}_cns.fasta" >> "${sample_name}_${ref_name}.fasta"
        
        # Contar bases no ambiguas
        non_ambiguous=$(count_non_ambiguous_bases "${sample_name}_${ref_name}.fasta")
    else
        echo "     ✗ NO PASA filtro de calidad (<80% posiciones con 50x)"
        non_ambiguous=0
    fi
    
    # Guardar resultado en archivo temporal
    echo "${sample_name},${ref_name},${avg_cov},${pct_50x},${pct_100x},${max_cov},${pass_filter},${non_ambiguous}" >> "../temp_results.csv"
    
    # Limpiar archivos temporales
    rm "Output_${ref_name}.aln.sam" "Output_${ref_name}.aln.bam" "SAMPLE_${ref_name}_cns.fasta" 2>/dev/null || true
    
    cd ..
}

# Función principal de procesamiento
process_sample() {
    local sample_dir=$1
    local sample_name=$2
    local sequence_type=$3
    local input_file=$4
    local R2_file=$5
    
    echo ""
    echo "=========================================="
    echo "  PROCESANDO: $sample_name"
    echo "=========================================="
    
    # Crear archivo temporal para resultados
    > temp_results.csv
    
    if [ "$serotype" == "ALL_DENV" ]; then
        # Ensamblar contra todos los serotipos
        for sero in "DENV_1" "DENV_2" "DENV_3" "DENV_4"; do
            if [ "$sequence_type" == "ILLUMINA" ]; then
                assemble_against_reference "$sero" "${DENV_REFS[$sero]}" "$input_file" "$sample_name" "true" "$R2_file"
            else
                assemble_against_reference "$sero" "${DENV_REFS[$sero]}" "$input_file" "$sample_name" "false" ""
            fi
        done
        
        # Analizar resultados y determinar mejor serotipo
        echo ""
        echo "  Analizando resultados..."
        
        best_serotype=""
        best_non_ambiguous=0
        best_coverage=0
        
        while IFS=',' read -r samp sero avg_cov pct_50x pct_100x max_cov pass non_amb; do
            if [ "$pass" == "YES" ]; then
                if [ "$non_amb" -gt "$best_non_ambiguous" ]; then
                    best_non_ambiguous=$non_amb
                    best_serotype=$sero
                    best_coverage=$avg_cov
                fi
            fi
        done < temp_results.csv
        
        if [ -n "$best_serotype" ]; then
            echo "  ✓ MEJOR SEROTIPO: $best_serotype"
            echo "     Bases no ambiguas: $best_non_ambiguous"
            echo "     Cobertura promedio: ${best_coverage}x"
            
            # Guardar en tabla resumen
            echo "${sample_name},${best_serotype},${best_coverage},${best_non_ambiguous}" >> "../Multi_Serotype_Results/serotype_summary.csv"
            
            # Copiar el mejor consenso al directorio principal
            cp "${sample_name}_${best_serotype}/${sample_name}_${best_serotype}.fasta" "../${sample_name}_${best_serotype}.fasta"
        else
            echo "  ✗ Ningún serotipo pasó el filtro de calidad"
            echo "${sample_name},NONE,0,0" >> "../Multi_Serotype_Results/serotype_summary.csv"
        fi
        
    else
        # Ensamblar contra un solo serotipo
        if [ "$sequence_type" == "ILLUMINA" ]; then
            assemble_against_reference "$serotype" "$fasta_file" "$input_file" "$sample_name" "true" "$R2_file"
        else
            assemble_against_reference "$serotype" "$fasta_file" "$input_file" "$sample_name" "false" ""
        fi
    fi
    
    rm temp_results.csv
}

# Inicializar tabla resumen
if [ "$serotype" == "ALL_DENV" ]; then
    echo "Sample,Best_Serotype,Avg_Coverage,Non_Ambiguous_Bases" > "Multi_Serotype_Results/serotype_summary.csv"
fi

# PROCESAMIENTO PRINCIPAL
if [ "$sequence_type" == "NANO" ]; then
    # Procesamiento Nanopore
    for file in $(ls -d */); do
        cd "$file" || continue
        sample_name=${PWD##*/}
        
        if ls *.fastq.gz 1> /dev/null 2>&1; then
            cat $(ls *.fastq.gz) > "${sample_name}.fastq.gz"
            gzip -df "${sample_name}.fastq.gz"
            muestra=$(ls *.fastq)
        elif ls *.fastq 1> /dev/null 2>&1; then
            muestra=$(ls *.fastq)
        else
            echo "No se encontraron archivos fastq en $file"
            cd ..
            continue
        fi
        
        # Remoción de primers si disponibles
        if [ -n "$PRIMERS_FASTA" ] && [ -f "$muestra" ]; then
            echo "Recortando primers..."
            mv "$muestra" "${muestra}.original"
            cutadapt -g file:"$PRIMERS_FASTA" -a file:"$PRIMERS_FASTA" \
                -o "$muestra" "${muestra}.original" --minimum-length 50 -j $threads
            rm "${muestra}.original"
        fi
        
        process_sample "." "$sample_name" "NANO" "$muestra" ""
        
        cd ..
    done
    
elif [ "$sequence_type" == "ILLUMINA" ]; then
    # Indexar referencias
    if [ "$serotype" == "ALL_DENV" ]; then
        for sero in "DENV_1" "DENV_2" "DENV_3" "DENV_4"; do
            bwa index "${DENV_REFS[$sero]}"
        done
    else
        bwa index "$fasta_file"
    fi
    
    # Procesamiento Illumina
    while IFS= read -r -d '' sample_dir; do
        cd "$sample_dir" || continue
        R1_file=$(find . -maxdepth 1 -name '*_R1_001.fastq.gz' -print -quit)
        R2_file=$(find . -maxdepth 1 -name '*_R2_001.fastq.gz' -print -quit)
        
        if [ -n "$R1_file" ] && [ -n "$R2_file" ]; then
            sample_name=$(basename "$R1_file" "_R1_001.fastq.gz")
            
            # FastQC inicial
            mkdir -p fastqc_results
            fastqc -t $threads -o fastqc_results "$R1_file" "$R2_file"
            
            # Trimming
            R1_trimmed="${sample_name}_R1_paired_trimmed.fastq.gz"
            R2_trimmed="${sample_name}_R2_paired_trimmed.fastq.gz"
            R1_unpaired="${sample_name}_R1_unpaired_trimmed.fastq.gz"
            R2_unpaired="${sample_name}_R2_unpaired_trimmed.fastq.gz"
            
            if [ -n "$ADAPTERS_FASTA" ]; then
                trimmomatic PE -threads $threads \
                    "$R1_file" "$R2_file" \
                    "$R1_trimmed" "$R1_unpaired" \
                    "$R2_trimmed" "$R2_unpaired" \
                    ILLUMINACLIP:"../$ADAPTERS_FASTA":2:30:10 \
                    SLIDINGWINDOW:4:20 MINLEN:50
            else
                trimmomatic PE -threads $threads \
                    "$R1_file" "$R2_file" \
                    "$R1_trimmed" "$R1_unpaired" \
                    "$R2_trimmed" "$R2_unpaired" \
                    SLIDINGWINDOW:4:20 MINLEN:50
            fi
            
            # Remoción de primers
            if [ -n "$PRIMERS_FASTA" ]; then
                R1_final="${sample_name}_R1_final_trimmed.fastq.gz"
                R2_final="${sample_name}_R2_final_trimmed.fastq.gz"
                
                cutadapt -g file:"$PRIMERS_FASTA" -G file:"$PRIMERS_FASTA" \
                    -o "$R1_final" -p "$R2_final" \
                    "$R1_trimmed" "$R2_trimmed" \
                    --minimum-length 50 -j $threads
                
                R1_for_mapping="$R1_final"
                R2_for_mapping="$R2_final"
            else
                R1_for_mapping="$R1_trimmed"
                R2_for_mapping="$R2_trimmed"
            fi
            
            process_sample "." "$sample_name" "ILLUMINA" "$R1_for_mapping" "$R2_for_mapping"
            
            # Limpiar archivos temporales
            rm "$R1_unpaired" "$R2_unpaired" 2>/dev/null || true
        fi
        
        cd ..
    done < <(find . -type d -mindepth 1 -maxdepth 1 -print0)
    
else
    echo "Tipo de secuenciación no válido. Use 'NANO' o 'ILLUMINA'."
    exit 1
fi

# GENERAR REPORTES FINALES
echo ""
echo "=========================================="
echo "  GENERANDO REPORTES FINALES"
echo "=========================================="

# Generar multifasta
if [ "$serotype" == "ALL_DENV" ]; then
    cat $(find . -name "*_DENV_*.fasta" -maxdepth 1) > "all_consensus_MULTI_DENV.fasta" 2>/dev/null || true
else
    cat $(find . -name "*.fasta" -maxdepth 1) > "all_consensus_$serotype.fasta" 2>/dev/null || true
fi

# Mostrar tabla resumen de serotipos
if [ "$serotype" == "ALL_DENV" ] && [ -f "Multi_Serotype_Results/serotype_summary.csv" ]; then
    echo ""
    echo "╔═══════════════════════════════════════════════════════════════╗"
    echo "║         TABLA RESUMEN DE SEROTIPOS POR MUESTRA               ║"
    echo "╚═══════════════════════════════════════════════════════════════╝"
    echo ""
    column -t -s',' "Multi_Serotype_Results/serotype_summary.csv" | head -20
    echo ""
    total_samples=$(tail -n +2 "Multi_Serotype_Results/serotype_summary.csv" | wc -l)
    passed_samples=$(tail -n +2 "Multi_Serotype_Results/serotype_summary.csv" | grep -v ",NONE," | wc -l)
    echo "Total de muestras: $total_samples"
    echo "Muestras que pasaron filtro: $passed_samples"
fi

# MultiQC final
multiqc -f -o "QC_Reports/Final_QC_Report" QC_Reports/ 2>/dev/null || true

# Limpiar
rm -f adapters.fasta primers.fasta 2>/dev/null

# Generar informe HTML
if [ -f "$current_dir/CONSENSO_D/Scripts/quality_report.py" ]; then
    python3 "$current_dir/CONSENSO_D/Scripts/quality_report.py" "$dir"
fi

# Mensaje final
GREEN='\033[0;32m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
MAGENTA='\033[0;35m'
BOLD='\033[1m'
NC='\033[0m'

cat << EOF

    ${GREEN}╔═══════════════════════════════════════════════════════════╗${NC}
    ${GREEN}║                                                           ║${NC}
    ${GREEN}║${NC}              ${BOLD}${CYAN}✓  P R O C E S O   C O M P L E T O${NC}          ${GREEN}║${NC}
    ${GREEN}║                                                           ║${NC}
    ${GREEN}╚═══════════════════════════════════════════════════════════╝${NC}

EOF

echo -e "    ${CYAN}📊 RESULTADOS${NC}"
echo -e "    ${BLUE}─────────────────────────────────────────────────────────${NC}"

if [ "$serotype" == "ALL_DENV" ]; then
    echo -e "       ${YELLOW}Consensos multi-serotipo:${NC} all_consensus_MULTI_DENV.fasta"
    echo -e "       ${YELLOW}Tabla de serotipos:${NC} Multi_Serotype_Results/serotype_summary.csv"
else
    echo -e "       ${YELLOW}Consenso multi-fasta:${NC} all_consensus_${serotype}.fasta"
fi

echo -e "       ${YELLOW}Ubicación:${NC} $dir"
echo -e "       ${GREEN}Filtro aplicado:${NC} ≥50x cobertura en ≥80% posiciones"
echo ""
echo -e "    ${CYAN}📁 CALIDAD${NC}"
echo -e "    ${BLUE}─────────────────────────────────────────────────────────${NC}"
echo -e "       ${YELLOW}Reportes:${NC} QC_Reports/"
echo -e "       ${YELLOW}Dashboard:${NC} reporte_calidad.html"
echo ""
echo -e "    ${GREEN}═══════════════════════════════════════════════════════════${NC}"
echo -e "     ${BOLD}Finalizado:${NC} ${MAGENTA}$(date '+%Y-%m-%d %H:%M:%S')${NC}"
echo -e "    ${GREEN}═══════════════════════════════════════════════════════════${NC}"
echo ""

# Mover log final
FINAL_LOG_NAME="CONSENSO_log.txt"
if [ -f "$LOG_FILE_TEMP" ]; then
    mv "$LOG_FILE_TEMP" "$dir/$FINAL_LOG_NAME"
fi
