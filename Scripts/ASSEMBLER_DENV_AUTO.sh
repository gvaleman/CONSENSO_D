#!/bin/bash

# CONSENSO_D_AUTO.sh - Ensamblaje con Clasificación Automática de Serotipos
# Versión: 1.1-auto
# Basado en CONSENSO_D v0.05-beta

# Salir a la carpeta home para asegurar que las rutas relativas sean correctas
cd

# Obtener los argumentos - CAMBIO: Ya no necesitamos serotipo como parámetro
sequence_type=$1
dir=$2
primer_bed_file=$3

# Obtener la ruta del directorio actual
current_dir=$PWD

# Redirigir la salida a un archivo temporal
LOG_FILE_TEMP="$HOME/proceso_auto_log_temp.txt"
exec > >(tee -a "$LOG_FILE_TEMP") 2>&1

# Calcular número óptimo de hilos
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

# Definir referencias para todos los serotipos
declare -A REFERENCES
REFERENCES["DENV_1"]="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_1.fasta"
REFERENCES["DENV_2"]="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_2.fasta"
REFERENCES["DENV_3"]="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_3.fasta"
REFERENCES["DENV_4"]="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_4.fasta"

declare -A MODELS
MODELS["DENV_1"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV1/dengue_lineage_classifier.joblib.gz"
MODELS["DENV_2"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV2/dengue_lineage_classifier.joblib.gz"
MODELS["DENV_3"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV3/dengue_lineage_classifier.joblib.gz"
MODELS["DENV_4"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV4/dengue_lineage_classifier.joblib.gz"

declare -A ENCODERS
ENCODERS["DENV_1"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV1/label_encoder.joblib"
ENCODERS["DENV_2"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV2/label_encoder.joblib"
ENCODERS["DENV_3"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV3/label_encoder.joblib"
ENCODERS["DENV_4"]="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV4/label_encoder.joblib"

# Colores
CYAN="\\e[36m"
YELLOW="\\e[33m"
GREEN="\\e[32m"
RED="\\e[31m"
WHITE="\\e[97m"
GRAY="\\e[90m"
BOLD="\\e[1m"
RESET="\\e[0m"
BLUE="\\e[34m"
MAGENTA="\\e[35m"

echo -e "${CYAN}====${RESET}"
echo -e "${BOLD}${WHITE}  CONSENSO AUTO ${GRAY}v1.1-auto${RESET}"
echo -e "${CYAN}====${RESET}"
echo ""
echo -e "  ${BOLD}Ensamblaje con Clasificación Automática${RESET}"
echo -e "  Modo: ${WHITE}${BOLD}Detección Automática de Serotipos${RESET}"
echo -e "  Clasificación de linaje: ${GREEN}${BOLD}Automática${RESET}"
echo ""
echo -e "${CYAN}----${RESET}"
echo -e "  ${GRAY}Instituto De Ciencias Sostenibles${RESET}"
echo -e "${CYAN}====${RESET}"

if [ ! -d "$dir" ]; then
    echo "No se puede acceder al directorio: $dir"
    exit 1
fi

cd "$dir"

# Crear directorios para resultados
mkdir -p QC_Reports
mkdir -p Lineage_Classification
mkdir -p Serotype_Classification

# Preparar archivos de adaptadores y primers
if [ -f "$current_dir/CONSENSO_D/Adapters/adapters.csv" ]; then
    awk -F',' 'NR>1 && $1!="" && $2!="" {print ">"$1"\\n"$2}' "$current_dir/CONSENSO_D/Adapters/adapters.csv" > adapters.fasta
    ADAPTERS_FASTA="adapters.fasta"
else
    ADAPTERS_FASTA=""
fi

# Preparar archivo de primers
PRIMERS_FASTA=""
primer_file_path=$3

if [ -n "$primer_file_path" ] && [ "$primer_file_path" != "none" ] && [ -f "$primer_file_path" ]; then
    echo "Usando archivo de primers proporcionado: $primer_file_path"
    PRIMERS_FASTA="$primer_file_path"
    echo "Primers listos para ser usados."
elif [ -n "$primer_file_path" ] && [ "$primer_file_path" != "none" ]; then
    echo "Error: No se pudo encontrar el archivo de primers en la ruta: $primer_file_path"
fi

# NUEVA FUNCIÓN: Clasificación automática de serotipos
classify_serotype_automatically() {
    local sample_name="$1"
    local R1_for_mapping="$2"
    local R2_for_mapping="$3"
    local muestra="$4"  # Para Nanopore

    echo "" >&2
    echo -e "${CYAN}🔍 CLASIFICACIÓN AUTOMÁTICA DE SEROTIPOS${RESET}" >&2
    echo "   Analizando muestra: $sample_name" >&2

    # Crear archivo de métricas
    echo "Serotipo,Cobertura_Porcentaje,Profundidad_Promedio,Lecturas_Mapeadas,Posiciones_Cubiertas" > "${sample_name}_serotype_metrics.csv"

    declare -A coverage_results
    declare -A depth_results
    declare -A mapped_reads

    # Mapear contra cada serotipo
    for serotype in DENV_1 DENV_2 DENV_3 DENV_4; do
        ref_file="${REFERENCES[$serotype]}"

        if [ ! -f "$ref_file" ]; then
            echo "   ⚠️  Referencia no encontrada para $serotype: $ref_file" >&2
            continue
        fi

        echo "   Mapeando contra $serotype..." >&2

        # Mapeo según tipo de secuenciación
        if [ "$sequence_type" == "ILLUMINA" ]; then
            # Indexar referencia si no existe
            if [ ! -f "${ref_file}.bwt" ]; then
                bwa index "$ref_file"
            fi

            (bwa mem -t $threads "$ref_file" "$R1_for_mapping" "$R2_for_mapping" | \
            samtools view -@ $threads -bS -F 4 - | \
            samtools sort -@ $threads - > "${sample_name}_${serotype}_temp.bam")
        else
            # Nanopore
            (minimap2 -t $threads -a "$ref_file" "$muestra" | \
            samtools view -@ $threads -bS -F 4 - | \
            samtools sort -@ $threads - > "${sample_name}_${serotype}_temp.bam")
        fi

        if [ -s "${sample_name}_${serotype}_temp.bam" ]; then
            samtools index "${sample_name}_${serotype}_temp.bam"

            # Calcular métricas
            depth_file="${sample_name}_${serotype}_depth.txt"
            samtools depth "${sample_name}_${serotype}_temp.bam" > "$depth_file"

            if [ -s "$depth_file" ]; then
                # Profundidad promedio
                avg_depth=$(awk '{sum+=$3; count++} END {if(count>0) printf "%.2f", sum/count; else print "0"}' "$depth_file")

                # Longitud del genoma de referencia
                genome_length=$(samtools view -H "${sample_name}_${serotype}_temp.bam" | \
                               grep -E "^@SQ" | awk '{print $3}' | cut -d: -f2)

                # Posiciones cubiertas (profundidad > 0)
                covered_positions=$(wc -l < "$depth_file")

                # Porcentaje de cobertura
                if [ "$genome_length" -gt 0 ]; then
                    coverage_percentage=$(echo "scale=2; $covered_positions * 100 / $genome_length" | bc)
                else
                    coverage_percentage="0"
                fi

                # Lecturas mapeadas
                mapped_count=$(samtools view -c -F 4 "${sample_name}_${serotype}_temp.bam")

                # Guardar resultados
                coverage_results[$serotype]=$coverage_percentage
                depth_results[$serotype]=$avg_depth
                mapped_reads[$serotype]=$mapped_count

                echo "$serotype,$coverage_percentage,$avg_depth,$mapped_count,$covered_positions" >> "${sample_name}_serotype_metrics.csv"

                echo "     $serotype: ${coverage_percentage}% cobertura, profundidad: ${avg_depth}x, lecturas: $mapped_count" >&2
            else
                echo "     $serotype: Sin cobertura" >&2
                echo "$serotype,0,0,0,0" >> "${sample_name}_serotype_metrics.csv"
            fi

            rm -f "$depth_file"
        else
            echo "     $serotype: Sin mapeo" >&2
            echo "$serotype,0,0,0,0" >> "${sample_name}_serotype_metrics.csv"
        fi

        rm -f "${sample_name}_${serotype}_temp.bam" "${sample_name}_${serotype}_temp.bam.bai"
    done

    # Determinar serotipo ganador
    best_serotype=""
    best_coverage=0
    second_best_coverage=0

    for serotype in DENV_1 DENV_2 DENV_3 DENV_4; do
        current_coverage=${coverage_results[$serotype]:-0}
        if (( $(echo "$current_coverage > $best_coverage" | bc -l) )); then
            second_best_coverage=$best_coverage
            best_coverage=$current_coverage
            best_serotype=$serotype
        elif (( $(echo "$current_coverage > $second_best_coverage" | bc -l) )); then
            second_best_coverage=$current_coverage
        fi
    done

    # Evaluar calidad de la clasificación
    classification_confidence="LOW"
    coinfection_detected=false

    if (( $(echo "$best_coverage >= 70" | bc -l) )); then
        if (( $(echo "$second_best_coverage < 30" | bc -l) )); then
            classification_confidence="HIGH"
        elif (( $(echo "$second_best_coverage < 50" | bc -l) )); then
            classification_confidence="MEDIUM"
        else
            classification_confidence="LOW"
            coinfection_detected=true
        fi
    elif (( $(echo "$best_coverage >= 50" | bc -l) )); then
        classification_confidence="MEDIUM"
    fi

    # Detectar co-infecciones
    coinfection_serotypes=""
    if [ "$coinfection_detected" = true ]; then
        for serotype in DENV_1 DENV_2 DENV_3 DENV_4; do
            current_coverage=${coverage_results[$serotype]:-0}
            if (( $(echo "$current_coverage >= 30" | bc -l) )); then
                if [ -n "$coinfection_serotypes" ]; then
                    coinfection_serotypes="${coinfection_serotypes}+${serotype}"
                else
                    coinfection_serotypes="$serotype"
                fi
            fi
        done
    fi

    # Guardar resultados de clasificación
    classification_file="${sample_name}_serotype_classification.txt"
    {
        echo "SAMPLE_NAME: $sample_name"
        echo "DETECTED_SEROTYPE: $best_serotype"
        echo "COVERAGE_PERCENTAGE: $best_coverage"
        echo "CLASSIFICATION_CONFIDENCE: $classification_confidence"
        echo "COINFECTION_DETECTED: $coinfection_detected"
        echo "COINFECTION_SEROTYPES: $coinfection_serotypes"
        echo "DEPTH_AVERAGE: ${depth_results[$best_serotype]:-0}"
        echo "MAPPED_READS: ${mapped_reads[$best_serotype]:-0}"
    } > "$classification_file"

    # Mostrar y devolver resultados
    echo "" >&2
    if [ "$coinfection_detected" = true ] && [ -n "$coinfection_serotypes" ]; then
        echo -e "   ${RED}⚠️  POSIBLE CO-INFECCIÓN: $coinfection_serotypes${RESET}" >&2
        echo -e "   ${YELLOW}📊 Cobertura mejor serotipo: ${best_coverage}% (${best_serotype})${RESET}" >&2
        echo -e "   ${YELLOW}🎯 Confianza: $classification_confidence${RESET}" >&2

        cp "$classification_file" "../Serotype_Classification/${sample_name}_serotype_classification.txt"
        cp "${sample_name}_serotype_metrics.csv" "../Serotype_Classification/${sample_name}_serotype_metrics.csv"

        echo "$coinfection_serotypes" | tr '+' ' '
        return 0
    elif [ -n "$best_serotype" ] && (( $(echo "$best_coverage >= 50" | bc -l) )); then
        echo -e "   ${GREEN}✅ SEROTIPO DETECTADO: ${BOLD}$best_serotype${RESET}" >&2
        echo -e "   ${YELLOW}📊 Cobertura: ${best_coverage}%${RESET}" >&2
        echo -e "   ${YELLOW}🎯 Confianza: $classification_confidence${RESET}" >&2

        cp "$classification_file" "../Serotype_Classification/${sample_name}_serotype_classification.txt"
        cp "${sample_name}_serotype_metrics.csv" "../Serotype_Classification/${sample_name}_serotype_metrics.csv"

        echo "$best_serotype"
        return 0
    else
        echo -e "   ${RED}❌ NO SE PUDO DETERMINAR SEROTIPO${RESET}" >&2
        echo -e "   ${YELLOW}📊 Mejor cobertura: ${best_coverage}% ($best_serotype)${RESET}" >&2

        cp "$classification_file" "../Serotype_Classification/${sample_name}_serotype_classification.txt"
        cp "${sample_name}_serotype_metrics.csv" "../Serotype_Classification/${sample_name}_serotype_metrics.csv"

        echo "UNKNOWN"
        return 1
    fi
}

# Función para clasificación de linajes
run_lineage_classification() {
    local working_dir="$1"
    local detected_serotype="$2"
    local sample_name="$3"

    if [ "$detected_serotype" = "UNKNOWN" ] || [ -z "$detected_serotype" ]; then
        echo "   Saltando clasificación de linajes (serotipo desconocido)"
        return 1
    fi

    model_path="${MODELS[$detected_serotype]}"
    encoder_path="${ENCODERS[$detected_serotype]}"

    if [ -f "$model_path" ] && [ -f "$encoder_path" ] && [ -f "$classifier_script" ]; then
        echo ""
        echo -e "${CYAN}🧬 CLASIFICACIÓN DE LINAJES PARA $detected_serotype${RESET}"
        echo "   Muestra: $sample_name"

        log_file="$working_dir/Lineage_Classification/${sample_name}_${detected_serotype}_classification.log"
        output_file="$working_dir/Lineage_Classification/${sample_name}_${detected_serotype}_classification.csv"

        vcf_file=$(find "$working_dir" -name "${sample_name}_${detected_serotype}_normalized.vcf.gz" -type f | head -1)

        if [ -f "$vcf_file" ]; then
            python3 "$classifier_script" \
                --vcf "$vcf_file" \
                --model "$model_path" \
                --encoder "$encoder_path" \
                --serotype "$detected_serotype" \
                --output "$output_file" > "$log_file" 2>&1

            if [ $? -eq 0 ] && [ -f "$output_file" ]; then
                echo -e "   ${GREEN}✅ Clasificación completada${RESET}"

                if [ -s "$output_file" ]; then
                    lineage=$(tail -n 1 "$output_file" | cut -d',' -f3)
                    confidence=$(tail -n 1 "$output_file" | cut -d',' -f4)
                    echo -e "   ${YELLOW}🏷️  Linaje: $lineage${RESET}"
                    echo -e "   ${YELLOW}🎯 Confianza: $confidence${RESET}"
                fi
                return 0
            else
                echo -e "   ${RED}❌ Error en clasificación de linajes${RESET}"
                return 1
            fi
        else
            echo -e "   ${RED}❌ VCF no encontrado para $sample_name ($detected_serotype)${RESET}"
            return 1
        fi
    else
        echo "   Clasificación de linajes no disponible para $detected_serotype"
        return 1
    fi
}

# PROCESAMIENTO PRINCIPAL
if [ "$sequence_type" == "NANO" ]; then
    echo "Contando muestras..." >&2
    num_samples=$(find . -type d -mindepth 1 -maxdepth 1 | wc -l)
    current_sample=0
    for file in $(ls -d */); do
        current_sample=$((current_sample + 1))
        cd "$file" || continue
        sample_name=${PWD##*/}

        if ls *.fastq.gz 1> /dev/null 2>&1; then
            echo -e "${BOLD}Procesando muestra $current_sample/$num_samples: $sample_name${RESET}"
            echo "----"
            cat $(ls *.fastq.gz) > "${sample_name}.fastq.gz"
            gzip -df "${sample_name}.fastq.gz"
            muestra=$(ls *.fastq)
        elif ls *.fastq 1> /dev/null 2>&1; then
            echo -e "${BOLD}Procesando muestra $current_sample/$num_samples: $sample_name${RESET}"
            echo "----"
            muestra=$(ls *.fastq)
        else
            echo "No se encontraron archivos fastq en $sample_name"
            cd ..
            continue
        fi

        if [ -n "$PRIMERS_FASTA" ] && [ -f "$muestra" ]; then
            echo "Recortando primers con cutadapt..."
            mv "$muestra" "${muestra}.original"
            cutadapt -g file:"$PRIMERS_FASTA" -a file:"$PRIMERS_FASTA" -o "$muestra" "${muestra}.original" --minimum-length 50 -j $threads
            rm "${muestra}.original"
        fi

        detected_serotypes=$(classify_serotype_automatically "$sample_name" "" "" "$muestra")

        if [ "$detected_serotypes" = "UNKNOWN" ]; then
            echo -e "${RED}⚠️  Saltando ensamblaje para $sample_name (serotipo no determinado)${RESET}"
            cd ..
            continue
        fi

        for detected_serotype in $detected_serotypes; do
            fasta_file="${REFERENCES[$detected_serotype]}"

            echo ""
            echo -e "${CYAN}🔧 ENSAMBLAJE ESPECÍFICO PARA $detected_serotype${RESET}"
            echo "   Referencia: $(basename "$fasta_file")"

            echo "   Mapeando contra referencia específica..."
            minimap2 -t $threads -a "$fasta_file" "$muestra" > Output.aln.sam

            samtools view -@ $threads -S -b Output.aln.sam > Output.aln.bam
            samtools sort -@ $threads Output.aln.bam -o Output.sorted.bam
            samtools index Output.sorted.bam

            echo "   Generando métricas de calidad..."
            samtools flagstat Output.sorted.bam > "${sample_name}_${detected_serotype}.flagstat"
            samtools stats Output.sorted.bam > "${sample_name}_${detected_serotype}.stats"

            mkdir -p "${sample_name}_${detected_serotype}_plots"
            plot-bamstats -p "${sample_name}_${detected_serotype}_plots/" "${sample_name}_${detected_serotype}.stats"

            multiqc -f -o "${sample_name}_${detected_serotype}_qc_report" "${sample_name}_${detected_serotype}.stats" "${sample_name}_${detected_serotype}.flagstat"

            samtools depth Output.sorted.bam > "${sample_name}_${detected_serotype}.coverage"

            echo "set terminal png size 1200,600
set output '${sample_name}_${detected_serotype}_coverage_plot.png'
set title 'Cobertura de secuenciacion para ${sample_name} (${detected_serotype})'
set xlabel 'Posicion en el genoma'
set ylabel 'Profundidad de cobertura'
set grid
plot '${sample_name}_${detected_serotype}.coverage' using 2:3 with lines linecolor rgb 'blue' title 'Cobertura'" > coverage_plot.gnuplot

            gnuplot coverage_plot.gnuplot

            echo "   Generando secuencia consenso..."
            bcftools mpileup -Ou -f "$fasta_file" Output.sorted.bam | bcftools call -c -Oz -o "${sample_name}_${detected_serotype}_calls.vcf.gz"
            bcftools norm -f "$fasta_file" "${sample_name}_${detected_serotype}_calls.vcf.gz" -Oz -o "${sample_name}_${detected_serotype}_normalized.vcf.gz"
            bcftools view -i 'QUAL>10' "${sample_name}_${detected_serotype}_normalized.vcf.gz" | vcfutils.pl vcf2fq > SAMPLE_cns.fastq
            seqtk seq -aQ64 -q10 -n N SAMPLE_cns.fastq > SAMPLE_cns.fasta

            echo ">${sample_name}_${detected_serotype}" > "${sample_name}_${detected_serotype}.fasta"
            tail -n +2 SAMPLE_cns.fasta >> "${sample_name}_${detected_serotype}.fasta"

            run_lineage_classification "$PWD/.." "$detected_serotype" "$sample_name"

            cp -r "${sample_name}_${detected_serotype}_plots" "../QC_Reports/"
            cp -r "${sample_name}_${detected_serotype}_qc_report" "../QC_Reports/"
            cp "${sample_name}_${detected_serotype}_coverage_plot.png" "../QC_Reports/"

            rm -f SAMPLE_cns.fastq SAMPLE_cns.fasta coverage_plot.gnuplot
            rm -f Output.aln.sam Output.aln.bam
        done

        echo -e "${GREEN}✅ Procesamiento completado para $sample_name${RESET}"
        echo "----"

        cd ..
    done

elif [ "$sequence_type" == "ILLUMINA" ]; then
    echo "Contando muestras..." >&2
    num_samples=$(find . -type d -mindepth 1 -maxdepth 1 -print0 | tr -d -c '\0' | wc -c)
    current_sample=0
    while IFS= read -r -d '' sample_dir; do
        current_sample=$((current_sample + 1))
        cd "$sample_dir" || continue
        R1_file=$(find . -maxdepth 1 -name '*_R1_001.fastq.gz' -print -quit)
        R2_file=$(find . -maxdepth 1 -name '*_R2_001.fastq.gz' -print -quit)

        if [ -n "$R1_file" ] && [ -n "$R2_file" ]; then
            sample_name=$(basename "$R1_file" "_R1_001.fastq.gz")
            echo -e "${BOLD}Procesando muestra $current_sample/$num_samples: $sample_name${RESET}"
            echo "----"

            mkdir -p fastqc_results
            fastqc -t $threads -o fastqc_results "$R1_file" "$R2_file"

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

            if [ -n "$PRIMERS_FASTA" ]; then
                R1_final="${sample_name}_R1_final_trimmed.fastq.gz"
                R2_final="${sample_name}_R2_final_trimmed.fastq.gz"

                cutadapt -g file:"$PRIMERS_FASTA" -G file:"$PRIMERS_FASTA" \
                    -o "$R1_final" -p "$R2_final" \
                    "$R1_trimmed" "$R2_trimmed" \
                    --minimum-length 50 \
                    -j $threads

                R1_for_mapping="$R1_final"
                R2_for_mapping="$R2_final"
            else
                R1_for_mapping="$R1_trimmed"
                R2_for_mapping="$R2_trimmed"
            fi

            fastqc -t $threads -o fastqc_results "$R1_for_mapping" "$R2_for_mapping"

            detected_serotypes=$(classify_serotype_automatically "$sample_name" "$R1_for_mapping" "$R2_for_mapping" "")

            if [ "$detected_serotypes" = "UNKNOWN" ]; then
                echo -e "${RED}⚠️  Saltando ensamblaje para $sample_name (serotipo no determinado)${RESET}"
                cd ..
                continue
            fi

            for detected_serotype in $detected_serotypes; do
                fasta_file="${REFERENCES[$detected_serotype]}"

                echo ""
                echo -e "${CYAN}🔧 ENSAMBLAJE ESPECÍFICO PARA $detected_serotype${RESET}"
                echo "   Referencia: $(basename "$fasta_file")"

                bwa index "$fasta_file"

                echo "   Mapeando contra referencia específica..."
                bwa mem -t $threads "$fasta_file" "$R1_for_mapping" "$R2_for_mapping" > "Output_${sample_name}.aln.sam"

                samtools view -@ $threads -S -b "Output_${sample_name}.aln.sam" > "Output_${sample_name}.aln.bam"
                samtools sort -@ $threads "Output_${sample_name}.aln.bam" -o "Output_${sample_name}.sorted.bam"
                samtools index "Output_${sample_name}.sorted.bam"

                if [ -s "Output_${sample_name}.sorted.bam" ]; then
                    echo "   Generando métricas de calidad..."
                    samtools flagstat "Output_${sample_name}.sorted.bam" > "${sample_name}_${detected_serotype}.flagstat"
                    samtools stats "Output_${sample_name}.sorted.bam" > "${sample_name}_${detected_serotype}.stats"

                    mkdir -p "${sample_name}_${detected_serotype}_plots"
                    plot-bamstats -p "${sample_name}_${detected_serotype}_plots/" "${sample_name}_${detected_serotype}.stats"

                    multiqc -f -o "${sample_name}_${detected_serotype}_qc_report" "${sample_name}_${detected_serotype}.stats" "${sample_name}_${detected_serotype}.flagstat" fastqc_results

                    samtools depth "Output_${sample_name}.sorted.bam" > "${sample_name}_${detected_serotype}.coverage"

                    echo "set terminal png size 1200,600
set output '${sample_name}_${detected_serotype}_coverage_plot.png'
set title 'Cobertura de secuenciacion para ${sample_name} (${detected_serotype})'
set xlabel 'Posicion en el genoma'
set ylabel 'Profundidad de cobertura'
set grid
plot '${sample_name}_${detected_serotype}.coverage' using 2:3 with lines linecolor rgb 'blue' title 'Cobertura'" > coverage_plot.gnuplot

                    gnuplot coverage_plot.gnuplot

                    echo "   Generando secuencia consenso..."
                    bcftools mpileup -Ou -f "$fasta_file" "Output_${sample_name}.sorted.bam" | \
                    bcftools call -c -Oz -o "${sample_name}_${detected_serotype}_calls.vcf.gz"

                    bcftools norm -f "$fasta_file" "${sample_name}_${detected_serotype}_calls.vcf.gz" \
                    -Oz -o "${sample_name}_${detected_serotype}_normalized.vcf.gz"

                    bcftools view -i 'QUAL>30' "${sample_name}_${detected_serotype}_normalized.vcf.gz" | \
                    vcfutils.pl vcf2fq | \
                    seqtk seq -aQ64 -q30 -n N > "SAMPLE_${sample_name}_cns.fasta"

                    echo ">${sample_name}_${detected_serotype}" > "${sample_name}_${detected_serotype}.fasta"
                    tail -n +2 "SAMPLE_${sample_name}_cns.fasta" >> "${sample_name}_${detected_serotype}.fasta"

                    run_lineage_classification "$PWD/.." "$detected_serotype" "$sample_name"

                    cp -r "${sample_name}_${detected_serotype}_plots" "../QC_Reports/"
                    cp -r "${sample_name}_${detected_serotype}_qc_report" "../QC_Reports/"
                    cp "${sample_name}_${detected_serotype}_coverage_plot.png" "../QC_Reports/"
                else
                    echo -e "${RED}⚠️  BAM vacío para $sample_name ($detected_serotype)${RESET}"
                fi

                rm -f "Output_${sample_name}.aln.sam" "Output_${sample_name}.aln.bam"
                rm -f "SAMPLE_${sample_name}_cns.fasta" coverage_plot.gnuplot
            done

            rm -f "$R1_unpaired" "$R2_unpaired"

            echo -e "${GREEN}✅ Procesamiento completado para $sample_name${RESET}"
            echo "----"
        else
            echo "No se encontraron archivos R1 y R2 en $(basename "$sample_dir")"
        fi

        cd ..
    done < <(find . -type d -mindepth 1 -maxdepth 1 -print0)

else
    echo "Tipo de secuenciación no válido. Use 'NANO' o 'ILLUMINA'."
    exit 1
fi

# GENERAR REPORTES FINALES
echo ""
echo -e "${CYAN}📊 GENERANDO REPORTES FINALES${RESET}"

# Generar archivo multifasta con todos los consensos
echo "Consolidando secuencias consenso..."
fasta_files=$(find . -maxdepth 2 -name "*.fasta" -type f)
if [ -n "$fasta_files" ]; then
    cat $fasta_files > "all_consensus_AUTO.fasta"
else
    echo "No consensus files found to concatenate." >&2
    touch "all_consensus_AUTO.fasta"
fi

# Generar resumen de clasificaciones de serotipos
echo "Generando resumen de clasificaciones..."
{
    echo "Sample_Name,Detected_Serotype,Coverage_Percentage,Classification_Confidence,Coinfection_Detected,Coinfection_Serotypes"
    for classification_file in Serotype_Classification/*_serotype_classification.txt; do
        if [ -f "$classification_file" ]; then
            sample=$(grep "SAMPLE_NAME:" "$classification_file" | cut -d: -f2 | tr -d ' ')
            serotype=$(grep "DETECTED_SEROTYPE:" "$classification_file" | cut -d: -f2 | tr -d ' ')
            coverage=$(grep "COVERAGE_PERCENTAGE:" "$classification_file" | cut -d: -f2 | tr -d ' ')
            confidence=$(grep "CLASSIFICATION_CONFIDENCE:" "$classification_file" | cut -d: -f2 | tr -d ' ')
            coinfection=$(grep "COINFECTION_DETECTED:" "$classification_file" | cut -d: -f2 | tr -d ' ')
            coinfection_types=$(grep "COINFECTION_SEROTYPES:" "$classification_file" | cut -d: -f2 | tr -d ' ')

            echo "$sample,$serotype,$coverage,$confidence,$coinfection,$coinfection_types"
        fi
    done
} > "serotype_classification_summary.csv"

# Generar informe final con MultiQC
multiqc -f -o "QC_Reports/Final_QC_Report_AUTO" QC_Reports/ Serotype_Classification/ Lineage_Classification/

# Consolidar reportes y logs de linajes
echo "Consolidando reportes de linajes..." >&2
FINAL_LINEAGE_CSV="Lineage_Classification/lineage_classifications.csv"
FINAL_LINEAGE_LOG="Lineage_Classification/lineage_classification.log"

# Consolidate CSVs
{
    header_file=$(find Lineage_Classification -name "*_classification.csv" -print -quit)
    if [ -n "$header_file" ] && [ -s "$header_file" ]; then
        head -n 1 "$header_file"
        find Lineage_Classification -name "*_classification.csv" -exec tail -q -n +2 {} +
    else
        echo "Sample,Serotype,Lineage,Confidence"
    fi
} > "$FINAL_LINEAGE_CSV"

# Consolidate Logs
echo "Consolidando logs de linajes..." >&2
cat Lineage_Classification/*_classification.log > "$FINAL_LINEAGE_LOG" 2>/dev/null || true

# Cleanup intermediate files
echo "Limpiando archivos intermedios de linaje..." >&2
rm -f Lineage_Classification/*_classification.csv
rm -f Lineage_Classification/*_classification.log

# Limpiar archivos temporales
rm -f adapters.fasta primers.fasta 2>/dev/null

# Generar informe HTML personalizado
echo "Generando informe interactivo..."
python3 "$current_dir/CONSENSO_D/Scripts/quality_report.py" "$dir" 2>/dev/null || echo "Script de reporte no disponible"

# MENSAJE FINAL
echo ""
echo -e "${GREEN}╔════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║                                        ║${NC}"
echo -e "${GREEN}║${NC}    ${BOLD}${CYAN}✓  P R O C E S O   C O M P L E T O${NC}    ${GREEN}║${NC}"
echo -e "${GREEN}║                                        ║${NC}"
echo -e "${GREEN}╚════════════════════════════════════════╝${NC}"
echo ""
echo -e "    ${CYAN}🧬 CLASIFICACIÓN AUTOMÁTICA${NC}"
echo -e "    ${BLUE}────────────────────────${NC}"
echo -e "    ${YELLOW}Resumen de serotipos:${NC} serotype_classification_summary.csv"
echo -e "    ${YELLOW}Clasificaciones individuales:${NC} Serotype_Classification/"
echo ""
echo -e "    ${CYAN}📊 RESULTADOS${NC}"
echo -e "    ${BLUE}────────────${NC}"
echo -e "    ${YELLOW}Consensos:${NC} all_consensus_AUTO.fasta"
echo -e "    ${YELLOW}Ubicación:${NC} $dir"
echo ""
echo -e "    ${CYAN}📁 CALIDAD${NC}"
echo -e "    ${BLUE}──────────${NC}"
echo -e "    ${YELLOW}Reportes por muestra:${NC} QC_Reports/"
echo -e "    ${YELLOW}Reporte consolidado:${NC} Final_QC_Report_AUTO/"
echo ""
echo -e "    ${CYAN}🧬 LINAJES${NC}"
echo -e "    ${BLUE}─────────${NC}"
echo -e "    ${YELLOW}Reporte de clasificación:${NC} lineage_classifications.csv"
echo -e "    ${YELLOW}Log de clasificación:${NC} lineage_classification.log"
echo -e "    ${YELLOW}Ubicación:${NC} Lineage_Classification/"
echo ""

# Mostrar estadísticas finales
if [ -f "serotype_classification_summary.csv" ]; then
    echo -e "    ${CYAN}📈 ESTADÍSTICAS${NC}"
    echo -e "    ${BLUE}──────────────${NC}"

    total_samples=$(tail -n +2 "serotype_classification_summary.csv" | wc -l)
    echo -e "    ${YELLOW}Total de muestras:${NC} $total_samples"

    if [ $total_samples -gt 0 ]; then
        echo -e "    ${YELLOW}Distribución por serotipo:${NC}"
        tail -n +2 "serotype_classification_summary.csv" | cut -d',' -f2 | sort | uniq -c | sort -nr | while read count serotype; do
            echo -e "      $serotype: $count muestra(s)"
        done

        coinfections=$(tail -n +2 "serotype_classification_summary.csv" | cut -d',' -f5 | grep -c "true" || echo "0")
        if [ $coinfections -gt 0 ]; then
            echo -e "    ${RED}⚠️  Co-infecciones detectadas: $coinfections${NC}"
        fi
    fi
    echo ""
fi

echo -e "    ${CYAN}📝 LOGS:${NC} CONSENSO_AUTO_log.txt"
echo ""
echo -e "    ${GREEN}════${NC}"
echo -e "    ${BOLD}Finalizado:${NC} ${MAGENTA}$(date '+%Y-%m-%d %H:%M:%S')${NC}"
echo -e "    ${GREEN}════${NC}"

# Guardar logs
FINAL_LOG_NAME="CONSENSO_AUTO_log.txt"
if [ -f "$LOG_FILE_TEMP" ]; then
    mv "$LOG_FILE_TEMP" "$dir/$FINAL_LOG_NAME"
    echo -e "${GREEN}✓${NC} Log guardado en: ${YELLOW}$dir/$FINAL_LOG_NAME${NC}"
fi

echo ""
echo -e "${GREEN}✓${NC} Pipeline automático finalizado exitosamente."
echo -e "Revise los archivos de clasificación para detalles específicos."
