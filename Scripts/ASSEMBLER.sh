#!/bin/bash


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
threads=${optimal_threads:-1}  # Si es 0, usar 1 hilo
threads=${threads%%.*}  # Eliminar decimales (redondeo hacia abajo)

echo "Usando $threads hilos para el procesamiento"

# Definir rutas del clasificador según el serotipo
classification_available=false
model_path=""
encoder_path=""

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

case "$serotype" in
    "DENV_1")
        fasta_file="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_1.fasta"
        model_path="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV1/dengue_lineage_classifier.joblib.gz"
        encoder_path="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV1/label_encoder.joblib"
        ;;
    "DENV_2")
        fasta_file="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_2.fasta"
        model_path="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV2/dengue_lineage_classifier.joblib.gz"
        encoder_path="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV2/label_encoder.joblib"
        ;;
    "DENV_3")
        fasta_file="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_3.fasta"
        model_path="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV3/dengue_lineage_classifier.joblib.gz"
        encoder_path="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV3/label_encoder.joblib"
        ;;
    "DENV_4")
        fasta_file="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_4.fasta"
        model_path="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV4/dengue_lineage_classifier.joblib.gz"
        encoder_path="$current_dir/CONSENSO_D/Viral-Branch/VB_DENV4/label_encoder.joblib"
        ;;
    "SARS_COV_2")
        fasta_file="$current_dir/CONSENSO_D/Ref_DENV/REF_NC_045512_SARS_COV_2.fasta"
        ;;
    "RABV")
        fasta_file="$current_dir/CONSENSO_D/Ref_DENV/RABV_Reference.fasta"
        ;;
    "N_RABV")
        fasta_file="$current_dir/CONSENSO_D/Ref_DENV/Lisa_RABV_N.fasta"
        ;;
    *)
        echo "Error: serotipo / Virus no válido"
        echo "

Intentar con las opciones
DENV_1
DENV_2
DENV_3
DENV_4
SARS_COV_2
RABV
N_RABV"
        exit 1
        ;;
esac

# Verificar si la clasificación está disponible para este serotipo
if [[ "$serotype" =~ ^DENV_[1-4]$ ]] && [ -f "$model_path" ] && [ -f "$encoder_path" ] && [ -f "$classifier_script" ]; then
    classification_available=true
    echo "🧬 Clasificación de linajes disponible para $serotype"
else
    if [[ "$serotype" =~ ^DENV_[1-4]$ ]]; then
        echo "ℹ️  Clasificación no disponible para $serotype"
        if [ ! -f "$model_path" ]; then
            echo "   ❌ Modelo no encontrado: $(basename "$model_path")"
        fi
        if [ ! -f "$encoder_path" ]; then
            echo "   ❌ Encoder no encontrado: $(basename "$encoder_path")"
        fi
        if [ ! -f "$classifier_script" ]; then
            echo "   ❌ Script clasificador no encontrado: $(basename "$classifier_script")"
        fi
    fi
fi

echo "-----------------------------"
sleep 0.3
echo "-- --"
sleep 0.3
echo "-- --"
sleep 0.3
echo "-- ENSAMBLANDO --"
sleep 0.3
echo "-- --"
echo "-- $serotype --"
sleep 0.3
echo "-- --"
if [ "$classification_available" = true ]; then
echo "-- CON CLASIFICACIÓN --"
fi
sleep 0.3
echo "-- --"
sleep 0.3
echo "-----------------------------"
sleep 1
echo "-- Instituto De ciencias Sostenibles --"
sleep 0.1
echo "------------------------------------------"

if [ ! -d "$dir" ]; then
    echo "No se puede acceder al directorio: $dir"
    exit 1
fi

cd "$dir"

# Crear directorios para resultados en la raíz del directorio de trabajo
mkdir -p QC_Reports
mkdir -p Lineage_Classification

# Preparar archivos de adaptadores y primers si existen
if [ -f "$current_dir/CONSENSO_D/Adapters/adapters.csv" ]; then
    awk -F',' 'NR>1 && $1!="" && $2!="" {print ">"$1"\n"$2}' "$current_dir/CONSENSO_D/Adapters/adapters.csv" > adapters.fasta
    ADAPTERS_FASTA="adapters.fasta"
else
    ADAPTERS_FASTA=""
fi

# Preparar archivo de primers a partir del argumento
PRIMERS_FASTA=""
primer_file_path=$4

if [ -n "$primer_file_path" ] && [ "$primer_file_path" != "none" ] && [ -f "$primer_file_path" ]; then
    echo "🧬 Usando archivo de primers proporcionado: $primer_file_path"
    PRIMERS_FASTA="$primer_file_path"
    echo "✅ Primers listos para ser usados."
elif [ -n "$primer_file_path" ] && [ "$primer_file_path" != "none" ]; then
    echo "❌ Error: No se pudo encontrar el archivo de primers en la ruta: $primer_file_path"
fi

# Función para buscar y procesar VCFs automáticamente
run_batch_classification() {
    local working_dir="$1"

    if [ "$classification_available" = true ]; then
        echo ""
        echo "🧬 BUSCANDO VCFs PARA CLASIFICACIÓN BATCH..."

        # Definir el nombre del archivo de log para la clasificación
        log_file="$working_dir/Lineage_Classification/classification.log"
        output_file="$working_dir/Lineage_Classification/batch_${serotype}_classifications.csv"

        # Ejecutar el clasificador en modo batch y redirigir la salida al log
        echo "Iniciando clasificación. Logs guardados en: $(basename "$log_file")"

        # CAMBIO CRUCIAL: Usar --search-dir en lugar de --vcf-dir y pasar el directorio de trabajo
        python3 "$classifier_script" \
            --search-dir "$working_dir" \
            --model "$model_path" \
            --encoder "$encoder_path" \
            --serotype "$serotype" \
            --output "$output_file" \
            --pattern "*_normalized.vcf.gz" > "$log_file" 2>&1

        classification_success=$?

        if [ $classification_success -eq 0 ]; then
            echo "✅ Clasificación batch completada exitosamente"
            echo "   📊 Resultados guardados en: $(basename "$output_file")"

            # Mostrar resumen rápido si el archivo existe
            if [ -f "$output_file" ]; then
                echo ""
                echo "📋 RESUMEN RÁPIDO DE CLASIFICACIONES:"
                echo "   Total de muestras procesadas: $(tail -n +2 "$output_file" | wc -l)"

                # Resumen por linajes
                echo "   Distribución por linajes:"
                tail -n +2 "$output_file" | cut -d',' -f3 | sort | uniq -c | sort -nr | head -5 | while read count lineage; do
                    echo "      $lineage: $count muestra(s)"
                done

                # Resumen por confianza
                echo "   Distribución por confianza:"
                tail -n +2 "$output_file" | cut -d',' -f4 | sort | uniq -c | sort -nr | while read count confidence; do
                    echo "      $confidence: $count muestra(s)"
                done
            fi
        else
            echo "❌ Error en la clasificación batch"
            echo "   🔍 Para más detalles, revise el log: $log_file"
        fi

        return $classification_success
    fi

    return 1
}

if [ "$sequence_type" == "NANO" ]; then
    # Procesamiento para secuenciación Nanopore (código original)
    for file in $(ls -d */); do
        cd "$file" || continue
        sample_name=${PWD##*/}

        # Procesar Nanopore
        if ls *.fastq.gz 1> /dev/null 2>&1; then
            echo "Procesando archivo : $(basename "$file")"
            printf '%s\n' "-----------------------------------------------------"
            echo "         $(basename "$file")             "
            printf '%s\n' "       Iniciando Proceso de ensamblaje         "
            printf '%s\n' "-----------------------------------------------------"
            printf '%s\n' " "

            printf '%s\n' "Descomprimiendo..."
            cat $(ls *.fastq.gz) > "${sample_name}.fastq.gz"
            gzip -df "${sample_name}.fastq.gz"
            printf '%s\n' "Descompresión realizada con éxito"
            muestra=$(ls *.fastq)
        elif ls *.fastq 1> /dev/null 2>&1; then
            echo "Procesando archivo : $(basename "$file")"
            echo "----------------------------------------------------------------------------"
            echo "⚠  Archivos .fastq ya descomprimidos en la carpeta \"$(basename "$file")\". Continuando y sobrescribiendo...  ⚠"
            echo "----------------------------------------------------------------------------"
            muestra=$(ls *.fastq)
        else
            printf '%s\n' " "
            printf '%s\n' "-----------------------------------------------------"
            echo "         $(basename "$file")             "
            printf '%s\n' "                             "
            printf '%s\n' "-----------------------------------------------------"
            printf '%s\n' " "
            echo "  "
            echo "⚠  No se encontraron archivos fastq.gz ni .fastq en la carpeta \"$(basename "$file")\"  ⚠"
            echo "  "
            cd ..
            continue
        fi

        # Remoción de primers con cutadapt si están disponibles
        if [ -n "$PRIMERS_FASTA" ] && [ -f "$muestra" ]; then
            echo "🧬 Recortando primers con cutadapt..."
            mv "$muestra" "${muestra}.original"
            cutadapt -g file:"$PRIMERS_FASTA" -a file:"$PRIMERS_FASTA" -o "$muestra" "${muestra}.original" --minimum-length 50 -j $threads
            rm "${muestra}.original"
        fi

        # Mapeo para Nanopore con multihilo
        printf '%s\n' "  "
        printf '%s\n' "     Mapeando contra referencia con $threads hilos...       "
        minimap2 -t $threads -a "$fasta_file" "$muestra" > Output.aln.sam

        # Proceso común para ambos tipos de datos
        samtools view -@ $threads -S -b Output.aln.sam > Output.aln.bam
        samtools sort -@ $threads Output.aln.bam -o Output.sorted.bam
        samtools index Output.sorted.bam

        printf '%s\n' "      "
        printf '%s\n' "Verificando calidad del mapeo..."
        samtools flagstat Output.sorted.bam > "${sample_name}.flagstat"
        samtools stats Output.sorted.bam > "${sample_name}.stats"

        # Generar gráficos de calidad con plot-bamstats
        mkdir -p "${sample_name}_plots"
        plot-bamstats -p "${sample_name}_plots/" "${sample_name}.stats"

        # Generar informe HTML de calidad con MultiQC
        multiqc -f -o "${sample_name}_qc_report" "${sample_name}.stats" "${sample_name}.flagstat"

        # Generar archivo de cobertura
        samtools depth Output.sorted.bam > "${sample_name}.coverage"
        awk '{sum+=$3; count++} END {print "Cobertura promedio:", sum/count}' "${sample_name}.coverage"

        # Generar gráfico de cobertura con gnuplot
        echo "set terminal png size 1200,600
         set output '${sample_name}_coverage_plot.png'
         set title 'Cobertura de secuenciación para ${sample_name}'
         set xlabel 'Posición en el genoma'
         set ylabel 'Profundidad de cobertura'
         set grid
         plot '${sample_name}.coverage' using 2:3 with lines linecolor rgb 'blue' title 'Cobertura'" > coverage_plot.gnuplot

        gnuplot coverage_plot.gnuplot

        printf '%s\n' "      "
        printf '%s\n' "::::::::::CREANDO SECUENCIA CONSENSO::::::::::::"
        # Generar archivo VCF y secuencia consenso con multihilo
        bcftools mpileup -Ou -f "$fasta_file" Output.sorted.bam | bcftools call -c -Oz -o "${sample_name}_calls.vcf.gz"
        bcftools norm -f "$fasta_file" "${sample_name}_calls.vcf.gz" -Oz -o "${sample_name}_normalized.vcf.gz"
        bcftools view -i 'QUAL>8' "${sample_name}_normalized.vcf.gz" | vcfutils.pl vcf2fq > SAMPLE_cns.fastq
        seqtk seq -aQ64 -q08 -n N SAMPLE_cns.fastq > SAMPLE_cns.fasta
        echo ">${sample_name}" > "${sample_name}.fasta"
        tail -n +2 SAMPLE_cns.fasta >> "${sample_name}.fasta"

        # Copiar informes de calidad al directorio principal
        cp -r "${sample_name}_plots" "../QC_Reports/${sample_name}_plots" 2>/dev/null || true
        cp -r "${sample_name}_qc_report" "../QC_Reports/${sample_name}_qc_report" 2>/dev/null || true
        cp "${sample_name}_coverage_plot.png" "../QC_Reports/${sample_name}_coverage_plot.png" 2>/dev/null || true

        # Limpiar archivos temporales pero mantener los importantes
        rm SAMPLE_cns.fastq SAMPLE_cns.fasta coverage_plot.gnuplot
        rm Output.aln.sam Output.aln.bam

        printf '%s\n' "----------------------------------------"
        printf '%s\n' "        PROCESO TERMINADO          "
        printf '%s\n' "----------------------------------------"

        cd ..
    done

elif [ "$sequence_type" == "ILLUMINA" ]; then
    # Procesamiento para secuenciación Illumina con mejoras del Script 1

    # 1. Indexar el genoma de referencia (sólo se hace una vez)
    bwa index "$fasta_file"

    # Usar find para procesar cada directorio
    while IFS= read -r -d '' sample_dir; do
        cd "$sample_dir" || continue
        R1_file=$(find . -maxdepth 1 -name '*_R1_001.fastq.gz' -print -quit)
        R2_file=$(find . -maxdepth 1 -name '*_R2_001.fastq.gz' -print -quit)

        if [ -n "$R1_file" ] && [ -n "$R2_file" ]; then
            sample_name=$(basename "$R1_file" "_R1_001.fastq.gz")
            echo "Procesando archivo : $sample_name"
            echo "Mapeando con Illumina para la muestra $sample_name usando $threads hilos..."

            # Generar informe de calidad de lecturas crudas con FastQC
            mkdir -p fastqc_results
            fastqc -t $threads -o fastqc_results "$R1_file" "$R2_file"

            # Paso de limpieza: Trimmomatic con adaptadores si están disponibles
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

            # Remoción de primers con cutadapt si están disponibles
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

            # FastQC post-trimming
            fastqc -t $threads -o fastqc_results "$R1_for_mapping" "$R2_for_mapping"

            # 2. Mapeo con BWA-MEM con multihilo
            bwa mem -t $threads "$fasta_file" "$R1_for_mapping" "$R2_for_mapping" > Output_${sample_name}.aln.sam

            # Proceso común para ambos tipos de datos con multihilo
            samtools view -@ $threads -S -b Output_${sample_name}.aln.sam > Output_${sample_name}.aln.bam
            samtools sort -@ $threads Output_${sample_name}.aln.bam -o Output_${sample_name}.sorted.bam
            samtools index Output_${sample_name}.sorted.bam

            # CRUCIAL: Verificar si el archivo BAM no está vacío antes de continuar.
            if [ -s "Output_${sample_name}.sorted.bam" ]; then
                printf '%s\n' "      "
                printf '%s\n' "Verificando calidad del mapeo..."
                samtools flagstat Output_${sample_name}.sorted.bam > "${sample_name}.flagstat"
                samtools stats Output_${sample_name}.sorted.bam > "${sample_name}.stats"

                # Generar gráficos de calidad con plot-bamstats (igual que en Nanopore)
                mkdir -p "${sample_name}_plots"
                plot-bamstats -p "${sample_name}_plots/" "${sample_name}.stats"

                # Generar informe HTML de calidad con MultiQC
                multiqc -f -o "${sample_name}_qc_report" "${sample_name}.stats" "${sample_name}.flagstat" fastqc_results

                # Generar archivo de cobertura
                samtools depth Output_${sample_name}.sorted.bam > "${sample_name}.coverage"

                # Generar gráfico de cobertura con gnuplot (igual que en Nanopore)
                echo "set terminal png size 1200,600
                 set output '${sample_name}_coverage_plot.png'
                 set title 'Cobertura de secuenciación para ${sample_name}'
                 set xlabel 'Posición en el genoma'
                 set ylabel 'Profundidad de cobertura'
                 set grid
                 plot '${sample_name}.coverage' using 2:3 with lines linecolor rgb 'blue' title 'Cobertura'" > coverage_plot.gnuplot

                gnuplot coverage_plot.gnuplot

                printf '%s\n' "      "
                printf '%s\n' "::::::::::CREANDO SECUENCIA CONSENSO::::::::::::"
                # Generar archivo VCF y secuencia consenso con multihilo - MANTENIENDO QUAL>30 del Script 2
                bcftools mpileup -Ou -f "$fasta_file" Output_${sample_name}.sorted.bam | bcftools call -c -Oz -o "${sample_name}_calls.vcf.gz"
                bcftools norm -f "$fasta_file" "${sample_name}_calls.vcf.gz" -Oz -o "${sample_name}_normalized.vcf.gz"
                #bcftools view -i 'QUAL>30' "${sample_name}_normalized.vcf.gz" | vcfutils.pl vcf2fq > SAMPLE_${sample_name}_cns.fastq
                bcftools view "${sample_name}_normalized.vcf.gz" | vcfutils.pl vcf2fq > SAMPLE_${sample_name}_cns.fastq
                #seqtk seq -aQ64 -q30 -n N SAMPLE_${sample_name}_cns.fastq > SAMPLE_${sample_name}_cns.fasta
                seqtk seq -aQ64 -n N SAMPLE_${sample_name}_cns.fastq > SAMPLE_${sample_name}_cns.fasta
                echo ">${sample_name}" > "${sample_name}.fasta"
                tail -n +2 SAMPLE_${sample_name}_cns.fasta >> "${sample_name}.fasta"


                # Copiar informes de calidad al directorio principal
                cp -r "${sample_name}_plots" "../QC_Reports/${sample_name}_plots" 2>/dev/null || true
                cp -r "${sample_name}_qc_report" "../QC_Reports/${sample_name}_qc_report" 2>/dev/null || true
                cp "${sample_name}_coverage_plot.png" "../QC_Reports/${sample_name}_coverage_plot.png" 2>/dev/null || true
            else
                echo "⚠ WARNING: Output BAM file is empty for sample $sample_name. Skipping consensus generation."
            fi

            # Limpiar archivos temporales pero mantener los importantes
            rm Output_${sample_name}.aln.sam Output_${sample_name}.aln.bam
            rm SAMPLE_${sample_name}_cns.fastq SAMPLE_${sample_name}_cns.fasta coverage_plot.gnuplot 2>/dev/null || true
            rm "$R1_unpaired" "$R2_unpaired" 2>/dev/null || true

        else
            echo "No se encontraron archivos R1 y R2 adecuados en $sample_dir"
        fi

        cd ..
    done < <(find . -type d -mindepth 1 -maxdepth 1 -print0)

else
    echo "Tipo de secuenciación no válido. Use 'NANO' o 'ILLUMINA'."
    exit 1
fi

# ==================== PROCESAMIENTO AUTOMÁTICO DE CLASIFICACIÓN ====================
# Ejecutar clasificación batch automáticamente después del ensamblaje
run_batch_classification "$dir"

# Generar archivo multifasta con todos los consensos
printf '%s\n' " "
printf '%s\n' " "
cat $(find . -name "*.fasta") > "all_consensus_$serotype.fasta"

# Generar informe final combinado con MultiQC
multiqc -f -o "QC_Reports/Final_QC_Report_$serotype" QC_Reports/

# Limpiar archivos temporales
rm -f adapters.fasta primers.fasta 2>/dev/null

#Informe
echo "*****************************************"   # ✅ Correcto
echo "****       GENERANDO INFORME         ****"   # ✅ Correcto
echo "*****************************************"   # ✅ Correcto

python3 "$current_dir/CONSENSO_D/Scripts/quality_report.py" "$dir"
if [ $? -eq 0 ]; then
    echo "✅ Informe guardado en: $(pwd)/reporte_calidad.html"
else
    echo "❌ Error generando el informe"
fi


echo " ____ ____ ___ __ ___ _____ ___
| | \ / \ / ] / ]/ / / \
| o ) D ) | / / / [( _ | |
| /| /| O |/ / | ]_ || O |
| | | | / _ | [ / \ || |
| | | . \ \ || |\ || |
|| ||_|_/ _||_| _| ___/

| || || \ / || | | || | / || \ / \
| | | | | _ || o || | | | |/ || o || \ | |
| |_ | | | | || || |___ | | | || || D || O |
| _] | | | | || _ || | | | | / || _ || || |
| | | | | | || | || | | | | || |
|| || |
|| |||||||||||||_____||||_|||||___| _/

"
echo "Se ha generado un archivo multi-fasta llamado: "
echo " --> all_consensus_$serotype.fasta"
echo " en la ubicación: $dir "
echo ""
echo "Los informes de calidad están disponibles en: $dir/QC_Reports/"
echo "Informe final combinado: $dir/QC_Reports/Final_QC_Report_$serotype"

if [ "$classification_available" = true ]; then
    echo ""
    echo "🧬 Los resultados de clasificación están disponibles en: $dir/Lineage_Classification/"
    if [ -f "$dir/Lineage_Classification/batch_${serotype}_classifications.csv" ]; then
        echo "📊 Clasificación automática: $dir/Lineage_Classification/batch_${serotype}_classifications.csv"
        echo "📄 Log de clasificación: $dir/Lineage_Classification/classification.log"
    fi
fi


# Guardar los logs en el log file
# Definir el nombre del archivo de log final
FINAL_LOG_NAME="CONSENSO_log.txt"

# Guardar los logs en el log file
if [ -f "$LOG_FILE_TEMP" ]; then
    mv "$LOG_FILE_TEMP" "$dir/$FINAL_LOG_NAME"
    echo "✅ Log movido a: $dir/$FINAL_LOG_NAME"


echo "✅ El pipeline ha finalizado exitosamente."
echo "Los registros completos del proceso se guardaron en:"
echo "$dir/CONSENSO_log.txt"
