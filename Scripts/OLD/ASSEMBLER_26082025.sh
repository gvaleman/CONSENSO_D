#!/bin/bash

# Salir a la carpeta home
cd

# Obtener el tipo de secuenciación
sequence_type=$1

# Obtener el serotipo a ensamblar
serotype=$2

# Obtener la ruta del directorio actual
current_dir=$PWD

# Calcular número óptimo de hilos (80% de los disponibles, mínimo 1)
available_threads=$(nproc)
optimal_threads=$(echo "scale=0; ($available_threads * 0.8)/1" | bc)
threads=${optimal_threads:-1}  # Si es 0, usar 1 hilo
threads=${threads%%.*}  # Eliminar decimales (redondeo hacia abajo)

echo "Usando $threads hilos para el procesamiento"

# Definir la secuencia de referencia según el serotipo especificado
if [ "$serotype" == "DENV_1" ]; then
    fasta_file="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_1.fasta"
elif [ "$serotype" == "DENV_2" ]; then
    fasta_file="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_2.fasta"
elif [ "$serotype" == "DENV_3" ]; then
    fasta_file="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_3.fasta"
elif [ "$serotype" == "DENV_4" ]; then
    fasta_file="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_4.fasta"
elif [ "$serotype" == "SARS_COV_2" ]; then
    fasta_file="$current_dir/CONSENSO_D/Ref_DENV/REF_NC_045512_SARS_COV_2.fasta"
elif [ "$serotype" == "RABV" ]; then
    fasta_file="$current_dir/CONSENSO_D/Ref_DENV/RABV_Reference.fasta"
elif [ "$serotype" == "N_RABV" ]; then
    fasta_file="$current_dir/CONSENSO_D/Ref_DENV/Lisa_RABV_N.fasta"
else
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
sleep 0.3
echo "-- $serotype --"
sleep 0.3
echo "-- --"
sleep 0.3
echo "-- --"
sleep 0.3
echo "-----------------------------"
sleep 1
echo "-- Instituto De ciencias Sostenibles --"
sleep 0.1
echo "------------------------------------------"

dir=$3
if [ ! -d "$dir" ]; then
    echo "No se puede acceder al directorio: $dir"
    exit 1
fi

cd "$dir"

# Crear directorio para resultados QC
mkdir -p QC_Reports

if [ "$sequence_type" == "NANO" ]; then
    # Procesamiento para secuenciación Nanopore
    for file in $(ls -d */); do
        cd "$file" || continue
        sample_name=${PWD##*/}
        
        # Procesar Nanopore
        if ls *.fastq.gz 1> /dev/null 2>&1; then
            echo "Procesando archivo : $(basename "$file")"
            printf '%s\n' "-----------------------------------------------------"
            echo "                $(basename "$file")                  "
            printf '%s\n' "          Iniciando Proceso de ensamblaje            "
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
            echo "                $(basename "$file")                  "
            printf '%s\n' "                                                     "
            printf '%s\n' "-----------------------------------------------------"
            printf '%s\n' " "
            echo "  "
            echo "⚠  No se encontraron archivos fastq.gz ni .fastq en la carpeta \"$(basename "$file")\"  ⚠"
            echo "  "
            cd ..
            continue
        fi

        # Mapeo para Nanopore con multihilo
        printf '%s\n' "  "
        printf '%s\n' "      Mapeando contra referencia con $threads hilos...      "
        minimap2 -t $threads -a "$fasta_file" "$muestra" > Output.aln.sam

        # Proceso común para ambos tipos de datos
        samtools view -@ $threads -S -b Output.aln.sam > Output.aln.bam
        samtools sort -@ $threads Output.aln.bam -o Output.sorted.bam
        samtools index Output.sorted.bam

        printf '%s\n' "     "
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

        printf '%s\n' "     "
        printf '%s\n' "::::::::::CREANDO SECUENCIA CONSENSO::::::::::::"
        # Generar archivo VCF y secuencia consenso con multihilo
        bcftools mpileup -Ou -f "$fasta_file" Output.sorted.bam | bcftools call -c -Oz -o "${sample_name}_calls.vcf.gz"
        bcftools norm -f "$fasta_file" "${sample_name}_calls.vcf.gz" -Oz -o "${sample_name}_normalized.vcf.gz"
        bcftools view -i 'QUAL>8' "${sample_name}_normalized.vcf.gz" | vcfutils.pl vcf2fq > SAMPLE_cns.fastq
        seqtk seq -aQ64 -q08 -n N SAMPLE_cns.fastq > SAMPLE_cns.fasta
        echo ">${sample_name}" > "${sample_name}.fasta"
        tail -n +2 SAMPLE_cns.fasta >> "${sample_name}.fasta"
        
        # Copiar informes de calidad al directorio principal
        cp -r "${sample_name}_plots" "../QC_Reports/${sample_name}_plots"
        cp -r "${sample_name}_qc_report" "../QC_Reports/${sample_name}_qc_report"
        cp "${sample_name}_coverage_plot.png" "../QC_Reports/${sample_name}_coverage_plot.png"
        
        # Limpiar archivos temporales pero mantener los importantes
        rm SAMPLE_cns.fastq SAMPLE_cns.fasta coverage_plot.gnuplot
        rm Output.aln.sam Output.aln.bam

        printf '%s\n' "----------------------------------------"
        printf '%s\n' "           PROCESO TERMINADO            "
        printf '%s\n' "----------------------------------------"
        cd ..
    done

elif [ "$sequence_type" == "ILLUMINA" ]; then
    # Procesamiento para secuenciación Illumina
    
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

            # 2. Mapeo con BWA-MEM con multihilo
            bwa mem -t $threads "$fasta_file" "$R1_file" "$R2_file" > Output_${sample_name}.aln.sam

            # Proceso común para ambos tipos de datos con multihilo
            samtools view -@ $threads -S -b Output_${sample_name}.aln.sam > Output_${sample_name}.aln.bam
            samtools sort -@ $threads Output_${sample_name}.aln.bam -o Output_${sample_name}.sorted.bam
            samtools index Output_${sample_name}.sorted.bam

            # CRUCIAL: Verificar si el archivo BAM no está vacío antes de continuar.
            if [ -s "Output_${sample_name}.sorted.bam" ]; then
                printf '%s\n' "     "
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

                printf '%s\n' "     "
                printf '%s\n' "::::::::::CREANDO SECUENCIA CONSENSO::::::::::::"
                # Generar archivo VCF y secuencia consenso con multihilo
                bcftools mpileup -Ou -f "$fasta_file" Output_${sample_name}.sorted.bam | bcftools call -c -Oz -o "${sample_name}_calls.vcf.gz"
                bcftools norm -f "$fasta_file" "${sample_name}_calls.vcf.gz" -Oz -o "${sample_name}_normalized.vcf.gz"
                bcftools view -i 'QUAL>20' "${sample_name}_normalized.vcf.gz" | vcfutils.pl vcf2fq > SAMPLE_${sample_name}_cns.fastq
                seqtk seq -aQ64 -q20 -n N SAMPLE_${sample_name}_cns.fastq > SAMPLE_${sample_name}_cns.fasta
                echo ">${sample_name}" > "${sample_name}.fasta"
                tail -n +2 SAMPLE_${sample_name}_cns.fasta >> "${sample_name}.fasta"
                
                # Copiar informes de calidad al directorio principal
                cp -r "${sample_name}_plots" "../QC_Reports/${sample_name}_plots"
                cp -r "${sample_name}_qc_report" "../QC_Reports/${sample_name}_qc_report"
                cp "${sample_name}_coverage_plot.png" "../QC_Reports/${sample_name}_coverage_plot.png"
            else
                echo "⚠ WARNING: Output BAM file is empty for sample $sample_name. Skipping consensus generation."
            fi

            # Limpiar archivos temporales pero mantener los importantes
            rm Output_${sample_name}.aln.sam Output_${sample_name}.aln.bam
            rm SAMPLE_${sample_name}_cns.fastq SAMPLE_${sample_name}_cns.fasta coverage_plot.gnuplot

        else
            echo "No se encontraron archivos R1 y R2 adecuados en $sample_dir"
        fi

        cd ..
    done < <(find . -type d -mindepth 1 -maxdepth 1 -print0)

else
    echo "Tipo de secuenciación no válido. Use 'NANO' o 'ILLUMINA'."
    exit 1
fi

# Generar archivo multifasta con todos los consensos
printf '%s\n' " "
printf '%s\n' " "
cat $(find . -name "*.fasta") > "all_consensus_$serotype.fasta"

# Generar informe final combinado con MultiQC
multiqc -f -o "QC_Reports/Final_QC_Report_$serotype" QC_Reports/

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
| | | | | | || | || | | | | || | || || |
|| |||||||||_____||||_|||||___| _/

"
echo "Se ha generado un archivo multi-fasta llamado: "
echo " --> all_consensus_$serotype.fasta"
echo " en la ubicación: $dir "
echo ""
echo "Los informes de calidad están disponibles en: $dir/QC_Reports/"
echo "Informe final combinado: $dir/QC_Reports/Final_QC_Report_$serotype"
