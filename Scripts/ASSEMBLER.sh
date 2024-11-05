#!/bin/bash
#Salir a la carpeta homa
cd

# Obtener el tipo de secuenciación
sequence_type=$1

# Obtener el serotipo a ensamblar y la ruta de la secuencia de referencia
serotype=$2

# Obtener la ruta del directorio actual
current_dir=$PWD

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
else
    echo "Error: serotipo / Virus no válido"
    echo "
    Intentar con las opciones
              DENV_1
              DENV_2
              DENV_3
              DENV_4
              SARS_COV_2
              RABV "
    exit 1
fi

echo "--*--*--*--*--*--*--*--*--*--*--*--*--*---"
sleep 0.3
echo "--                                      --"
sleep 0.3
echo "--                                      --"
sleep 0.3
echo "--             ENSAMBLANDO              --"
sleep 0.3
echo "--                                      --"
sleep 0.3
echo "--             $serotype                --"
sleep 0.3
echo "--                                      --"
sleep 0.3
echo "--                                      --"
sleep 0.3
echo "--*--*--*--*--*--*--*--*--*--*--*--*--*---"
sleep 1
echo "--  Instituto De ciencias Sostenibles   --"
sleep 0.1
echo "------------------------------------------"

dir=$3
if [ ! -d "$dir" ]; then
    echo "No se puede acceder al directorio: $dir"
    exit 1
fi

cd "$dir"

if [ "$sequence_type" == "NANO" ]; then
    # El bloque para NANO permanece igual

    for file in $(ls -d */); do
        cd "$file" || continue
        sample_name=${PWD##*/}
        
        # Procesar Nanopore
        if ls *.fastq.gz 1> /dev/null 2>&1; then
            printf '%s\n' " "
            printf '%s\n' " "
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

        # Mapeo para Nanopore
        printf '%s\n' "  "
        printf '%s\n' "      Mapeando contra referencia...      "
        minimap2 -a "$fasta_file" "$muestra" > Output.aln.sam

        # Proceso común para ambos tipos de datos
        samtools view -S -b Output.aln.sam > Output.aln.bam
        samtools sort Output.aln.bam -o Output.sorted.bam
        samtools index Output.sorted.bam

        printf '%s\n' "     "
        printf '%s\n' "Verificando calidad del mapeo..."
        samtools flagstat Output.sorted.bam > "${sample_name}.flagstat"
        samtools stats Output.sorted.bam > "${sample_name}.stats"
        plot-bamstats -p "${sample_name}_plots/" "${sample_name}.stats"

        # Generar archivo de cobertura
        samtools depth Output.sorted.bam > "${sample_name}.coverage"

        printf '%s\n' "     "
        printf '%s\n' "::::::::::CREANDO SECUENCIA CONSENSO::::::::::::"
        # Generar archivo VCF y secuencia consenso
        bcftools mpileup -Ou -f "$fasta_file" Output.sorted.bam | bcftools call -c -Oz -o calls.vcf.gz
        bcftools norm -f "$fasta_file" calls.vcf.gz -Oz -o normalized.vcf.gz
        bcftools view -i 'QUAL>20' normalized.vcf.gz | vcfutils.pl vcf2fq > SAMPLE_cns.fastq
        seqtk seq -aQ64 -q20 -n N SAMPLE_cns.fastq > SAMPLE_cns.fasta
        echo ">${sample_name}" > "${sample_name}.fasta"
        tail -n +2 SAMPLE_cns.fasta >> "${sample_name}.fasta"
        rm SAMPLE_cns.fastq SAMPLE_cns.fasta
        rm calls.vcf.gz normalized.vcf.gz Output.aln.sam Output.aln.bam Output.sorted.bam Output.sorted.bam.bai

        printf '%s\n' "----------------------------------------"
        printf '%s\n' "           PROCESO TERMINADO            "
        printf '%s\n' "----------------------------------------"
        cd ..
    done

elif [ "$sequence_type" == "ILLUMINA" ]; then
    # Procesar archivos Illumina en subdirectorios
    while IFS= read -r -d '' sample_dir; do
        cd "$sample_dir" || continue
        R1_file=$(find . -maxdepth 1 -name '*_R1_001.fastq.gz' -print -quit)
        R2_file=$(find . -maxdepth 1 -name '*_R2_001.fastq.gz' -print -quit)

        if [ -n "$R1_file" ] && [ -n "$R2_file" ]; then
            sample_name=$(basename "$R1_file" "_R1_001.fastq.gz")
            echo "Mapeando con Illumina para la muestra $sample_name..."
            minimap2 -ax sr "$fasta_file" "$R1_file" "$R2_file" > Output_${sample_name}.aln.sam

            # Proceso común para ambos tipos de datos
            samtools view -S -b Output_${sample_name}.aln.sam > Output_${sample_name}.aln.bam
            samtools sort Output_${sample_name}.aln.bam -o Output_${sample_name}.sorted.bam
            samtools index Output_${sample_name}.sorted.bam

            printf '%s\n' "     "
            printf '%s\n' "Verificando calidad del mapeo..."
            samtools flagstat Output_${sample_name}.sorted.bam > "${sample_name}.flagstat"
            samtools stats Output_${sample_name}.sorted.bam > "${sample_name}.stats"
            plot-bamstats -p "${sample_name}_plots/" "${sample_name}.stats"

            # Generar archivo de cobertura
            samtools depth Output_${sample_name}.sorted.bam > "${sample_name}.coverage"

            printf '%s\n' "     "
            printf '%s\n' "::::::::::CREANDO SECUENCIA CONSENSO::::::::::::"
            # Generar archivo VCF y secuencia consenso
            bcftools mpileup -Ou -f "$fasta_file" Output_${sample_name}.sorted.bam | bcftools call -c -Oz -o calls_${sample_name}.vcf.gz
            bcftools norm -f "$fasta_file" calls_${sample_name}.vcf.gz -Oz -o normalized_${sample_name}.vcf.gz
            bcftools view -i 'QUAL>20' normalized_${sample_name}.vcf.gz | vcfutils.pl vcf2fq > SAMPLE_${sample_name}_cns.fastq
            seqtk seq -aQ64 -q20 -n N SAMPLE_${sample_name}_cns.fastq > SAMPLE_${sample_name}_cns.fasta
            echo ">${sample_name}" > "${sample_name}.fasta"
            tail -n +2 SAMPLE_${sample_name}_cns.fasta >> "${sample_name}.fasta"
            rm SAMPLE_${sample_name}_cns.fastq SAMPLE_${sample_name}_cns.fasta
            rm calls_${sample_name}.vcf.gz normalized_${sample_name}.vcf.gz Output_${sample_name}.aln.sam Output_${sample_name}.aln.bam Output_${sample_name}.sorted.bam Output_${sample_name}.sorted.bam.bai
        else
            echo "No se encontraron archivos R1 y R2 adecuados en $sample_dir"
        fi

        cd ..
    done < <(find . -type d -mindepth 1 -maxdepth 1 -print0)

else
    echo "Tipo de secuenciación no válido. Use 'NANO' o 'ILLUMINA'."
    exit 1
fi

# ... (el resto del código final permanece igual)

printf '%s\n' "     "
printf '%s\n' "     "
cat $(find . -name "*.fasta") > "all_consensus_$serotype.fasta"
echo "	 ____  ____   ___     __    ___  _____  ___                         
	|    \|    \ /   \   /  ]  /  _]/ ___/ /   \                        
	|  o  )  D  )     | /  /  /  [_(   \_ |     |                       
	|   _/|    /|  O  |/  /  |    _]\__  ||  O  |                       
	|  |  |    \|     /   \_ |   [_ /  \ ||     |                       
	|  |  |  .  \     \     ||     |\    ||     |                       
	|__|  |__|\_|\___/ \____||_____| \___| \___/                        
                                                                    
 _____  ____  ____    ____  _      ____  _____   ____  ___     ___  
|     ||    ||    \  /    || |    |    ||     | /    ||   \   /   \ 
|   __| |  | |  _  ||  o  || |     |  | |__/  ||  o  ||    \ |     |
|  |_   |  | |  |  ||     || |___  |  | |   __||     ||  D  ||  O  |
|   _]  |  | |  |  ||  _  ||     | |  | |  /  ||  _  ||     ||     |
|  |    |  | |  |  ||  |  ||     | |  | |     ||  |  ||     ||     |
|__|   |____||__|__||__|__||_____||____||_____||__|__||_____| \___/ 
                                                                   
"
echo "Se ha generado un archivo multi-fasta llamado: "
echo " -->    *all_consensus_$serotype.fasta*"
echo " en la ubicación: $dir "
