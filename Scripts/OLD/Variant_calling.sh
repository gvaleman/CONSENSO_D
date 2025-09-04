#!/bin/bash

# Obtener el serotipo y la ruta del archivo fasta
serotype=$1
input_fasta=$2

# Obtener la ruta del directorio actual
current_dir=$(dirname "$input_fasta")

# Definir la secuencia de referencia según el serotipo especificado
if [ "$serotype" == "DENV_1" ]; then
    reference_fasta="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_1.fasta"
elif [ "$serotype" == "DENV_2" ]; then
    reference_fasta="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_2.fasta"
elif [ "$serotype" == "DENV_3" ]; then
    reference_fasta="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_3.fasta"
elif [ "$serotype" == "DENV_4" ]; then
    reference_fasta="$current_dir/CONSENSO_D/Ref_DENV/Reference_DV_4.fasta"
elif [ "$serotype" == "RABV" ]; then
    reference_fasta="$current_dir/CONSENSO_D/Ref_DENV/RABV_Reference.fasta"
else
    echo "Error: serotipo / Virus no válido"
    echo "
    Intentar con las opciones
              DENV_1
              DENV_2
              DENV_3
              DENV_4
              RABV "
    exit 1
fi

echo "--*--*--*--*--*--*--*--*--*--*--*--*--*---"
echo "--                                      --"
echo "--                                      --"
echo "--             LLAMADO VARIANTES        --"
echo "--                                      --"
echo "--               $serotype              --"
echo "--                                      --"
echo "--                                      --"
echo "--*--*--*--*--*--*--*--*--*--*--*--*--*---"
echo "--  Instituto De ciencias Sostenibles   --"
echo "------------------------------------------"

# Crear el directorio de resultados
result_dir="$current_dir/variant_calling_analysis"
mkdir -p $result_dir/vcfs

# Comprobación del archivo fasta de entrada 
if [ -z "$input_fasta" ]; then
    echo "Error: No se proporcionó un archivo fasta de entrada."
    exit 1
fi

# Procesar cada secuencia en el archivo multiFASTA
seq_count=0
while read line; do
    if [[ $line == ">"* ]]; then
        if [[ $seq_count -gt 0 ]]; then
            # Procesar la secuencia leída
            echo "Mapeando contra la referencia..."
            minimap2 -a $reference_fasta "$current_seq.fasta" > "$result_dir/aln_$seq_count.sam"
            samtools view -S -b "$result_dir/aln_$seq_count.sam" > "$result_dir/aln_$seq_count.bam"
            samtools sort "$result_dir/aln_$seq_count.bam" -o "$result_dir/sorted_$seq_count.bam"
            samtools index "$result_dir/sorted_$seq_count.bam"

            # Generar archivo VCF
            echo "Llamado de variantes para la secuencia $seq_count..."
            bcftools mpileup -Ou -f $reference_fasta "$result_dir/sorted_$seq_count.bam" | bcftools call -mv -Oz -o "$result_dir/vcfs/calls_$seq_count.vcf.gz"
        fi
        seq_count=$((seq_count + 1))
        current_seq="$result_dir/sequence_$seq_count"
        echo $line > "$current_seq.fasta"
    else
        echo $line >> "$current_seq.fasta"
    fi
done < $input_fasta

# Asegúrate de procesar la última secuencia
if [[ -n "$current_seq" ]]; then
    echo "Mapeando contra la referencia..."
    minimap2 -a $reference_fasta "$current_seq.fasta" > "$result_dir/aln_$seq_count.sam"
    samtools view -S -b "$result_dir/aln_$seq_count.sam" > "$result_dir/aln_$seq_count.bam"
    samtools sort "$result_dir/aln_$seq_count.bam" -o "$result_dir/sorted_$seq_count.bam"
    samtools index "$result_dir/sorted_$seq_count.bam"

    # Generar archivo VCF
    echo "Llamado de variantes para la secuencia $seq_count..."
    bcftools mpileup -Ou -f $reference_fasta "$result_dir/sorted_$seq_count.bam" | bcftools call -mv -Oz -o "$result_dir/vcfs/calls_$seq_count.vcf.gz"
fi

echo "Proceso de llamado de variantes completado. Los archivos VCF se han guardado en la carpeta '$result_dir/vcfs/'."
