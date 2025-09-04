from Bio import SeqIO

def dna_to_amino_v2(cadena, frame=0):
    codones_dict = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TAT": "Y", "TAC": "Y", "TAA": "-", "TAG": "-",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "-", "TGG": "W",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }

    def traducir_cadena(seq, frame):
        if seq is None or len(seq) < 1:
            return "-"  # Manejar NA o cadenas vacías

        aminoacidos = []
        for i in range(frame, len(seq) - 2, 3):
            triplete = seq[i:i + 3]
            aa = codones_dict.get(triplete, "-")  # Obtener el aminoácido o "-" si no se encuentra
            aminoacidos.append(aa)

        return "".join(aminoacidos)

    return traducir_cadena(cadena[frame:], frame)

def translate_fasta_with_cds_v2(fasta_file_path, cds_start, cds_end):
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        # Extraer solo la región codificante
        dna_sequence = str(record.seq[cds_start-1:cds_end])  # Ajustar a base 0

        # Usar la función dna_to_amino_v2 para traducir la secuencia
        protein_sequence = dna_to_amino_v2(dna_sequence, 0)

        print(f"ID de la secuencia: {record.id}")
        print(f"Descripción: {record.description}")
        print(f"Secuencia de proteínas: {protein_sequence}\n")

# Usar la función con la región codificante conocida
fasta_file_path = '/home/ics2/CONSENSO_D/Ref_DENV/Reference_DV_3.fasta'
cds_start = 95
cds_end = 10267
translate_fasta_with_cds_v2(fasta_file_path, cds_start, cds_end)




"""
Documentación para Interpretación de Anotaciones GenBank
Descripción General

Los archivos GenBank contienen información detallada sobre las secuencias de nucleótidos de genes, genes completos, y regiones genómicas, junto con anotaciones que describen características funcionales, como regiones codificantes, estructuras secundarias y sitios regulatorios. Estas anotaciones son esenciales para el análisis bioinformático y el procesamiento de datos de secuencias.
Componentes Clave de las Anotaciones
1. 5'UTR (Región No Traducida 5')

    Descripción: La región 5'UTR abarca desde el inicio de la transcripción hasta justo antes del codón de inicio de la traducción. Es no codificante y no se traduce en proteínas.
    Función: Regula la traducción, estabilidad y localización del ARN mensajero (ARNm).
    Uso en Análisis: No se traduce. La traducción de proteínas comienza después de esta región.

2. Stem-loops (Lazos de Tallo) y Sitios Regulatorios

    Stem-loop A (SLA): Ubicado en la región 2..69, es una estructura secundaria que puede estar involucrada en la regulación de la transcripción o traducción.

    Oligo U Track Spacer: Un elemento regulatorio en la región 70..78, cuya función específica puede variar. Generalmente está asociado con la modulación de la transcripción o la interacción con proteínas.

    5' upstream AUG region (UAR) y Stem-loop B (SLB): Ubicado en la región 79..94, este elemento puede influir en la iniciación de la traducción y en la eficiencia de la unión del ribosoma.

3. CDS (Región Codificante)

    Descripción: La Región Codificante del ADN (CDS) indica la porción de la secuencia que se traduce en proteínas. En este caso, abarca desde la base 95 hasta la 10267.
    Importancia: Esta es la sección que se traduce para obtener la secuencia de aminoácidos de una proteína.
    Codon_start: Indica el marco de lectura, donde 1 significa que la traducción comienza desde el primer nucleótido de la CDS.
    Uso en Análisis: Asegura que la traducción comience en el lugar correcto y en el marco correcto para producir la secuencia de proteínas deseada.

4. Gene y Locus Information

    Gene: Identifica el gen asociado con la secuencia, en este caso, gene="POLY", que se refiere al gen de la poliproteína.
    Locus Tag y Sinónimos: Proveen identificadores únicos y alternativos para el gen, útiles para referencias cruzadas y búsquedas en bases de datos.
    db_xref: Referencias cruzadas a bases de datos externas que contienen más información sobre el gen.

Proceso de Uso en el Análisis de Secuencias

    Extracción de la Región Codificante: Utiliza las posiciones del CDS (95..10267) para extraer solo la parte de la secuencia que necesitas traducir.

    Traducción de la Secuencia: Usa un marco de lectura que comience en el primer nucleótido del CDS (codon_start=1) para asegurar una traducción precisa.

    Interpretación del Resultado: Compara la secuencia de aminoácidos resultante con las referencias de bases de datos (como GenBank) para validar.

    
"""