import os
import gzip
import pandas as pd
from Bio import SeqIO
from Mutaciones import dna_to_amino_v2

# Convertimos las posiciones de nucleótidos a posiciones de aminoácidos
def nucleotide_to_amino_region(region):
    return {
        "start": (region["start"] - 1) // 3 + 1,  # Convertir posición inicial
        "end": (region["end"] - 1) // 3 + 1,      # Convertir posición final
        "gene": region["gene"],
        "product": region["product"],
    }

REGIONES = [
    {"start": 95, "end": 436, "gene": "poly_", "product": "ancC:"},
    {"start": 95, "end": 394, "gene": "poly_", "product": "C:"},
    {"start": 437, "end": 934, "gene": "poly_", "product": "prM:"},
    {"start": 437, "end": 709, "gene": "poly_", "product": "pr:"},
    {"start": 710, "end": 934, "gene": "poly_", "product": "M:"},
    {"start": 935, "end": 2413, "gene": "poly_", "product": "E:"},
    {"start": 2414, "end": 3469, "gene": "poly_", "product": "NS1:"},
    {"start": 3470, "end": 4123, "gene": "poly_", "product": "NS2A:"},
    {"start": 4124, "end": 4513, "gene": "poly_", "product": "NS2B:"},
    {"start": 4514, "end": 6370, "gene": "poly_", "product": "NS3:"},
    {"start": 6371, "end": 6751, "gene": "poly_", "product": "NS4A:"},
    {"start": 6752, "end": 6820, "gene": "poly_", "product": "2K:"},
    {"start": 6821, "end": 7564, "gene": "poly_", "product": "NS4B:"},
    {"start": 7565, "end": 10264, "gene": "poly_", "product": "NS5:"},
]

# Convertimos todas las regiones a posiciones de aminoácidos
REGIONES_AA = [nucleotide_to_amino_region(region) for region in REGIONES]

def parse_vcf(vcf_file):
    try:
        with gzip.open(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                yield pos, ref, alt
    except Exception as e:
        print(f"Error al procesar el archivo VCF {vcf_file}: {e}")

def nucleotide_to_amino_position(nucleotide_pos, cds_start):
    return (nucleotide_pos - cds_start) // 3

def get_amino_acid_from_codon(codon, codon_dict):
    return codon_dict.get(codon, '-')

def find_region(amino_position):
    for region in REGIONES_AA:
        if region["start"] <= amino_position <= region["end"]:
            return f"{region['gene']}{region['product']}"
    return "Unknown Region"

def generate_mutations(vcf_file, fasta_file_path, cds_start, codon_dict):
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        # Obtener la secuencia de ADN a partir del inicio del CDS
        dna_sequence = str(record.seq[cds_start-1:])
        protein_sequence = dna_to_amino_v2(dna_sequence, 0)
    
    mutations = set()  # Utilizar un conjunto para evitar duplicados de mutaciones
    for nucleotide_pos, ref, alt in parse_vcf(vcf_file):
        # Calcular la posición de aminoácido relativa al CDS
        aa_position = nucleotide_to_amino_position(nucleotide_pos, cds_start)
        codon_start = aa_position * 3
        codon_index = (nucleotide_pos - cds_start) % 3  # Calcular la posición dentro del codón
        
        # Imprimir valores para depuración
        print(f"Procesando pos: {nucleotide_pos}, ref: {ref}, alt: {alt}, aa_position: {aa_position}, codon_start: {codon_start}, codon_index: {codon_index}")
        
        # Asegúrate que estamos dentro del rango codificante
        if 0 <= codon_start < len(dna_sequence) - 2:
            original_codon = dna_sequence[codon_start:codon_start+3]
            
            if original_codon[codon_index] == ref:
                mutated_codon = (original_codon[:codon_index] +
                                 alt +
                                 original_codon[codon_index + 1:])
                
                # Obtener el aminoácido mutado
                mutated_aa = get_amino_acid_from_codon(mutated_codon, codon_dict)
                original_aa = protein_sequence[aa_position]
                
                if original_aa != mutated_aa:
                    # Cambiar el formato a AA_originalPosiciónAA_nuevo
                    mutation_str = f"{original_aa}{aa_position+1}{mutated_aa}"
                    region_name = find_region(aa_position + 1)
                    formatted_mutation = f"{region_name}{mutation_str}"
                    mutations.add(formatted_mutation)

    return list(mutations)  # Convertir el conjunto de vuelta a lista para devolver

def load_mapping_as_df(mapping_file_path):
    df = pd.read_csv(mapping_file_path, sep=" -> ", header=None, engine="python", names=["original", "vcf"])
    df['original'] = df['original'].str.strip()
    df['vcf'] = df['vcf'].str.strip()
    
    print("Contenido del DataFrame de mapeo:")
    print(df)
    
    return df

def analyze_directory(main_directory_path, fasta_file_path, codons_dict, cds_start):
    mapping_file_path = os.path.join(main_directory_path, 'vcf_mapping.txt')
    vcf_directory_path = os.path.join(main_directory_path, "vcfs")
    vcf_mapping_df = load_mapping_as_df(mapping_file_path)
    
    data = []
    for file in os.listdir(vcf_directory_path):
        if file.endswith(".vcf.gz"):
            vcf_file_path = os.path.join(vcf_directory_path, file)
            vcf_name_with_ext = f">{file}".lstrip('>')
            mutations = generate_mutations(vcf_file_path, fasta_file_path, cds_start, codons_dict)
            
            # Imprimir información para depuración
            print(f"Archivo VCF procesado: {file}")
            print(f"Mutaciones generadas: {mutations}")
            
            data.append({"vcf": f">{vcf_name_with_ext}", "mutaciones": ";".join(mutations) if mutations else "No mutations found"})

    mutations_df = pd.DataFrame(data)
    print("Contenido del DataFrame de mutaciones:")
    print(mutations_df)
    
    result_df = pd.merge(vcf_mapping_df, mutations_df, on='vcf', how='left')
    
    # Cambiar el nombre de la columna 'original' a 'nombre' y eliminar el símbolo '>'
    result_df.rename(columns={'original': 'nombre'}, inplace=True)
    result_df['nombre'] = result_df['nombre'].str.lstrip('>')
    
    result_df[['nombre', 'mutaciones']].to_csv("mutaciones_resultado.csv", index=False)
    print("Resultados exportados a mutaciones_resultado.csv")
    print(result_df)

# Configuraciones iniciales
fasta_file_path = '/home/ics2/CONSENSO_D/Ref_DENV/Reference_DV_3.fasta'
cds_start = 95
codons_dict = {
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

main_directory_path = '/home/ics2/Dengue-GLUE/variant_calling_analysis'

analyze_directory(main_directory_path, fasta_file_path, codons_dict, cds_start)