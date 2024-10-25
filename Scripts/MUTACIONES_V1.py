import gzip
from Bio import SeqIO
from Mutaciones import dna_to_amino_v2

# Definimos las regiones y productos basados en el diccionario del GenBank
REGIONES = [
    {"start": 95, "end": 436, "gene": "poly", "product": "anchored capsid protein ancC"},
    {"start": 95, "end": 394, "gene": "poly", "product": "capsid protein C"},
    {"start": 437, "end": 934, "gene": "poly", "product": "membrane glycoprotein precursor prM"},
    {"start": 437, "end": 709, "gene": "poly", "product": "protein pr"},
    {"start": 710, "end": 934, "gene": "poly", "product": "membrane glycoprotein M"},
    {"start": 935, "end": 2413, "gene": "poly", "product": "envelope protein E"},
    {"start": 2414, "end": 3469, "gene": "poly", "product": "nonstructural protein NS1"},
    {"start": 3470, "end": 4123, "gene": "poly", "product": "nonstructural protein NS2A"},
    {"start": 4124, "end": 4513, "gene": "poly", "product": "nonstructural protein NS2B"},
    {"start": 4514, "end": 6370, "gene": "poly", "product": "nonstructural protein NS3"},
    {"start": 6371, "end": 6751, "gene": "poly", "product": "nonstructural protein NS4A"},
    {"start": 6752, "end": 6820, "gene": "poly", "product": "protein 2K"},
    {"start": 6821, "end": 7564, "gene": "poly", "product": "nonstructural protein NS4B"},
    {"start": 7565, "end": 10264, "gene": "poly", "product": "RNA-dependent RNA polymerase NS5"},
]

REGIONES = [
    {"start": 95, "end": 436, "gene": "poly_","product": "ancC:"},
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


def parse_vcf(vcf_file):
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            yield pos, ref, alt

def nucleotide_to_amino_position(nucleotide_pos, cds_start):
    return (nucleotide_pos - cds_start) // 3

def get_amino_acid_from_codon(codon, codon_dict):
    return codon_dict.get(codon, '-')

def find_region(amino_position):
    for region in REGIONES:
        if region["start"] <= amino_position <= region["end"]:
            return f"{region['gene']} {region['product']}"
    return "Unknown Region"

def generate_mutations(vcf_file, fasta_file_path, cds_start, codon_dict):
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        dna_sequence = str(record.seq[cds_start-1:])
        protein_sequence = dna_to_amino_v2(dna_sequence, 0)
    
    mutations = []
    for nucleotide_pos, ref, alt in parse_vcf(vcf_file):
        aa_position = nucleotide_to_amino_position(nucleotide_pos, cds_start)
        if aa_position >= 0 and aa_position < len(protein_sequence):
            original_aa = protein_sequence[aa_position]
            
            mutated_codon = (dna_sequence[aa_position*3:aa_position*3+3].replace(ref, alt))
            mutated_aa = get_amino_acid_from_codon(mutated_codon, codon_dict)
            
            if original_aa != mutated_aa:
                mutation_str = f"{original_aa}{aa_position+1}{mutated_aa}"
                region_name = find_region(aa_position + 1)
                # Eliminar espacios al formatear las cadenas
                formatted_mutation = f"{region_name.replace(' ', '')}{mutation_str}"
                mutations.append(formatted_mutation)

    return mutations

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

vcf_file_path = '/home/ics2/Dengue-GLUE/variant_calling_analysis/vcfs/D3_G_NI_336.11A1_MZ008476_2014.vcf.gz'
fasta_file_path = '/home/ics2/CONSENSO_D/Ref_DENV/Reference_DV_3.fasta'
cds_start = 95

# Archivo VCF y Fasta
vcf_file_path = '/home/ics2/Dengue-GLUE/variant_calling_analysis/vcfs/>D3_G_NI_336.11A1_MZ008476_2014.vcf.gz' 
fasta_file_path = '/home/ics2/CONSENSO_D/Ref_DENV/Reference_DV_3.fasta'
cds_start = 95

# Generar mutaciones
mutations = generate_mutations(vcf_file_path, fasta_file_path, cds_start, codons_dict)
for mutation in mutations:
    print(mutation)