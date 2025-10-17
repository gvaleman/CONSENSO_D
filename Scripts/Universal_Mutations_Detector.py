#!/usr/bin/env python3
"""
Universal Mutations Detector - CONSENSO System
Detecta mutaciones en formato Nextclade para múltiples virus
Compatible con DENV1-4, SARS-CoV-2, RABV

Autor: Sistema CONSENSO
Versión: 1.0
"""

import os
import sys
import gzip
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional
from dataclasses import dataclass
from Bio import SeqIO
import pandas as pd


# ==================== CONFIGURACIÓN DE LOGGING ====================
def setup_logging(output_dir: str, verbose: bool = False):
    """Configura el sistema de logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    log_file = os.path.join(output_dir, "mutations_detection.log")
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


# ==================== MAPAS GENÓMICOS ====================
@dataclass
class GenomicRegion:
    """Región genómica con información de producto"""
    start: int
    end: int
    gene: str
    product: str
    
    def __post_init__(self):
        """Convierte posiciones de nucleótidos a aminoácidos (1-indexed)"""
        self.aa_start = (self.start - 1) // 3 + 1
        self.aa_end = (self.end - 1) // 3 + 1


class GenomicMaps:
    """Mapas genómicos de todos los virus soportados"""
    
    MAPS = {
        'DENV_1': {
            'cds_start': 95,
            'regions': [
                GenomicRegion(95, 436, "POLY", "ancC"),
                GenomicRegion(95, 394, "POLY", "C"),
                GenomicRegion(437, 934, "POLY", "prM"),
                GenomicRegion(437, 709, "POLY", "pr"),
                GenomicRegion(710, 934, "POLY", "M"),
                GenomicRegion(935, 2419, "POLY", "E"),
                GenomicRegion(2420, 3475, "POLY", "NS1"),
                GenomicRegion(3476, 4129, "POLY", "NS2A"),
                GenomicRegion(4130, 4519, "POLY", "NS2B"),
                GenomicRegion(4520, 6376, "POLY", "NS3"),
                GenomicRegion(6377, 6757, "POLY", "NS4A"),
                GenomicRegion(6758, 6826, "POLY", "2K"),
                GenomicRegion(6827, 7573, "POLY", "NS4B"),
                GenomicRegion(7574, 10270, "POLY", "NS5"),
            ]
        },
        'DENV_2': {
            'cds_start': 97,
            'regions': [
                GenomicRegion(97, 438, "POLY", "ancC"),
                GenomicRegion(97, 396, "POLY", "C"),
                GenomicRegion(439, 936, "POLY", "prM"),
                GenomicRegion(439, 711, "POLY", "pr"),
                GenomicRegion(712, 936, "POLY", "M"),
                GenomicRegion(937, 2421, "POLY", "E"),
                GenomicRegion(2422, 3477, "POLY", "NS1"),
                GenomicRegion(3478, 4131, "POLY", "NS2A"),
                GenomicRegion(4132, 4521, "POLY", "NS2B"),
                GenomicRegion(4522, 6375, "POLY", "NS3"),
                GenomicRegion(6376, 6756, "POLY", "NS4A"),
                GenomicRegion(6757, 6825, "POLY", "2K"),
                GenomicRegion(6826, 7569, "POLY", "NS4B"),
                GenomicRegion(7570, 10269, "POLY", "NS5"),
            ]
        },
        'DENV_3': {
            'cds_start': 95,
            'regions': [
                GenomicRegion(95, 436, "POLY", "ancC"),
                GenomicRegion(95, 394, "POLY", "C"),
                GenomicRegion(437, 934, "POLY", "prM"),
                GenomicRegion(437, 709, "POLY", "pr"),
                GenomicRegion(710, 934, "POLY", "M"),
                GenomicRegion(935, 2413, "POLY", "E"),
                GenomicRegion(2414, 3469, "POLY", "NS1"),
                GenomicRegion(3470, 4123, "POLY", "NS2A"),
                GenomicRegion(4124, 4513, "POLY", "NS2B"),
                GenomicRegion(4514, 6370, "POLY", "NS3"),
                GenomicRegion(6371, 6751, "POLY", "NS4A"),
                GenomicRegion(6752, 6820, "POLY", "2K"),
                GenomicRegion(6821, 7564, "POLY", "NS4B"),
                GenomicRegion(7565, 10264, "POLY", "NS5"),
            ]
        },
        'DENV_4': {
            'cds_start': 102,
            'regions': [
                GenomicRegion(102, 440, "POLY", "ancC"),
                GenomicRegion(102, 398, "POLY", "C"),
                GenomicRegion(441, 938, "POLY", "prM"),
                GenomicRegion(441, 713, "POLY", "pr"),
                GenomicRegion(714, 938, "POLY", "M"),
                GenomicRegion(939, 2423, "POLY", "E"),
                GenomicRegion(2424, 3479, "POLY", "NS1"),
                GenomicRegion(3480, 4133, "POLY", "NS2A"),
                GenomicRegion(4134, 4523, "POLY", "NS2B"),
                GenomicRegion(4524, 6377, "POLY", "NS3"),
                GenomicRegion(6378, 6758, "POLY", "NS4A"),
                GenomicRegion(6759, 6827, "POLY", "2K"),
                GenomicRegion(6828, 7562, "POLY", "NS4B"),
                GenomicRegion(7563, 10262, "POLY", "NS5"),
            ]
        },
        'SARS_COV_2': {
            'cds_start': 266,
            'regions': [
                # ORF1ab polyprotein (nsp1-16)
                GenomicRegion(266, 805, "ORF1ab", "nsp1"),
                GenomicRegion(806, 2719, "ORF1ab", "nsp2"),
                GenomicRegion(2720, 8554, "ORF1ab", "nsp3"),
                GenomicRegion(8555, 10054, "ORF1ab", "nsp4"),
                GenomicRegion(10055, 10972, "ORF1ab", "nsp5"),
                GenomicRegion(10973, 11842, "ORF1ab", "nsp6"),
                GenomicRegion(11843, 12091, "ORF1ab", "nsp7"),
                GenomicRegion(12092, 12685, "ORF1ab", "nsp8"),
                GenomicRegion(12686, 13024, "ORF1ab", "nsp9"),
                GenomicRegion(13025, 13441, "ORF1ab", "nsp10"),
                GenomicRegion(13442, 16236, "ORF1ab", "nsp12"),
                GenomicRegion(16237, 18039, "ORF1ab", "nsp13"),
                GenomicRegion(18040, 19620, "ORF1ab", "nsp14"),
                GenomicRegion(19621, 20658, "ORF1ab", "nsp15"),
                GenomicRegion(20659, 21552, "ORF1ab", "nsp16"),
                # Structural and accessory proteins
                GenomicRegion(21563, 25384, "S", "S"),           # Spike protein
                GenomicRegion(25393, 26220, "ORF3a", "ORF3a"),   # ORF3a protein
                GenomicRegion(26245, 26472, "E", "E"),           # Envelope protein
                GenomicRegion(26523, 27191, "M", "M"),           # Membrane protein
                GenomicRegion(27202, 27387, "ORF6", "ORF6"),     # ORF6 protein
                GenomicRegion(27394, 27759, "ORF7a", "ORF7a"),   # ORF7a protein
                GenomicRegion(27756, 27887, "ORF7b", "ORF7b"),   # ORF7b protein
                GenomicRegion(27894, 28259, "ORF8", "ORF8"),     # ORF8 protein
                GenomicRegion(28274, 29533, "N", "N"),           # Nucleocapsid protein
                GenomicRegion(29558, 29674, "ORF10", "ORF10"),   # ORF10 protein
            ]
        },
        'RABV': {
            'cds_start': 71,
            'regions': [
                GenomicRegion(71, 1423, "N", "N"),              # Nucleoprotein (N gene)
                GenomicRegion(1514, 2407, "P", "P"),            # Phosphoprotein (P gene)
                GenomicRegion(2496, 3104, "M", "M"),            # Matrix protein (M gene, M2 region)
                GenomicRegion(3318, 4892, "G", "G"),            # Glycoprotein (G gene)
                GenomicRegion(5418, 11846, "L", "L"),           # Large protein / Polymerase (L gene)
            ]
        },
        'N_RABV': {
            'cds_start': 1,
            'regions': [
                GenomicRegion(1, 1353, "N", "nucleoprotein"),
            ]
        }
    }
    
    @classmethod
    def get_map(cls, virus_type: str) -> Optional[Dict]:
        """Obtiene el mapa genómico para un virus específico"""
        return cls.MAPS.get(virus_type)
    
    @classmethod
    def get_supported_viruses(cls) -> List[str]:
        """Retorna lista de virus soportados"""
        return list(cls.MAPS.keys())


# ==================== TABLA DE CODONES ====================
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


# ==================== FUNCIONES PRINCIPALES ====================
def translate_dna_to_protein(dna_sequence: str, frame: int = 0) -> str:
    """
    Traduce secuencia de ADN a proteína
    
    Args:
        dna_sequence: Secuencia de ADN
        frame: Marco de lectura (0, 1, 2)
    
    Returns:
        Secuencia de proteína
    """
    protein = []
    dna_sequence = dna_sequence.upper()
    
    for i in range(frame, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, 'X')
            protein.append(aa)
    
    return ''.join(protein)


def parse_vcf(vcf_file: str, logger: logging.Logger) -> List[Tuple[int, str, str]]:
    """
    Parse archivo VCF (gzipped o no)
    
    Args:
        vcf_file: Ruta al archivo VCF
        logger: Logger para mensajes
    
    Returns:
        Lista de tuplas (posición, ref, alt)
    """
    variants = []
    
    try:
        # Detectar si está comprimido
        if vcf_file.endswith('.gz'):
            file_handle = gzip.open(vcf_file, 'rt')
        else:
            file_handle = open(vcf_file, 'r')
        
        with file_handle as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    continue
                
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                
                # Solo procesar SNPs (no indels por ahora)
                if len(ref) == 1 and len(alt) == 1:
                    variants.append((pos, ref, alt))
        
        logger.info(f"  ✓ Parsed {len(variants)} SNPs from VCF")
        return variants
    
    except Exception as e:
        logger.error(f"  ✗ Error parsing VCF {vcf_file}: {e}")
        return []


def find_region_for_position(aa_position: int, regions: List[GenomicRegion]) -> Optional[GenomicRegion]:
    """
    Encuentra la región genómica para una posición de aminoácido
    
    Args:
        aa_position: Posición del aminoácido (1-indexed)
        regions: Lista de regiones genómicas
    
    Returns:
        GenomicRegion o None si no se encuentra
    """
    for region in regions:
        if region.aa_start <= aa_position <= region.aa_end:
            return region
    return None


def detect_mutations(
    vcf_file: str,
    reference_fasta: str,
    virus_type: str,
    logger: logging.Logger
) -> Set[str]:
    """
    Detecta mutaciones de aminoácidos desde un VCF
    
    Args:
        vcf_file: Ruta al archivo VCF
        reference_fasta: Ruta al FASTA de referencia
        virus_type: Tipo de virus (DENV_1, DENV_2, etc.)
        logger: Logger
    
    Returns:
        Set de mutaciones en formato Nextclade (ej: "E:D614G")
    """
    mutations = set()
    
    # Obtener mapa genómico
    genome_map = GenomicMaps.get_map(virus_type)
    if not genome_map:
        logger.error(f"  ✗ No genomic map for {virus_type}")
        return mutations
    
    cds_start = genome_map['cds_start']
    regions = genome_map['regions']
    
    # Cargar secuencia de referencia
    try:
        reference_record = next(SeqIO.parse(reference_fasta, "fasta"))
        reference_sequence = str(reference_record.seq).upper()
        logger.info(f"  ✓ Loaded reference: {len(reference_sequence)} bp")
    except Exception as e:
        logger.error(f"  ✗ Error loading reference FASTA: {e}")
        return mutations
    
    # Parsear VCF
    variants = parse_vcf(vcf_file, logger)
    
    if not variants:
        logger.warning("  ⚠ No variants found in VCF")
        return mutations
    
    # Determinar si es poliproteína continua (DENGUE) o genes separados (SARS-CoV-2, RABV)
    is_polyprotein = _is_continuous_polyprotein(regions)
    
    if is_polyprotein:
        logger.info(f"  ℹ Detected continuous polyprotein (e.g., DENGUE)")
        mutations = _process_continuous_polyprotein(
            variants, reference_sequence, regions, cds_start, logger
        )
    else:
        logger.info(f"  ℹ Detected separated genes (e.g., SARS-CoV-2, RABV)")
        mutations = _process_separated_genes(
            variants, reference_sequence, regions, logger
        )
    
    logger.info(f"  ✓ Detected {len(mutations)} amino acid changes")
    return mutations


def _is_continuous_polyprotein(regions: List[GenomicRegion]) -> bool:
    """
    Determina si las regiones forman una poliproteína continua
    
    Args:
        regions: Lista de regiones genómicas
    
    Returns:
        True si es poliproteína continua, False si son genes separados
    """
    # Ordenar regiones por inicio
    sorted_regions = sorted(regions, key=lambda r: r.start)
    
    # Si hay gaps grandes entre regiones (>100 bp), probablemente son genes separados
    for i in range(len(sorted_regions) - 1):
        gap = sorted_regions[i + 1].start - sorted_regions[i].end
        if gap > 100:  # Gap significativo
            return False
    
    return True


def _process_continuous_polyprotein(
    variants: List[Tuple[int, str, str]],
    reference_sequence: str,
    regions: List[GenomicRegion],
    cds_start: int,
    logger: logging.Logger
) -> Set[str]:
    """
    Procesa mutaciones en poliproteína continua (DENGUE)
    """
    mutations = set()
    
    # Encontrar el final de la región codificante
    max_coding_end = max(region.end for region in regions)
    
    # Extraer solo la secuencia codificante (no incluir 3'UTR)
    cds_sequence = reference_sequence[cds_start - 1:max_coding_end]
    reference_protein = translate_dna_to_protein(cds_sequence, frame=0)
    logger.info(f"  ✓ Coding sequence: {len(cds_sequence)} bp ({cds_start}-{max_coding_end})")
    logger.info(f"  ✓ Reference protein: {len(reference_protein)} aa")
    
    # Procesar cada variante
    for nuc_pos, ref_base, alt_base in variants:
        # Verificar que la variante está en la región codificante
        if nuc_pos < cds_start or nuc_pos > max_coding_end:
            logger.debug(f"  ⊘ Skipping pos {nuc_pos}: outside coding region ({cds_start}-{max_coding_end})")
            continue
        
        # CORRECCIÓN CRÍTICA: Calcular posiciones correctamente
        # Posición relativa al inicio del CDS (0-indexed)
        relative_pos = nuc_pos - cds_start
        
        # Posición del aminoácido (0-indexed para acceso a array)
        aa_index = relative_pos // 3
        
        # Posición dentro del codón (0, 1, 2)
        codon_position = relative_pos % 3
        
        # Verificar límites
        if aa_index >= len(reference_protein):
            continue
        
        # Calcular inicio del codón en la secuencia CDS
        codon_start_in_cds = aa_index * 3
        
        # Verificar que tenemos suficiente secuencia
        if codon_start_in_cds + 3 > len(cds_sequence):
            continue
        
        # Extraer codón original
        original_codon = cds_sequence[codon_start_in_cds:codon_start_in_cds + 3]
        
        # Verificar que el nucleótido de referencia coincide
        if original_codon[codon_position] != ref_base:
            logger.debug(
                f"  ⚠ Mismatch at pos {nuc_pos}: "
                f"VCF says {ref_base}, reference has {original_codon[codon_position]}"
            )
        
        # Crear codón mutado
        mutated_codon = (
            original_codon[:codon_position] +
            alt_base +
            original_codon[codon_position + 1:]
        )
        
        # Traducir codones
        original_aa = CODON_TABLE.get(original_codon, 'X')
        mutated_aa = CODON_TABLE.get(mutated_codon, 'X')
        
        # Solo reportar si hay cambio de aminoácido
        if original_aa != mutated_aa:
            # Posición del aminoácido en coordenadas genómicas (1-indexed)
            aa_position_genomic = aa_index + 1
            
            # Encontrar región
            region = find_region_for_position(aa_position_genomic, regions)
            
            if region:
                # Calcular posición relativa dentro de la región
                aa_position_in_region = aa_position_genomic - region.aa_start + 1
                
                # Formato Nextclade: Producto:AminoOriginalPosiciónAminoMutado
                mutation_str = f"{region.product}:{original_aa}{aa_position_in_region}{mutated_aa}"
                mutations.add(mutation_str)
                
                logger.debug(
                    f"  → Mutation: {mutation_str} "
                    f"(nuc pos {nuc_pos}, codon {original_codon}→{mutated_codon})"
                )
    
    return mutations


def _process_separated_genes(
    variants: List[Tuple[int, str, str]],
    reference_sequence: str,
    regions: List[GenomicRegion],
    logger: logging.Logger
) -> Set[str]:
    """
    Procesa mutaciones en genes separados (SARS-CoV-2, RABV)
    Cada región se procesa independientemente
    """
    mutations = set()
    
    # Procesar cada región independientemente
    for region in regions:
        # Extraer secuencia de esta región específica
        region_seq = reference_sequence[region.start - 1:region.end]
        region_protein = translate_dna_to_protein(region_seq, frame=0)
        
        logger.debug(
            f"  Processing {region.product}: {region.start}-{region.end} "
            f"({len(region_seq)} bp, {len(region_protein)} aa)"
        )
        
        # Filtrar variantes que están en esta región
        region_variants = [
            (pos, ref, alt) for pos, ref, alt in variants
            if region.start <= pos <= region.end
        ]
        
        if not region_variants:
            continue
        
        # Procesar cada variante en esta región
        for nuc_pos, ref_base, alt_base in region_variants:
            # Posición relativa al inicio de ESTA región (0-indexed)
            relative_pos = nuc_pos - region.start
            
            # Posición del aminoácido dentro de esta región (0-indexed)
            aa_index = relative_pos // 3
            
            # Posición dentro del codón (0, 1, 2)
            codon_position = relative_pos % 3
            
            # Verificar límites
            if aa_index >= len(region_protein):
                continue
            
            # Calcular inicio del codón en la secuencia de la región
            codon_start_in_region = aa_index * 3
            
            # Verificar que tenemos suficiente secuencia
            if codon_start_in_region + 3 > len(region_seq):
                continue
            
            # Extraer codón original
            original_codon = region_seq[codon_start_in_region:codon_start_in_region + 3]
            
            # Verificar que el nucleótido de referencia coincide
            if original_codon[codon_position] != ref_base:
                logger.debug(
                    f"  ⚠ Mismatch at pos {nuc_pos} in {region.product}: "
                    f"VCF says {ref_base}, reference has {original_codon[codon_position]}"
                )
            
            # Crear codón mutado
            mutated_codon = (
                original_codon[:codon_position] +
                alt_base +
                original_codon[codon_position + 1:]
            )
            
            # Traducir codones
            original_aa = CODON_TABLE.get(original_codon, 'X')
            mutated_aa = CODON_TABLE.get(mutated_codon, 'X')
            
            # Solo reportar si hay cambio de aminoácido
            if original_aa != mutated_aa:
                # Posición del aminoácido dentro de la región (1-indexed)
                aa_position_in_region = aa_index + 1
                
                # Formato Nextclade: Producto:AminoOriginalPosiciónAminoMutado
                mutation_str = f"{region.product}:{original_aa}{aa_position_in_region}{mutated_aa}"
                mutations.add(mutation_str)
                
                logger.debug(
                    f"  → Mutation: {mutation_str} "
                    f"(nuc pos {nuc_pos}, codon {original_codon}→{mutated_codon})"
                )
    
    return mutations


def process_single_vcf(
    vcf_file: str,
    reference_fasta: str,
    virus_type: str,
    output_file: str,
    logger: logging.Logger
) -> bool:
    """
    Procesa un único VCF y guarda resultados
    
    Args:
        vcf_file: Ruta al VCF
        reference_fasta: Ruta al FASTA de referencia
        virus_type: Tipo de virus
        output_file: Archivo de salida
        logger: Logger
    
    Returns:
        True si fue exitoso
    """
    sample_name = Path(vcf_file).stem.replace('_normalized.vcf', '').replace('.vcf', '')
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Processing: {sample_name}")
    logger.info(f"Virus: {virus_type}")
    logger.info(f"{'='*60}")
    
    mutations = detect_mutations(vcf_file, reference_fasta, virus_type, logger)
    
    # Preparar output
    mutations_str = ";".join(sorted(mutations)) if mutations else "No mutations detected"
    
    # Guardar resultados
    try:
        with open(output_file, 'w') as f:
            f.write("sample,virus_type,mutations\n")
            f.write(f"{sample_name},{virus_type},{mutations_str}\n")
        
        logger.info(f"  ✓ Results saved to: {output_file}")
        return True
    
    except Exception as e:
        logger.error(f"  ✗ Error saving results: {e}")
        return False


def process_directory(
    directory: str,
    virus_type: str,
    reference_fasta: str,
    output_dir: str,
    logger: logging.Logger,
    pattern: str = "*_normalized.vcf.gz"
) -> int:
    """
    Procesa todos los VCFs en un directorio
    
    Args:
        directory: Directorio a procesar
        virus_type: Tipo de virus
        reference_fasta: FASTA de referencia
        output_dir: Directorio de salida
        logger: Logger
        pattern: Patrón de archivos VCF
    
    Returns:
        Número de archivos procesados
    """
    from glob import glob
    
    vcf_files = glob(os.path.join(directory, "**", pattern), recursive=True)
    
    if not vcf_files:
        logger.warning(f"No VCF files found matching pattern: {pattern}")
        return 0
    
    logger.info(f"\nFound {len(vcf_files)} VCF files to process")
    
    results = []
    
    for vcf_file in vcf_files:
        sample_name = Path(vcf_file).stem.replace('_normalized.vcf', '').replace('.vcf', '')
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing: {sample_name}")
        logger.info(f"{'='*60}")
        
        mutations = detect_mutations(vcf_file, reference_fasta, virus_type, logger)
        mutations_str = ";".join(sorted(mutations)) if mutations else "No mutations detected"
        
        results.append({
            'sample': sample_name,
            'virus_type': virus_type,
            'mutations': mutations_str,
            'mutation_count': len(mutations)
        })
    
    # Guardar resultados consolidados
    output_file = os.path.join(output_dir, f"mutations_{virus_type}.csv")
    df = pd.DataFrame(results)
    df.to_csv(output_file, index=False)
    
    logger.info(f"\n{'='*60}")
    logger.info(f"✓ Batch processing complete")
    logger.info(f"  Processed: {len(results)} samples")
    logger.info(f"  Results: {output_file}")
    logger.info(f"{'='*60}")
    
    return len(results)


# ==================== MAIN ====================
def main():
    parser = argparse.ArgumentParser(
        description='Universal Mutations Detector - CONSENSO System',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single VCF
  python Universal_Mutations_Detector.py \\
    --vcf sample1.vcf.gz \\
    --reference ref_denv3.fasta \\
    --virus DENV_3 \\
    --output mutations.csv
  
  # Batch processing
  python Universal_Mutations_Detector.py \\
    --directory /path/to/samples \\
    --reference ref_denv3.fasta \\
    --virus DENV_3 \\
    --output-dir Mutations_Results
  
  # Auto-detect mode (for automatic assembler)
  python Universal_Mutations_Detector.py \\
    --vcf sample1_DENV_2.vcf.gz \\
    --reference ref_denv2.fasta \\
    --auto-detect
        """
    )
    
    # Modo de entrada
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--vcf', help='Single VCF file to process')
    input_group.add_argument('--directory', help='Directory with VCF files (batch mode)')
    
    # Referencia
    parser.add_argument('--reference', required=True, help='Reference FASTA file')
    
    # Tipo de virus
    virus_group = parser.add_mutually_exclusive_group(required=True)
    virus_group.add_argument(
        '--virus',
        choices=GenomicMaps.get_supported_viruses(),
        help='Virus type'
    )
    virus_group.add_argument(
        '--auto-detect',
        action='store_true',
        help='Auto-detect virus from VCF filename (for automatic assembler)'
    )
    
    # Output
    parser.add_argument('--output', help='Output CSV file (single mode)')
    parser.add_argument('--output-dir', help='Output directory (batch mode)', default='Mutations_Results')
    
    # Opciones adicionales
    parser.add_argument('--pattern', default='*_normalized.vcf.gz', help='VCF file pattern (batch mode)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Crear directorio de salida
    output_dir = args.output_dir if args.directory else os.path.dirname(args.output) or '.'
    os.makedirs(output_dir, exist_ok=True)
    
    # Setup logging
    logger = setup_logging(output_dir, args.verbose)
    
    logger.info("="*60)
    logger.info("Universal Mutations Detector - CONSENSO System")
    logger.info("="*60)
    
    # Determinar tipo de virus
    if args.auto_detect:
        if args.vcf:
            # Extraer de nombre de archivo
            filename = os.path.basename(args.vcf)
            virus_type = None
            for v in GenomicMaps.get_supported_viruses():
                if v in filename:
                    virus_type = v
                    break
            
            if not virus_type:
                logger.error("Could not auto-detect virus type from filename")
                return 1
            
            logger.info(f"Auto-detected virus: {virus_type}")
        else:
            logger.error("Auto-detect mode requires --vcf")
            return 1
    else:
        virus_type = args.virus
    
    # Verificar referencia
    if not os.path.exists(args.reference):
        logger.error(f"Reference file not found: {args.reference}")
        return 1
    
    # Procesar
    if args.vcf:
        # Modo single
        if not os.path.exists(args.vcf):
            logger.error(f"VCF file not found: {args.vcf}")
            return 1
        
        output_file = args.output or os.path.join(output_dir, f"{Path(args.vcf).stem}_mutations.csv")
        
        success = process_single_vcf(
            args.vcf,
            args.reference,
            virus_type,
            output_file,
            logger
        )
        
        return 0 if success else 1
    
    else:
        # Modo batch
        if not os.path.isdir(args.directory):
            logger.error(f"Directory not found: {args.directory}")
            return 1
        
        processed = process_directory(
            args.directory,
            virus_type,
            args.reference,
            output_dir,
            logger,
            args.pattern
        )
        
        return 0 if processed > 0 else 1


if __name__ == '__main__':
    sys.exit(main())
