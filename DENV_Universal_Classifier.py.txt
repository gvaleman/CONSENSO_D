#!/usr/bin/env python3
"""
DENV Viral Branch - Clasificador Universal de Linajes (Versión Optimizada)
Procesa VCFs automáticamente desde directorios de ensamblaje y predice linajes DENV1-4
Autor: Instituto De Ciencias Sostenibles
Versión: 2.1 - Optimizada para búsqueda automática de VCFs
"""

import pandas as pd
import numpy as np
import joblib
from pathlib import Path
from cyvcf2 import VCF
import argparse
import sys
from datetime import datetime
import warnings
import glob
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

warnings.filterwarnings('ignore')

class DENVUniversalClassifier:
    """Clasificador universal de linajes DENV optimizado para procesamiento batch."""
    
    def __init__(self, model_path, encoder_path, serotype=None):
        """
        Inicializa el clasificador.
        
        Parameters:
        -----------
        model_path : str
            Ruta al modelo entrenado (.joblib)
        encoder_path : str
            Ruta al label encoder (.joblib)
        serotype : str, optional
            Serotipo DENV (DENV1, DENV2, DENV3, DENV4)
        """
        self.model_path = Path(model_path)
        self.encoder_path = Path(encoder_path)
        self.serotype = serotype
        self.model = None
        self.label_encoder = None
        self.feature_names = None
        
        # Umbrales de confianza
        self.HIGH_CONFIDENCE_THRESHOLD = 0.8
        self.MEDIUM_CONFIDENCE_THRESHOLD = 0.6
        self.AMBIGUITY_THRESHOLD = 0.2
        
        # Auto-detectar serotipo desde el path si no se especifica
        if self.serotype is None:
            self.serotype = self._detect_serotype_from_path()
    
    def _detect_serotype_from_path(self):
        """Detecta el serotipo automáticamente desde la ruta del modelo."""
        path_str = str(self.model_path)
        
        serotype_mapping = {
            'VB_DENV1': 'DENV1',
            'VB_DENV2': 'DENV2', 
            'VB_DENV3': 'DENV3',
            'VB_DENV4': 'DENV4',
            'DENV_1': 'DENV1',
            'DENV_2': 'DENV2',
            'DENV_3': 'DENV3',
            'DENV_4': 'DENV4'
        }
        
        for key, value in serotype_mapping.items():
            if key in path_str:
                return value
        
        return "DENV_UNKNOWN"
    
    @staticmethod
    def find_vcfs_recursive(directory, pattern="*_normalized.vcf.gz"):
        """
        Busca VCFs recursivamente en un directorio y sus subdirectorios.
        
        Parameters:
        -----------
        directory : str/Path
            Directorio raíz para buscar
        pattern : str
            Patrón de archivos a buscar
            
        Returns:
        --------
        list : Lista de tuplas (vcf_path, sample_name)
        """
        directory = Path(directory)
        vcf_files = []
        
        # Buscar recursivamente
        for vcf_path in directory.rglob(pattern):
            if vcf_path.is_file():
                # Determinar nombre de muestra basado en la estructura de directorios
                sample_name = vcf_path.parent.name
                
                # Si está en el directorio raíz, usar nombre del archivo
                if sample_name == directory.name:
                    sample_name = vcf_path.stem.replace('_normalized.vcf', '').replace('.gz', '')
                
                vcf_files.append((str(vcf_path), sample_name))
        
        return vcf_files
        
    def load_model(self):
        """Carga el modelo y el encoder."""
        
        if not self.model_path.exists():
            raise FileNotFoundError(f"Modelo no encontrado: {self.model_path}")
        
        if not self.encoder_path.exists():
            raise FileNotFoundError(f"Encoder no encontrado: {self.encoder_path}")
        
        try:
            self.model = joblib.load(self.model_path)
            self.label_encoder = joblib.load(self.encoder_path)
            
            print(f"Modelo cargado exitosamente")
            print(f"   Serotipo: {self.serotype}")
            print(f"   Modelo: {self.model_path.name}")
            print(f"   Encoder: {self.encoder_path.name}")
            print(f"   Linajes disponibles: {len(self.label_encoder.classes_)}")
            
            # Obtener nombres de características del modelo entrenado
            if hasattr(self.model, 'feature_names_in_'):
                self.feature_names = self.model.feature_names_in_
                print(f"   Características: {len(self.feature_names)}")
            else:
                print("   Advertencia: No se pudieron obtener los nombres de características")
            
        except Exception as e:
            raise Exception(f"Error cargando modelo: {e}")
    
    def process_vcf(self, vcf_path):
        """
        Procesa un archivo VCF y extrae variantes.
        
        Parameters:
        -----------
        vcf_path : str/Path
            Ruta al archivo VCF
            
        Returns:
        --------
        dict : Diccionario con variantes encontradas
        """
        vcf_path = Path(vcf_path)
        
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF no encontrado: {vcf_path}")
        
        variants = {}
        
        try:
            vcf = VCF(str(vcf_path))
            
            for variant in vcf:
                variant_key = f"{variant.CHROM}_{variant.POS}_{variant.REF}_{','.join(map(str, variant.ALT))}"
                variants[variant_key] = 1
                
        except Exception as e:
            raise Exception(f"Error procesando VCF {vcf_path}: {e}")
        
        return variants
    
    def prepare_feature_matrix(self, variants_dict, sample_name):
        """
        Prepara la matriz de características para predicción.
        
        Parameters:
        -----------
        variants_dict : dict
            Diccionario con variantes de la muestra
        sample_name : str
            Nombre de la muestra
            
        Returns:
        --------
        tuple : (pd.DataFrame, float) - Matriz de características y cobertura
        """
        
        if self.feature_names is None:
            raise Exception("Modelo no tiene nombres de características definidos")
        
        # Crear DataFrame con todas las características del modelo entrenado
        feature_matrix = pd.DataFrame(0, 
                                    index=[sample_name], 
                                    columns=self.feature_names,
                                    dtype=int)
        
        # Marcar variantes presentes
        variants_found = 0
        for variant in variants_dict:
            if variant in feature_matrix.columns:
                feature_matrix.loc[sample_name, variant] = 1
                variants_found += 1
        
        coverage = (variants_found / len(self.feature_names)) * 100
        
        return feature_matrix, coverage
    
    def predict_with_confidence(self, feature_matrix):
        """
        Realiza predicción con análisis de confianza detallado.
        
        Parameters:
        -----------
        feature_matrix : pd.DataFrame
            Matriz de características
            
        Returns:
        --------
        dict : Resultados de predicción con métricas de confianza
        """
        
        # Obtener probabilidades para todas las clases
        probabilities = self.model.predict_proba(feature_matrix)[0]
        
        # Predicción estándar (clase con mayor probabilidad)
        predicted_class_idx = np.argmax(probabilities)
        predicted_lineage = self.label_encoder.classes_[predicted_class_idx]
        
        # Métricas de confianza
        max_probability = np.max(probabilities)
        
        # Obtener top 3 predicciones
        top3_indices = np.argsort(probabilities)[-3:][::-1]
        top3_predictions = []
        
        for i, idx in enumerate(top3_indices):
            lineage = self.label_encoder.classes_[idx]
            prob = probabilities[idx]
            top3_predictions.append({
                'rank': i + 1,
                'lineage': lineage,
                'probability': prob,
                'percentage': prob * 100
            })
        
        # Calcular margen de confianza (diferencia entre 1° y 2°)
        if len(top3_predictions) >= 2:
            confidence_margin = top3_predictions[0]['probability'] - top3_predictions[1]['probability']
        else:
            confidence_margin = max_probability
        
        # Clasificar nivel de confianza
        if max_probability >= self.HIGH_CONFIDENCE_THRESHOLD:
            confidence_level = "ALTA"
        elif max_probability >= self.MEDIUM_CONFIDENCE_THRESHOLD:
            confidence_level = "MEDIA"
        else:
            confidence_level = "BAJA"
        
        # Detectar predicciones ambiguas
        is_ambiguous = confidence_margin < self.AMBIGUITY_THRESHOLD
        
        # Calcular entropía normalizada
        entropy = -np.sum(probabilities * np.log2(probabilities + 1e-10))
        max_entropy = np.log2(len(probabilities))
        normalized_entropy = entropy / max_entropy
        
        return {
            'serotype': self.serotype,
            'predicted_lineage': predicted_lineage,
            'max_probability': max_probability,
            'confidence_level': confidence_level,
            'confidence_margin': confidence_margin,
            'is_ambiguous': is_ambiguous,
            'normalized_entropy': normalized_entropy,
            'top3_predictions': top3_predictions
        }
    
    def classify_single_sample(self, vcf_path, sample_name=None):
        """
        Clasifica una muestra individual desde VCF.
        
        Parameters:
        -----------
        vcf_path : str/Path
            Ruta al archivo VCF
        sample_name : str, optional
            Nombre de la muestra
            
        Returns:
        --------
        dict : Resultados de clasificación
        """
        
        vcf_path = Path(vcf_path)
        
        if sample_name is None:
            sample_name = vcf_path.stem.replace('_normalized.vcf', '').replace('.gz', '')
        
        # Procesar VCF
        variants = self.process_vcf(vcf_path)
        
        # Preparar matriz de características
        feature_matrix, coverage = self.prepare_feature_matrix(variants, sample_name)
        
        # Realizar predicción
        results = self.predict_with_confidence(feature_matrix)
        
        # Agregar información adicional
        results.update({
            'sample_name': sample_name,
            'vcf_path': str(vcf_path),
            'variants_found': len(variants),
            'feature_coverage': coverage,
            'timestamp': datetime.now().isoformat()
        })
        
        return results
    
    def classify_batch(self, vcf_files, max_workers=None):
        """
        Clasifica múltiples VCFs en paralelo.
        
        Parameters:
        -----------
        vcf_files : list
            Lista de tuplas (vcf_path, sample_name)
        max_workers : int, optional
            Número de procesos paralelos
            
        Returns:
        --------
        list : Lista de resultados de clasificación
        """
        
        if max_workers is None:
            max_workers = min(len(vcf_files), multiprocessing.cpu_count())
        
        print(f"Procesando {len(vcf_files)} VCFs con {max_workers} procesos paralelos...")
        
        all_results = []
        
        # Procesamiento secuencial para evitar problemas con joblib y multiprocessing
        for i, (vcf_path, sample_name) in enumerate(vcf_files, 1):
            try:
                print(f"   Procesando {i}/{len(vcf_files)}: {sample_name}")
                result = self.classify_single_sample(vcf_path, sample_name)
                all_results.append(result)
                
            except Exception as e:
                print(f"   Error procesando {sample_name}: {e}")
                continue
        
        return all_results
    
    def save_batch_results(self, results, output_file):
        """
        Guarda resultados batch en archivo CSV.
        
        Parameters:
        -----------
        results : list
            Lista de resultados de clasificación
        output_file : str/Path
            Archivo de salida
        """
        
        if not results:
            print("No hay resultados para guardar")
            return
        
        # Preparar datos para CSV
        csv_data = []
        for result in results:
            row = {
                'Sample_Name': result['sample_name'],
                'Serotype': result['serotype'],
                'Predicted_Lineage': result['predicted_lineage'],
                'Confidence_Level': result['confidence_level'],
                'Max_Probability': result['max_probability'],
                'Confidence_Percentage': result['max_probability'] * 100,
                'Is_Ambiguous': result['is_ambiguous'],
                'Variants_Found': result['variants_found'],
                'Feature_Coverage': result['feature_coverage'],
                'Timestamp': result['timestamp']
            }
            
            # Agregar top 3 predicciones
            for i, pred in enumerate(result['top3_predictions']):
                row[f'Top{i+1}_Lineage'] = pred['lineage']
                row[f'Top{i+1}_Probability'] = pred['probability']
                row[f'Top{i+1}_Percentage'] = pred['percentage']
            
            csv_data.append(row)
        
        # Crear DataFrame y guardar
        df = pd.DataFrame(csv_data)
        df.to_csv(output_file, index=False)
        
        print(f"Resultados guardados en: {output_file}")
        
        # Mostrar estadísticas rápidas
        print(f"   Total de muestras: {len(results)}")
        
        # Distribución por confianza
        confidence_counts = df['Confidence_Level'].value_counts()
        print("   Distribución por confianza:")
        for level, count in confidence_counts.items():
            percentage = (count / len(results)) * 100
            print(f"      {level}: {count} ({percentage:.1f}%)")
        
        # Top linajes
        lineage_counts = df['Predicted_Lineage'].value_counts().head(5)
        print("   Top 5 linajes identificados:")
        for lineage, count in lineage_counts.items():
            percentage = (count / len(results)) * 100
            print(f"      {lineage}: {count} ({percentage:.1f}%)")

def process_directory_batch(model_path, encoder_path, search_dir, serotype=None, output_file=None, pattern="*_normalized.vcf.gz"):
    """
    Función principal para procesamiento batch desde un directorio.
    
    Parameters:
    -----------
    model_path : str
        Ruta al modelo
    encoder_path : str
        Ruta al encoder
    search_dir : str
        Directorio donde buscar VCFs
    serotype : str, optional
        Serotipo a procesar
    output_file : str, optional
        Archivo de salida
    pattern : str
        Patrón de búsqueda para VCFs
    """
    
    # Inicializar clasificador
    classifier = DENVUniversalClassifier(model_path, encoder_path, serotype)
    classifier.load_model()
    
    # Buscar VCFs recursivamente
    print(f"Buscando VCFs en: {search_dir}")
    vcf_files = classifier.find_vcfs_recursive(search_dir, pattern)
    
    if not vcf_files:
        print(f"No se encontraron VCFs con patrón '{pattern}' en {search_dir}")
        return None
    
    print(f"Encontrados {len(vcf_files)} VCFs para procesar")
    for vcf_path, sample_name in vcf_files[:5]:  # Mostrar primeros 5
        print(f"   {sample_name}: {Path(vcf_path).name}")
    if len(vcf_files) > 5:
        print(f"   ... y {len(vcf_files) - 5} más")
    
    # Procesar en batch
    results = classifier.classify_batch(vcf_files)
    
    if not results:
        print("No se procesaron muestras exitosamente")
        return None
    
    # Generar nombre de archivo de salida si no se especifica
    if output_file is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        serotype_str = classifier.serotype.lower()
        output_file = f"batch_{serotype_str}_classifications_{timestamp}.csv"
    
    # Guardar resultados
    classifier.save_batch_results(results, output_file)
    
    return output_file

def main():
    """Función principal con interfaz de línea de comandos optimizada."""
    
    parser = argparse.ArgumentParser(
        description="DENV Universal Classifier - Versión optimizada para búsqueda automática de VCFs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:

1. Buscar y clasificar VCFs automáticamente en un directorio:
   python DENV_Universal_Classifier.py --search-dir ./samples/ --model model.joblib --encoder encoder.joblib

2. Especificar patrón de búsqueda personalizado:
   python DENV_Universal_Classifier.py --search-dir ./ --pattern "*calls.vcf.gz" --model model.joblib --encoder encoder.joblib

3. Clasificar un directorio específico (modo original):
   python DENV_Universal_Classifier.py --vcf-dir ./vcfs/ --model model.joblib --encoder encoder.joblib

4. Clasificar un VCF individual:
   python DENV_Universal_Classifier.py --vcf sample.vcf.gz --model model.joblib --encoder encoder.joblib
        """
    )
    
    # Argumentos principales
    parser.add_argument('--vcf', type=str, help='Archivo VCF individual para clasificar')
    parser.add_argument('--vcf-dir', type=str, help='Directorio con VCFs (búsqueda directa)')
    parser.add_argument('--search-dir', type=str, help='Directorio para búsqueda recursiva de VCFs')
    parser.add_argument('--model', type=str, required=True, help='Ruta al modelo entrenado (.joblib)')
    parser.add_argument('--encoder', type=str, required=True, help='Ruta al label encoder (.joblib)')
    
    # Argumentos opcionales
    parser.add_argument('--serotype', type=str, help='Serotipo DENV (DENV1, DENV2, DENV3, DENV4)')
    parser.add_argument('--output', type=str, help='Archivo de salida CSV')
    parser.add_argument('--pattern', type=str, default="*_normalized.vcf.gz", 
                       help='Patrón de búsqueda para VCFs (default: *_normalized.vcf.gz)')
    parser.add_argument('--sample-name', type=str, help='Nombre personalizado para muestra individual')
    parser.add_argument('--quiet', action='store_true', help='Modo silencioso')
    
    args = parser.parse_args()
    
    # Validar argumentos
    modes = [args.vcf, args.vcf_dir, args.search_dir]
    if sum(x is not None for x in modes) != 1:
        parser.error("Especifique exactamente una opción: --vcf, --vcf-dir, o --search-dir")
    
    try:
        if args.search_dir:
            # Modo de búsqueda recursiva (NUEVO - RECOMENDADO)
            output_file = process_directory_batch(
                args.model, args.encoder, args.search_dir, 
                args.serotype, args.output, args.pattern
            )
            if output_file:
                print(f"\nProcesamiento completado. Resultados en: {output_file}")
            
        elif args.vcf_dir:
            # Modo de directorio directo (ORIGINAL)
            classifier = DENVUniversalClassifier(args.model, args.encoder, args.serotype)
            classifier.load_model()
            
            vcf_dir = Path(args.vcf_dir)
            if not vcf_dir.exists():
                print(f"Directorio no encontrado: {vcf_dir}")
                sys.exit(1)
            
            # Buscar VCFs directamente en el directorio
            vcf_files = []
            for pattern in ["*.vcf", "*.vcf.gz"]:
                vcf_files.extend([(str(f), f.stem.replace('.vcf', '').replace('.gz', '')) 
                                for f in vcf_dir.glob(pattern)])
            
            if not vcf_files:
                print(f"No se encontraron VCFs en: {vcf_dir}")
                sys.exit(1)
            
            results = classifier.classify_batch(vcf_files)
            
            if results:
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                output_file = args.output or f"{classifier.serotype.lower()}_batch_{timestamp}.csv"
                classifier.save_batch_results(results, output_file)
            
        elif args.vcf:
            # Modo de VCF individual (ORIGINAL)
            classifier = DENVUniversalClassifier(args.model, args.encoder, args.serotype)
            classifier.load_model()
            
            result = classifier.classify_single_sample(args.vcf, args.sample_name)
            
            if not args.quiet:
                # Mostrar resultados detallados
                print(f"\nRESULTADOS DE CLASIFICACION {result['serotype']}")
                print(f"="*50)
                print(f"Muestra: {result['sample_name']}")
                print(f"Linaje: {result['predicted_lineage']}")
                print(f"Confianza: {result['confidence_level']} ({result['max_probability']:.3f})")
                print(f"Variantes: {result['variants_found']}")
                print(f"Cobertura: {result['feature_coverage']:.1f}%")
            
            # Guardar resultado
            output_file = args.output or f"{result['sample_name']}_classification.csv"
            classifier.save_batch_results([result], output_file)
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
