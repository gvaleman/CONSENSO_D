#!/usr/bin/env python3
import csv
import os
import re
import shutil

class SampleSheet:
    def __init__(self, csv_path=None):
        self.samples = []
        self.header = {}
        if csv_path:
            self.load(csv_path)
    
    def load(self, csv_path):
        """Carga un archivo CSV de sample sheet de Illumina"""
        self.samples = []
        self.header = {}
        
        try:
            with open(csv_path, 'r') as f:
                # Leer el archivo CSV
                lines = f.readlines()
                
                # Procesar secciones del sample sheet
                section = None
                header_data = {}
                sample_data = []
                
                for line in lines:
                    line = line.strip()
                    if not line:
                        continue
                    
                    # Detectar secciones
                    if line.startswith('['):
                        section = line.strip('[]')
                        continue
                    
                    # Procesar sección Header
                    if section == 'Header':
                        if ',' in line:
                            key, value = line.split(',', 1)
                            header_data[key.strip()] = value.strip()
                    
                    # Procesar sección Data
                    elif section == 'Data':
                        if line.startswith('Sample_ID') or line.startswith('SampleID'):
                            # Cabecera de columnas
                            headers = [h.strip() for h in line.split(',')]
                        else:
                            # Datos de muestra
                            values = [v.strip() for v in line.split(',')]
                            if len(values) >= len(headers):
                                sample = {headers[i]: values[i] for i in range(len(headers))}
                                sample_data.append(sample)
            
            self.header = header_data
            self.samples = sample_data
            return True
            
        except Exception as e:
            print(f"Error al cargar el sample sheet: {e}")
            return False
    
    def get_sample_names(self):
        """Retorna un diccionario con Sample_ID como clave y Sample_Name como valor"""
        sample_dict = {}
        for sample in self.samples:
            # Buscar las claves correctas (pueden variar según el formato)
            sample_id = sample.get('Sample_ID', sample.get('SampleID', ''))
            sample_name = sample.get('Sample_Name', sample.get('SampleName', sample_id))
            
            if sample_id:
                sample_dict[sample_id] = sample_name
        
        return sample_dict
    
    def rename_consensus_files(self, consensus_file, output_dir=None):
        """
        Renombra o divide el archivo de consenso según los nombres del sample sheet
        
        Args:
            consensus_file: Ruta al archivo all_consensus.fasta
            output_dir: Directorio de salida (opcional)
        
        Returns:
            Lista de archivos generados
        """
        if not os.path.exists(consensus_file):
            print(f"Error: No se encontró el archivo {consensus_file}")
            return []
        
        # Crear directorio de salida si no existe
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Obtener diccionario de nombres de muestra
        sample_names = self.get_sample_names()
        if not sample_names:
            print("Error: No se encontraron muestras en el sample sheet")
            return []
        
        # Leer el archivo de consenso
        with open(consensus_file, 'r') as f:
            content = f.read()
        
        # Dividir por secuencias (cabeceras FASTA)
        sequences = re.split(r'(>.*?\n)', content)
        if len(sequences) <= 1:
            print(f"Error: No se encontraron secuencias en {consensus_file}")
            return []
        
        # Agrupar cabeceras con secuencias
        fasta_entries = []
        current_header = None
        
        for item in sequences:
            if item.startswith('>'):
                current_header = item
            elif current_header:
                fasta_entries.append((current_header, item))
                current_header = None
        
        # Crear archivos individuales o renombrar
        output_files = []
        
        for header, sequence in fasta_entries:
            # Extraer ID de la cabecera
            header_id = header.strip().split()[0][1:]  # Quitar '>' y tomar hasta el primer espacio
            
            # Buscar nombre correspondiente en el sample sheet
            sample_name = None
            for sample_id, name in sample_names.items():
                # Verificar si el ID de la cabecera contiene el ID de la muestra
                if sample_id in header_id:
                    sample_name = name
                    break
            
            # Si no se encontró, usar el ID original
            if not sample_name:
                sample_name = header_id
            
            # Crear nuevo archivo o añadir al archivo de consenso renombrado
            if output_dir:
                output_file = os.path.join(output_dir, f"{sample_name}.fasta")
            else:
                output_file = f"{sample_name}.fasta"
            
            with open(output_file, 'w') as f:
                f.write(f">{sample_name}\n{sequence}")
            
            output_files.append(output_file)
            print(f"Creado archivo: {output_file}")
        
        return output_files

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Procesar sample sheet de Illumina y renombrar archivos de consenso")
    parser.add_argument("--sample_sheet", required=True, help="Archivo CSV de sample sheet")
    parser.add_argument("--consensus", required=True, help="Archivo FASTA de consenso")
    parser.add_argument("--output_dir", help="Directorio de salida para los archivos renombrados")
    
    args = parser.parse_args()
    
    # Procesar sample sheet
    sheet = SampleSheet(args.sample_sheet)
    
    # Renombrar archivos de consenso
    sheet.rename_consensus_files(args.consensus, args.output_dir)

if __name__ == "__main__":
    main()
