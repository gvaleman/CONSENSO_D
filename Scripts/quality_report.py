#!/usr/bin/env python3
"""
Script simple para generar reportes de calidad
Uso: python quality_report.py [directorio_muestras]
"""

import os
import sys
import glob
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot
import datetime

def main():
    # Directorio de trabajo
    work_dir = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
    
    print(f"Procesando: {work_dir}")
    
    # COBERTURA - buscar *.coverage
    print("\n=== COBERTURA ===")
    coverage_files = glob.glob(os.path.join(work_dir, "**", "*.coverage"), recursive=True)
    coverage_plots = ""
    coverage_stats = []
    
    for coverage_file in coverage_files:
        # Nombre de la muestra = nombre del archivo de cobertura sin la extensión
        sample_name = os.path.basename(coverage_file).replace(".coverage", "")
        
        try:
            # Leer el archivo de cobertura con los nombres correctos de columnas
            df = pd.read_csv(coverage_file, sep="\t", header=None, names=['Ref', 'Pos', 'Cov'])
            df['Cov'] = pd.to_numeric(df['Cov'], errors='coerce').fillna(0)
            df['Pos'] = pd.to_numeric(df['Pos'], errors='coerce').fillna(0)
            
            # Debug: imprimir info del archivo
            print(f"  {sample_name}: {len(df)} posiciones, rango {df['Pos'].min()}-{df['Pos'].max()}")
            
            # Crear gráfico de barras para la cobertura
            fig = go.Figure()
            fig.add_trace(go.Bar(
                x=df['Pos'], 
                y=df['Cov'], 
                name="Cobertura", 
                marker_color='blue',
                width=1  # Ancho de barras para que se vean todas las posiciones
            ))
            
            # Calcular rangos dinámicos basados en los datos
            pos_min = df['Pos'].min()
            pos_max = df['Pos'].max()
            cov_max = df['Cov'].max()
            
            # Configurar el diseño con rangos dinámicos
            fig.update_layout(
                title=f"{sample_name} (Posiciones: {pos_min:,}-{pos_max:,})", 
                xaxis_title="Posición Genómica", 
                yaxis_title="Profundidad de Cobertura", 
                height=500, 
                template="plotly_white",
                xaxis=dict(
                    range=[pos_min - 50, pos_max + 50],  # Pequeño margen
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='lightgray',
                    tickmode='linear',
                    dtick=max(1, (pos_max - pos_min) // 10)  # ~10 ticks en el eje X
                ),
                yaxis=dict(
                    range=[0, 1500],  # RANGO FIJO como solicitaste
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='lightgray'
                ),
                bargap=0,
                showlegend=False,
                margin=dict(l=60, r=60, t=80, b=80)
            )
            
            # Línea de referencia de 1000x
            fig.add_hline(y=1000, line_dash="dash", line_color="red", 
                         annotation_text="1000x", annotation_position="top right")
            
            plot_html = plot(fig, output_type='div', include_plotlyjs=True)
            coverage_plots += f'<div style="margin: 20px 0; padding: 20px; background: white; border-radius: 10px;"><h3>{sample_name}</h3>{plot_html}</div>'
            
            # Estadísticas para el resumen
            avg_cov = df['Cov'].mean()
            coverage_stats.append(avg_cov)
            
            print(f"✓ {sample_name}: {avg_cov:.0f}x promedio")
            
        except Exception as e:
            print(f"❌ Error en {sample_name}: {e}")
            import traceback
            traceback.print_exc()  # Para debug detallado
    
    # READS - buscar *general_stats.txt
    print("\n=== READS ===")
    stats_files = glob.glob(os.path.join(work_dir, "**", "*general_stats.txt"), recursive=True)
    reads_data = {}
    
    for stats_file in stats_files:
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(stats_file)))
        
        try:
            df = pd.read_csv(stats_file, sep="\t")
            # Buscar columnas con 'sequences' o 'reads'
            seq_cols = [col for col in df.columns if any(x in col.lower() for x in ['sequences', 'reads'])]
            
            if seq_cols:
                total_reads = df[seq_cols].sum().sum()
                # Convertir a millones si es muy pequeño (probablemente ya está en millones)
                if total_reads < 10:
                    total_reads = total_reads * 1
                
                reads_data[sample_name] = int(total_reads)
                print(f"✓ {sample_name}: {total_reads:,} reads")
        except Exception as e:
            print(f"❌ Error en {sample_name}: {e}")
    
    # Gráfico de reads
    reads_plot = ""
    if reads_data:
        fig = go.Figure(data=[go.Bar(x=list(reads_data.keys()), y=list(reads_data.values()), 
                                   marker_color='steelblue')])
        fig.update_layout(title="Reads Totales por Muestra", xaxis_title="Muestra", 
                         yaxis_title="Número de Reads", template="plotly_white",
                         xaxis={'tickangle': 45})
        reads_plot = plot(fig, output_type='div', include_plotlyjs=True)
    
    # Q-SCORES - buscar *fastqc.txt
    print("\n=== Q-SCORES ===")
    fastqc_files = glob.glob(os.path.join(work_dir, "**", "*fastqc.txt"), recursive=True)
    qscore_data = {}
    
    for fastqc_file in fastqc_files:
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(fastqc_file)))
        
        try:
            df = pd.read_csv(fastqc_file, sep="\t")
            
            # Buscar columnas de calidad
            quality_cols = [col for col in df.columns if 'quality' in col.lower() and 'mean' in col.lower()]
            
            if quality_cols and not df.empty:
                # Filtrar solo R1 si existe
                r1_data = df[df['Sample'].str.contains('R1', case=False, na=False)] if 'Sample' in df.columns else df
                target_df = r1_data if not r1_data.empty else df
                
                score = target_df[quality_cols[0]].iloc[0]
                if pd.notna(score):
                    qscore_data[sample_name] = float(score)
                    print(f"✓ {sample_name}: {score:.1f}")
        except Exception as e:
            print(f"❌ Error en {sample_name}: {e}")
    
    # Gráfico de Q-scores
    qscore_plot = ""
    if qscore_data:
        fig = go.Figure(data=[go.Bar(x=list(qscore_data.keys()), y=list(qscore_data.values()), 
                                   marker_color='green')])
        fig.update_layout(title="Q-score Promedio por Muestra", xaxis_title="Muestra", 
                         yaxis_title="Q-score", template="plotly_white",
                         xaxis={'tickangle': 45})
        qscore_plot = plot(fig, output_type='div', include_plotlyjs=True)
    
    # Estadísticas finales
    avg_coverage = sum(coverage_stats) / len(coverage_stats) if coverage_stats else 0
    samples_above_1000x = sum(1 for cov in coverage_stats if cov >= 1000)
    total_samples = len(coverage_files)
    current_year = datetime.datetime.now().year
    
    # HTML simple
    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Reporte de Calidad</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }}
        .header {{ background: #2c3e50; color: white; padding: 25px; border-radius: 15px; text-align: center; margin-bottom: 25px; }}
        .stats {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin: 20px 0; }}
        .stat {{ background: white; padding: 20px; border-radius: 10px; text-align: center; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        .stat-number {{ font-size: 2em; font-weight: bold; color: #3498db; }}
        .charts {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin: 20px 0; }}
        .chart {{ background: white; padding: 20px; border-radius: 10px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        .coverage {{ background: white; padding: 20px; border-radius: 10px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); margin: 20px 0; }}
        @media (max-width: 768px) {{ .charts {{ grid-template-columns: 1fr; }} }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Reporte de Calidad de Secuenciación Y Ensamblaje de Genomas</h1>
        <h2>CONSENSO</h2>
        <p>Análisis Automático - {pd.Timestamp.now().strftime('%d/%m/%Y %H:%M')}</p>
    </div>
    
    <div class="stats">
        <div class="stat">
            <div class="stat-number">{total_samples}</div>
            <div>Muestras</div>
        </div>
        <div class="stat">
            <div class="stat-number">{avg_coverage:.0f}x</div>
            <div>Cobertura Promedio</div>
        </div>
        <div class="stat">
            <div class="stat-number">{samples_above_1000x}</div>
            <div>Muestras ≥1000x</div>
        </div>
        <div class="stat">
            <div class="stat-number">{len(reads_data)}</div>
            <div>Con Datos FastQC</div>
        </div>
    </div>
    
    <div class="charts">
        <div class="chart">
            <h3>Reads Totales por Muestra</h3>
            {reads_plot if reads_plot else '<p style="text-align: center; color: #666;">No hay datos de reads disponibles</p>'}
        </div>
        <div class="chart">
            <h3>Q-scores Promedio</h3>
            {qscore_plot if qscore_plot else '<p style="text-align: center; color: #666;">No hay datos de Q-score disponibles</p>'}
        </div>
    </div>
    
    <div class="coverage">
        <h3>Cobertura Genómica por Muestra</h3>
        <p style="margin-bottom: 20px; color: #666;">La línea roja indica 1000x de cobertura (óptima para análisis viral)</p>
        {coverage_plots if coverage_plots else '<p style="text-align: center; color: #666;">No hay datos de cobertura disponibles</p>'}
    </div>

    <div style="text-align: center; margin-top: 20px; color: #666;">
        <p>Archivos procesados: {len(coverage_files)} coverage | {len(stats_files)} stats | {len(fastqc_files)} fastqc</p>
    </div>

    <div style="text-align: center; margin-top: 40px; color: #666; font-size: 0.8em;">
        <p>Laboratorio Nacional de Virología, Ministerio de Salud | Instituto de Ciencias Sostenibles, Managua Nicaragua., {current_year}</p>
    </div>

</body>
</html>"""
    
    # Guardar reporte
    output_file = os.path.join(work_dir, "reporte_calidad.html")
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html)
    
    print(f"\n✅ Reporte generado: {output_file}")

if __name__ == "__main__":
    main()