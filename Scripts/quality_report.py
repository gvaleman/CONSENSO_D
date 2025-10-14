#!/usr/bin/env python3
"""
Script para generar reportes de calidad de secuenciación
Uso: python quality_report.py [directorio_muestras]
"""

import os
import sys
import glob
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot
import datetime

def parse_stats_file(stats_file):
    """Parsea archivo .stats de samtools y extrae métricas clave"""
    data = {}
    try:
        with open(stats_file, 'r') as f:
            for line in f:
                if line.startswith('SN'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        key = parts[1].replace(':', '').strip()
                        value = parts[2].strip()
                        # Intentar convertir a número
                        try:
                            if '.' in value:
                                data[key] = float(value)
                            else:
                                data[key] = int(value)
                        except:
                            data[key] = value
    except Exception as e:
        print(f"Error parseando {stats_file}: {e}")
    return data

def main():
    work_dir = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
    print(f"Procesando: {work_dir}")

    # === COBERTURA ===
    print("\n=== COBERTURA ===")
    coverage_files = glob.glob(os.path.join(work_dir, "**", "*.coverage"), recursive=True)
    coverage_plots_cards = []
    coverage_stats = []

    for coverage_file in coverage_files:
        sample_name = os.path.basename(coverage_file).replace(".coverage", "")

        try:
            df = pd.read_csv(coverage_file, sep="\t", header=None, names=['Ref', 'Pos', 'Cov'])
            df['Cov'] = pd.to_numeric(df['Cov'], errors='coerce').fillna(0)
            df['Pos'] = pd.to_numeric(df['Pos'], errors='coerce').fillna(0)

            print(f"  {sample_name}: {len(df)} posiciones, rango {df['Pos'].min()}-{df['Pos'].max()}")

            fig = go.Figure()
            fig.add_trace(go.Bar(
                x=df['Pos'],
                y=df['Cov'],
                name="Cobertura",
                marker_color='blue',
                width=1
            ))

            pos_min = df['Pos'].min()
            pos_max = df['Pos'].max()

            fig.update_layout(
                title=f"{sample_name}",
                xaxis_title="Posición Genómica",
                yaxis_title="Profundidad",
                height=350,
                template="plotly_white",
                xaxis=dict(
                    range=[pos_min - 50, pos_max + 50],
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='lightgray',
                    tickmode='linear',
                    dtick=max(1, (pos_max - pos_min) // 10)
                ),
                yaxis=dict(
                    range=[0, 1500],
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='lightgray'
                ),
                bargap=0,
                showlegend=False,
                margin=dict(l=50, r=30, t=50, b=50)
            )

            fig.add_hline(y=1000, line_dash="dash", line_color="red",
                         annotation_text="1000x", annotation_position="top right")

            plot_html = plot(fig, output_type='div', include_plotlyjs=False)
            coverage_plots_cards.append(plot_html)

            avg_cov = df['Cov'].mean()
            coverage_stats.append(avg_cov)
            print(f"✓ {sample_name}: {avg_cov:.0f}x promedio")

        except Exception as e:
            print(f"❌ Error en {sample_name}: {e}")

    # === READS Y Q-SCORES desde archivos .stats ===
    print("\n=== READS Y Q-SCORES (desde .stats) ===")
    stats_files = glob.glob(os.path.join(work_dir, "**", "*.stats"), recursive=True)
    
    reads_mapped = {}
    reads_unmapped = {}
    reads_total = {}
    qscore_data = {}

    for stats_file in stats_files:
        # Extraer nombre de muestra del archivo
        basename = os.path.basename(stats_file)
        # Remover prefijos comunes y extensión
        sample_name = basename.replace('Output_', '').replace('.sorted.bam.stats', '').replace('.stats', '')
        
        data = parse_stats_file(stats_file)
        
        if data:
            # Reads totales
            total = data.get('raw total sequences', 0)
            mapped = data.get('reads mapped', 0)
            unmapped = data.get('reads unmapped', 0)
            qscore = data.get('average quality', None)
            
            if total > 0:
                reads_total[sample_name] = total
                reads_mapped[sample_name] = mapped
                reads_unmapped[sample_name] = unmapped
                print(f"✓ {sample_name}: {total:,} reads totales ({mapped:,} mapeadas, {unmapped:,} no mapeadas)")
            
            if qscore is not None:
                qscore_data[sample_name] = float(qscore)
                print(f"  Q-score: {qscore}")

    # Gráfico de reads con barras apiladas (mapeadas vs no mapeadas)
    reads_plot = ""
    if reads_total:
        samples = list(reads_total.keys())
        
        fig = go.Figure()
        
        # Barras apiladas
        fig.add_trace(go.Bar(
            name='Mapeadas',
            x=samples,
            y=[reads_mapped.get(s, 0) for s in samples],
            marker_color='steelblue',
            text=[f"{reads_mapped.get(s, 0):,}" for s in samples],
            textposition='inside',
            hovertemplate='<b>%{x}</b><br>Mapeadas: %{y:,}<extra></extra>'
        ))
        
        fig.add_trace(go.Bar(
            name='No Mapeadas',
            x=samples,
            y=[reads_unmapped.get(s, 0) for s in samples],
            marker_color='lightcoral',
            text=[f"{reads_unmapped.get(s, 0):,}" for s in samples],
            textposition='inside',
            hovertemplate='<b>%{x}</b><br>No Mapeadas: %{y:,}<extra></extra>'
        ))
        
        fig.update_layout(
            title="Reads Totales por Muestra",
            xaxis_title="Muestra",
            yaxis_title="Número de Reads",
            template="plotly_white",
            barmode='stack',
            xaxis={'tickangle': 45},
            height=350,
            showlegend=True,
            legend=dict(x=0.01, y=0.99, bgcolor='rgba(255,255,255,0.8)')
        )
        
        reads_plot = plot(fig, output_type='div', include_plotlyjs=False)

    # Gráfico de Q-scores
    qscore_plot = ""
    if qscore_data:
        fig = go.Figure(data=[go.Bar(
            x=list(qscore_data.keys()),
            y=list(qscore_data.values()),
            marker_color='green',
            text=[f"{v:.1f}" for v in qscore_data.values()],
            textposition='outside'
        )])
        
        fig.update_layout(
            title="Q-score Promedio por Muestra",
            xaxis_title="Muestra",
            yaxis_title="Q-score",
            template="plotly_white",
            xaxis={'tickangle': 45},
            height=350
        )
        
        qscore_plot = plot(fig, output_type='div', include_plotlyjs=False)

    # Estadísticas finales
    avg_coverage = sum(coverage_stats) / len(coverage_stats) if coverage_stats else 0
    samples_above_1000x = sum(1 for cov in coverage_stats if cov >= 1000)
    total_samples = len(coverage_files)
    current_year = datetime.datetime.now().year

    # Generar HTML con grid de 2 columnas para cobertura
    coverage_grid_html = ""
    for plot_html in coverage_plots_cards:
        coverage_grid_html += f'<div class="coverage-card">{plot_html}</div>\n'

    # HTML
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
        .coverage-section {{ margin: 20px 0; background: white; padding: 20px; border-radius: 10px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        .coverage-grid {{ display: grid; grid-template-columns: repeat(2, 1fr); gap: 20px; margin-top: 20px; }}
        .coverage-card {{ background: #f9f9f9; padding: 15px; border-radius: 8px; border: 1px solid #e0e0e0; }}
        @media (max-width: 1024px) {{
            .charts {{ grid-template-columns: 1fr; }}
            .coverage-grid {{ grid-template-columns: 1fr; }}
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Reporte de Calidad de Secuenciación y Ensamblaje de Genomas</h1>
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
            <div class="stat-number">{len(qscore_data)}</div>
            <div>Con Datos Q-score</div>
        </div>
    </div>

    <div class="charts">
        <div class="chart">
            <h3>Reads Totales por Muestra</h3>
            {reads_plot if reads_plot else '<p style="text-align: center; color: #666;">No hay datos disponibles</p>'}
        </div>
        <div class="chart">
            <h3>Q-scores Promedio</h3>
            {qscore_plot if qscore_plot else '<p style="text-align: center; color: #666;">No hay datos disponibles</p>'}
        </div>
    </div>

    <div class="coverage-section">
        <h3>Cobertura Genómica por Muestra</h3>
        <p style="color: #666;">La línea roja indica 1000x de cobertura (óptima para análisis viral)</p>
        <div class="coverage-grid">
            {coverage_grid_html if coverage_grid_html else '<p style="text-align: center; color: #666;">No hay datos disponibles</p>'}
        </div>
    </div>

    <div style="text-align: center; margin-top: 20px; color: #666;">
        <p>Archivos procesados: {len(coverage_files)} coverage | {len(stats_files)} stats</p>
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
    print(f"   - {total_samples} muestras procesadas")
    print(f"   - {len(reads_total)} con datos de reads")
    print(f"   - {len(qscore_data)} con datos de Q-score")

if __name__ == "__main__":
    main()