import subprocess

def ensamblar(serotype, path, tec_seq=None, seq_ref=None, primers=None, data=None):
    print(f"\nEnsamblando {serotype} en el directorio {path}")
    if tec_seq:
        print(f"Tecnología de Secuenciación: {tec_seq}")
    if seq_ref:
        print(f"Secuencia de Referencia: {seq_ref}")
    if primers:
        print(f"Primers: {primers}")
    if data:
        print(f"Datos adicionales: {data}")

    try:
        # Reemplaza este comando con la lógica de ensamblaje real
        subprocess.run(["echo", f"Ensamblando {serotype} en {path} con {tec_seq}, {seq_ref}, {primers}, {data}"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error durante la ejecución del comando: {e}")