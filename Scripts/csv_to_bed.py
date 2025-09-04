import sys
import os
import csv
import subprocess
import tempfile

# --- Configuration ---
# Get the directory of the current script (e.g., /path/to/CONSENSO_D/Scripts)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Get the parent directory of SCRIPT_DIR (e.g., /path/to/CONSENSO_D)
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

REFERENCE_MAP = {
    "DENV_1": os.path.join(PROJECT_ROOT, "Ref_DENV/Reference_DV_1.fasta"),
    "DENV_2": os.path.join(PROJECT_ROOT, "Ref_DENV/Reference_DV_2.fasta"),
    "DENV_3": os.path.join(PROJECT_ROOT, "Ref_DENV/Reference_DV_3.fasta"),
    "DENV_4": os.path.join(PROJECT_ROOT, "Ref_DENV/Reference_DV_4.fasta"),
    "SARS_COV_2": os.path.join(PROJECT_ROOT, "Ref_DENV/REF_NC_045512_SARS_COV_2.fasta"),
    "RABV": os.path.join(PROJECT_ROOT, "Ref_DENV/RABV_Reference.fasta"),
    "N_RABV": os.path.join(PROJECT_ROOT, "Ref_DENV/Lisa_RABV_N.fasta"),
}

def run_command(command):
    """Runs a shell command and prints its output in real time."""
    print(f"\nExecuting: {' '.join(command)}")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding='utf-8')
    for line in iter(process.stdout.readline, ''):
        print(line.strip())
    process.wait()
    return process.returncode

def parse_sam_and_write_bed(sam_content, bed_path, primer_lengths):
    """Parses SAM content and writes to a BED file."""
    print(f"\nParsing alignment results and writing to {os.path.basename(bed_path)}...")
    count = 0
    with open(bed_path, 'w') as f_bed:
        for line in sam_content.strip().split('\n'):
            if line.startswith('@'):
                continue

            fields = line.split('\t')
            if len(fields) < 10:
                continue

            query_name = fields[0]
            flag = int(fields[1])
            ref_name = fields[2]
            pos = int(fields[3])

            if ref_name == '*':
                print(f"- WARNING: Primer '{query_name}' did not align. Skipping.")
                continue

            start = pos - 1
            end = start + primer_lengths[query_name]
            strand = '-' if (flag & 16) else '+'
            score = 0 # Default score

            f_bed.write(f"{ref_name}\t{start}\t{end}\t{query_name}\t{score}\t{strand}\n")
            count += 1
    print(f"Successfully wrote {count} primers to BED file.")
    return count

def main(csv_path, virus_type):
    """Main function to convert CSV primers to BED format."""
    print(f"Starting conversion for {os.path.basename(csv_path)} with reference '{virus_type}'")

    # 1. Validate inputs
    if not os.path.exists(csv_path):
        print(f"Error: Input CSV file not found: {csv_path}", file=sys.stderr)
        return 1

    ref_fasta = REFERENCE_MAP.get(virus_type)
    if not ref_fasta or not os.path.exists(ref_fasta):
        print(f"Error: Reference FASTA for '{virus_type}' not found at '{ref_fasta}'", file=sys.stderr)
        return 1

    # 2. Define output path
    output_bed_path = os.path.splitext(csv_path)[0] + ".bed"

    # 3. Read CSV and create temporary FASTA
    primers = {}
    try:
        with open(csv_path, 'r', encoding='utf-8-sig') as f_csv:
            reader = csv.reader(f_csv)
            header = next(reader) # Skip header
            for row in reader:
                if len(row) >= 2 and row[0] and row[1]:
                    # Sanitize primer name for FASTA header
                    primer_name = row[0].strip().replace(' ', '_')
                    sequence = row[1].strip().upper()
                    primers[primer_name] = sequence
    except Exception as e:
        print(f"Error reading or parsing CSV file: {e}", file=sys.stderr)
        return 1

    if not primers:
        print("Error: No valid primers found in the CSV file.", file=sys.stderr)
        return 1

    print(f"Found {len(primers)} primers in CSV file.")

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta") as tmp_fasta:
        for name, seq in primers.items():
            tmp_fasta.write(f">{name}\n{seq}\n")
        tmp_fasta_path = tmp_fasta.name

    # 4. Index reference and align
    print(f"\nIndexing reference genome: {os.path.basename(ref_fasta)}")
    if run_command(["bwa", "index", ref_fasta]) != 0:
        print("Error: BWA indexing failed.", file=sys.stderr)
        os.unlink(tmp_fasta_path)
        return 1

    print(f"\nAligning primers to reference...")
    align_process = subprocess.run(["bwa", "mem", ref_fasta, tmp_fasta_path], capture_output=True, text=True)
    sam_content = align_process.stdout

    # --- DEBUGGING: Print SAM output ---
    print("\n--- BWA Alignment (SAM) Output ---")
    print(sam_content if sam_content.strip() else "(SAM output is empty)")
    print("--- End of SAM Output ---")

    # 5. Cleanup temporary FASTA
    os.unlink(tmp_fasta_path)

    if align_process.returncode != 0:
        print(f"Error: BWA alignment failed.\n{align_process.stderr}", file=sys.stderr)
        return 1

    # 6. Parse SAM and write BED
    primer_lengths = {name: len(seq) for name, seq in primers.items()}
    primers_written = parse_sam_and_write_bed(sam_content, output_bed_path, primer_lengths)

    if primers_written == 0:
        print("\nError: Alignment completed, but no primers could be mapped to the reference genome.", file=sys.stderr)
        print("Please check that your primer sequences in the CSV are correct and correspond to the selected reference genome.", file=sys.stderr)
        return 1

    print("\nConversion process finished successfully!")
    return 0




if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python csv_to_bed.py <input_csv_path> <virus_type>", file=sys.stderr)
        sys.exit(1)

    exit_code = main(sys.argv[1], sys.argv[2])
    sys.exit(exit_code)
