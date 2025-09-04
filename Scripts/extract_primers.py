import sys
import os

def parse_fasta(fasta_file):
    sequences = {}
    with open(fasta_file, 'r') as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    sequences[header] = ''.join(seq)
                header = line.split()[0][1:]
                seq = []
            else:
                seq.append(line)
        if header:
            sequences[header] = ''.join(seq)
    return sequences

def main():
    if len(sys.argv) != 4:
        print("Usage: python extract_primers.py <reference.fasta> <primers.bed> <output.fasta>")
        sys.exit(1)

    ref_fasta_file = sys.argv[1]
    primers_bed_file = sys.argv[2]
    output_fasta_file = sys.argv[3]

    if not os.path.exists(ref_fasta_file):
        print(f"Error: Reference FASTA file not found at {ref_fasta_file}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(primers_bed_file):
        print(f"Error: Primers BED file not found at {primers_bed_file}", file=sys.stderr)
        sys.exit(1)

    sequences = parse_fasta(ref_fasta_file)

    with open(primers_bed_file, 'r') as bed_f, open(output_fasta_file, 'w') as out_f:
        for line in bed_f:
            if line.strip() == '' or line.startswith('#'):
                continue
            try:
                fields = line.strip().split()
                chrom, start, end, name, score, strand = fields[:6]
                start, end = int(start), int(end)

                if chrom not in sequences:
                    print(f"Warning: Chromosome '{chrom}' from BED file not found in FASTA file. Skipping primer {name}.", file=sys.stderr)
                    continue

                ref_seq = sequences[chrom]
                primer_seq = ref_seq[start:end]

                if strand == '-':
                    complement = str.maketrans('ATCGN', 'TAGCN')
                    primer_seq = primer_seq.translate(complement)[::-1]

                out_f.write(f">{name}\n{primer_seq}\n")
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse line in BED file: {line.strip()}. Error: {e}. Skipping.", file=sys.stderr)
                continue

if __name__ == "__main__":
    main()
