import pandas as pd
from pathlib import Path

FASTA_FILE = "ncbi_sgene_good_unique_aligned.aln"

def main():

    path = Path(FASTA_FILE)

    sequences = []
    with open(path, "r") as fasta_file:
        sequence = {}
        accession = None
        for line in fasta_file:
            if '>' in line:
                if sequence and accession:
                    sequences.append({
                        "accession": accession,
                        "sgene_nucleotide": "".join(sequence)
                    })

                parts = line.split(" ")
                accession = parts[0][1:].strip()
                sequence = []
            else:
                sequence.append(line.strip())
    
    df = pd.DataFrame.from_dict(sequences)
    df.to_csv(path.with_suffix(".csv"), index=False)

if  __name__ == "__main__":
    main()