"""
This script converts genome and protein files in fasta format to csv files
"""

import sys
import time
import traceback
import json
import pandas as pd

GENOMIC_NUCLEOTIDE_FASTA = "genomic.fna"
GENOMIC_NUCLEOTIDE_CSV = "genomic.csv"
PROTEIN_FASTA = "protein.faa"
PROTEIN_CSV = "protein.csv"

SGENE_SYNONIMS = ["spike glycoprotein", "surface glycoprotein", "spike", "surface"]

DEBUG = False

def parse_genomic_nucleotide_fasta():
    """ Reads genomic data and returns dict

        Returns:
            data (dict) - {accession: "FULL GENOME IN FASTA", 
                           ...}
    """

    print(f"Reading GENOMIC_NUCLEOTIDE_FASTA {GENOMIC_NUCLEOTIDE_FASTA}")
    data = {}

    with open(GENOMIC_NUCLEOTIDE_FASTA, "r") as fasta_file:

        sequence = []
        accession = None

        for line in fasta_file:
            if '>' in line:
                if sequence and accession:
                    data[accession] = "".join(sequence)

                accession = line.split(" ")[0][1:]
                sequence = []
            else:
                sequence.append(line.strip())

    return data


def parse_protein_fasta():
    """ Reads protein data and returns dict

        Returns:
            data (dict) - {protein_accession: "SGENE PROTEIN IN FASTA", 
                           ...}
    """

    print(f"Reading PROTEIN_FASTA {PROTEIN_FASTA}")
    data = {}

    with open(PROTEIN_FASTA, "r") as fasta_file:

        sequence = []
        sequence_started = True
        accession = None

        for line in fasta_file:
            if '>' in line:

                if sequence and accession:
                    data[accession] = "".join(sequence)

                sequence = []
                sequence_started = any(syn in line for syn in SGENE_SYNONIMS)
                accession = line.split(":", 1)[0][1:]
            else:
                if sequence_started:
                    sequence.append(line.strip())
    return data

def dict_to_df(data):
    df = pd.DataFrame.from_dict(data, orient='index', columns=['sequence'])
    df.insert(0, 'accession', df.index)
    df = df.reset_index(drop=True)
    return df

def main():

    t0 = time.process_time()

    genomic = parse_genomic_nucleotide_fasta()
    dfg = dict_to_df(genomic)
    dfg.to_csv(GENOMIC_NUCLEOTIDE_CSV, index=False)
    print(f"Wrote to {GENOMIC_NUCLEOTIDE_CSV}")

    protein = parse_protein_fasta()
    dfp = dict_to_df(protein)
    dfp.to_csv(PROTEIN_CSV, index=False)
    print(f"Wrote to {PROTEIN_CSV}")

    print(f"Time taken {int(time.process_time() - t0)} sec")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        traceback.print_exc()
