"""
 This script reads: metadata, genome, protein csv files, 
 joins three datasets together, extracts sgene and
 saves:
    - full dataset in csv
    - sgene in fasta
"""

import sys
import time
import json
import traceback
import pandas as pd

METADATA_CSV = "metadata.csv"
GENOMIC_NUCLEOTIDE_CSV = "genomic.csv"
PROTEIN_CSV = "protein.csv"


def cut_sgene(row):
    begin = int(row['sgene_begin']) - 1
    end = int(row['sgene_end'])
    sgene = row['genome'][begin:end]
    return sgene


def main():

    t0 = time.process_time()

    df = pd.read_csv(METADATA_CSV)
    print(f"Metadata shape {df.shape}")

    dfg = pd.read_csv(GENOMIC_NUCLEOTIDE_CSV)
    print(f"Genomic shape {dfg.shape}")

    dfp = pd.read_csv(PROTEIN_CSV)
    print(f"Protein shape {dfg.shape}")

    # join genomic data to metadata
    df = df.join(dfg.set_index('accession'), on='accession')

    # join protein data to metadata
    df = df.join(dfp.set_index('accession'), on='protein_accession')
    df.rename(columns={"protein": "sgene_protein",
                       "protein_desc": "sgene_protein_desc"}, inplace=True)

    df.dropna(inplace=True)

    with pd.option_context('mode.chained_assignment', None):
        df['sgene_nucleotide'] = df.apply(cut_sgene, axis=1)

    #
    # Save
    #

    # Full dataframe
    df.to_csv("ncbi_full.csv", index=False)

    # Sgenes in fasta format
    with open("ncbi_sgene_full.fasta", "w") as fasta_file:
        for _, row in df.iterrows():
            fasta_file.write(f">{row['accession']}\n")
            fasta_file.write(f"{row['sgene_nucleotide']}\n")
    
    #
    # Exclude sgene nucleotides containing letters
    #
    letters = ['R', 'M', "S", "B", "H", "N", "Y", "K", "W", "D", "V"]
    df_good = df[~df['sgene_nucleotide'].str.contains("|".join(letters))]

    # Good dataframe
    df_good.to_csv("ncbi_good.csv", index=False)

    # Good Sgenes in fasta format
    with open("ncbi_sgene_good.fasta", "w") as fasta_file:
        for _, row in df_good.iterrows():
            fasta_file.write(f">{row['accession']}\n")
            fasta_file.write(f"{row['sgene_nucleotide']}\n")

    print(f"Time taken {int(time.process_time() - t0)} sec")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        traceback.print_exc()
