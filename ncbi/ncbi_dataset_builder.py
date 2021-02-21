"""
This script parses metadata, genome and protein files. Outputs csv file
containing following columns:

- accesion (str) - ID 
- sgene_begin (int) - whene sgene begins in full genome. **FASTA format starts from 1**
- sgene_end (int) - where sgene ends in full genome
- sgene_nucleotide (str) - sgene cut from full genome
- collection date (str) - when sample was collected. **Date format is NOT uniform**
- location (str) - country where sample was collected
- region (str) - continent where sample was collected
- protein_accession - protein ID
- sgene_protein (str) - sgene protein


Download data from https://www.ncbi.nlm.nih.gov/datasets/coronavirus/genomes/

python3 -m pip install -r requirements.txt

python3 ncbi_dataset_builder.py
"""

import sys
import time
import traceback
import json
import sqlite3
import pandas as pd

METADATA_JSON = "data_report.jsonl"
GENOMIC_NUCLEOTIDE_FASTA = "genomic.fna"
PROTEIN_FASTA = "protein.faa"

SGENE_SYNONIMS = ["spike glycoprotein", "surface glycoprotein",
                  "spike", "surface"]

DEBUG = False


class NCBIMetadataParser:

    def __init__(self, meta_str):
        self.meta = json.loads(meta_str)
        self.accession = self.meta["accession"]

    def find_sgene_annotation(self):
        """ Returns CDS dict of S gene. Raises Exception if not found
        """

        for gene in self.meta["annotation"]["genes"]:

            if gene["name"] != "S":
                continue

            for cds in gene["cds"]:

                # S gene has many names
                if cds["name"] not in SGENE_SYNONIMS:
                    continue

                return cds

        raise Exception(f"{self.accession} S Gene not found")

    def find_location(self):
        try:
            location = self.meta["location"]["geographicLocation"]
            region = self.meta["location"]["geographicRegion"]
            return location, region
        except KeyError:
            raise Exception(f"{self.accession} Missing 'location'")

    def find_collection_date(self):
        try:
            return self.meta["isolate"]["collectionDate"]
        except KeyError:
            raise Exception(f"{self.accession} Missing 'collection date'")

    def is_complete(self):
        if "completeness" not in self.meta:
            return False

        return self.meta["completeness"] == "COMPLETE"

    def is_annotated(self):
        return "annotation" in self.meta


def parse_metadata():
    """ Parses metadata and returns dict

        Returns
            data (dict) {accession: { protein accession, collection date
                                    sgene_begin, sgene_end, location, region },
                        ...}

    """

    print(f"READING METADATA_JSON {METADATA_JSON}")
    metadata = {}

    with open(METADATA_JSON, "r") as meta_file:
        for line_json in meta_file:
            try:
                parser = NCBIMetadataParser(line_json)

                if not parser.is_annotated():
                    raise Exception(f"{parser.accession} Not annotated")

                cds = parser.find_sgene_annotation()

                location, region = parser.find_location()

                metadata[parser.accession] = {
                    "protein_accession": cds["protein"]["accessionVersion"],
                    "collection_date": parser.find_collection_date(),
                    "sgene_begin": int(cds["nucleotide"]["range"][0]["begin"]),
                    "sgene_end": int(cds["nucleotide"]["range"][0]["end"]),
                    "location": location,
                    "region": region
                }

            except Exception as exc:
                print(f"{exc}")

                if DEBUG:
                    print(line_json)
                    traceback.print_exc()

    print(f"Metadata size: {len(metadata)} entries")

    return metadata


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

def main():

    t0 = time.process_time()

    dataset = parse_metadata()
    genomic = parse_genomic_nucleotide_fasta()
    protein = parse_protein_fasta()

    not_found = []  # accession ids not found in genomic dataset
    not_found_protein = []  # accession ids which protein is not found in protein dataset
    for accession, meta in dataset.items():

        if accession not in genomic:
            print(f"{accession} not found in genomic data")
            not_found.append(accession)
            continue

        if meta["protein_accession"] not in protein:
            print(f"{accession} | Protein {meta['protein_accession']} not found in protein data")
            not_found_protein.append(accession)
            continue

        begin = meta["sgene_begin"] - 1  # FASTA count starts from 1
        end = meta["sgene_end"]

        dataset[accession]["sgene_nucleotide"] = genomic[accession][begin:end]
        dataset[accession]["sgene_protein"] = protein[meta["protein_accession"]]

    # remove accessions from dataset not found in genomic dataset
    for accession in not_found:
        del dataset[accession]

    # remove accessions from dataset not found in protein dataset
    for accession in not_found_protein:
        del dataset[accession]

    df = pd.DataFrame.from_dict(dataset, orient='index')
    df.insert(0, 'accession', df.index)
    df = df.reset_index(drop=True)

    print(f"Dataframe shape: {df.shape}")

    df.to_csv("ncbi.csv", index=False)

    print(f"Time taken {int(time.process_time() - t0)} sec")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        traceback.print_exc()
