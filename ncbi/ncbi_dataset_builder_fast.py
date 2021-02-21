"""
This script parses metadata and genome files and creates sqlite database
containing following columns:

accesion (ID) 
sgene_begin
sgene_end
sgene - sgene cut from full genome
collection date - when sample was collected
location - country where sample was collected
region - continent where sample was collected


Download data from https://www.ncbi.nlm.nih.gov/datasets/coronavirus/genomes/

python3 -m pip install -r requirements.txt

python3 ncbi_dataset_builder_fast.py
"""

import sys
import time
import traceback
import json
import sqlite3
import pandas as pd

METADATA_JSON = "data_report.jsonl"
GENOMIC_NUCLEOTIDE_FASTA = "genomic.fna"

SGENE_SYNONIMS = ["spike glycoprotein", "surface glycoprotein",
                  "spike", "surface"]

DEBUG = False


class NCBIMetadataParser:

    def __init__(self, meta_str):
        self.meta = json.loads(meta_str)
        self.accession = self.meta["accession"]

    def find_sgene_range(self):

        for gene in self.meta["annotation"]["genes"]:

            if gene["name"] != "S":
                continue

            for cds in gene["cds"]:

                # S gene has many names
                if cds["name"] not in SGENE_SYNONIMS:
                    continue

                begin = int(cds["nucleotide"]["range"][0]["begin"])
                end = int(cds["nucleotide"]["range"][0]["end"])

            return begin, end
        
        raise Exception("Range not found")

    def find_location(self):
        try:
            location = self.meta["location"]["geographicLocation"]
            region = self.meta["location"]["geographicRegion"]
            return location, region
        except KeyError:
            raise Exception("Missing 'location'")

    def find_collection_date(self):
        try:
            return self.meta["isolate"]["collectionDate"]
        except KeyError:
            raise Exception("Missing 'collection date'")

    def is_complete(self):
        if "completeness" not in self.meta:
            return False

        return self.meta["completeness"] == "COMPLETE"

    def is_annotated(self):
        return "annotation" in self.meta


def parse_metadata():
    """ Parses metadata and writes it to sqlite database
        Parsed values: accession, sgene range, location, collection date

        Returns
            data (dict)

    """

    print("Parsing meta data")
    metadata = {}

    with open(METADATA_JSON, "r") as meta_file:
        for line_json in meta_file:
            accession = None
            try:
                parser = NCBIMetadataParser(line_json)

                accession = parser.accession # For DEBUG (print in exception)

                if not parser.is_annotated():
                    raise Exception("Not annotated")

                begin, end = parser.find_sgene_range()
                location, region = parser.find_location()

                metadata[accession] = {
                    "collection_date": parser.find_collection_date(),
                    "sgene_begin": begin,
                    "sgene_end": end,
                    "location": location,
                    "region": region
                }

            except Exception as exc:
                print(f"{accession} {exc}")

                if DEBUG:
                    print(line_json)
                    traceback.print_exc()

    print(f"Metadata size: {len(metadata)} entries")

    return metadata
    

def parse_genomic():

    print("Reading genomic data")
    data = {}

    with open("genomic.fna", "r") as fasta_file:
        
        sequence = []
        seq_id = None

        for line in fasta_file:
            if '>' in line:
                if sequence and seq_id:
                    data[seq_id] = "".join(sequence)

                seq_id = line.split(" ")[0][1:]
                sequence = []            
            else:
                sequence.append(line.strip())
    
    return data


def main():

    t0 = time.process_time()

    dataset = parse_metadata()
    genomic = parse_genomic()

    not_found = [] # accesion ids not found in genomic dataset
    for accession, meta in dataset.items():

        if accession not in genomic:
            print(f"{accession} not found in genomic data")
            not_found.append(accession)
            continue

        begin = meta["sgene_begin"] - 1 # FASTA count starts from 1
        end = meta["sgene_end"]

        dataset[accession]["sgene"] = genomic[accession][begin:end]

    # remove accessions from dataset not found in genomic dataset
    for accession in not_found:
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
