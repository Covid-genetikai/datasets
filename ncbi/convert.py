"""
This script converts 
    - fasta files (genome and protein)
    - jsonl file (metadata)
    to csv files
"""

import sys
import time
import traceback
import json
import pandas as pd

METADATA_JSON = "data_report.jsonl"
METADATA_CSV = "metadata.csv"
GENOMIC_NUCLEOTIDE_FASTA = "genomic.fna"
GENOMIC_NUCLEOTIDE_CSV = "genomic.csv"
PROTEIN_FASTA = "protein.faa"
PROTEIN_CSV = "protein.csv"

SGENE_SYNONIMS = ["spike glycoprotein",
                  "surface glycoprotein", "spike", "surface"]

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
        description = None
        for line in fasta_file:
            if '>' in line:
                if sequence and accession:
                    data[accession] = {
                        "genome": "".join(sequence),
                        "genome_desc": description
                    }

                parts = line.split(" ")
                accession = parts[0][1:]
                description = " ".join(parts[1:])

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
        description = None
        for line in fasta_file:
            if '>' in line:

                if sequence and accession:
                    data[accession] = {
                        "protein": "".join(sequence),
                        "protein_desc": description
                    }

                sequence = []
                sequence_started = any(syn in line for syn in SGENE_SYNONIMS)
                parts = line.split(":")
                accession = parts[0][1:]
                description = " ".join(parts[1:])
            else:
                if sequence_started:
                    sequence.append(line.strip())
    return data


def dict_to_df(data):
    df = pd.DataFrame.from_dict(data, orient='index')
    df.insert(0, 'accession', df.index)
    df = df.reset_index(drop=True)
    return df


def main():

    t0 = time.process_time()

    metadata = parse_metadata()
    dfm = dict_to_df(metadata)
    dfm.to_csv(METADATA_CSV, index=False)
    print(f"Wrote to {METADATA_CSV}")

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
