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

python3 ncbi_dataset_builder_slow.py
"""

import sys
import traceback
import json
import sqlite3
import fastaparser

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


def parse_metadata(sqlite_con):
    """ Parses metadata and writes it to sqlite database
        Parsed values: accession, sgene range, location, collection date

        Args:
            sqlite_con - sqlite connection object 
    """

    print("Parsing meta data")
    metadata = []

    cursor = sqlite_con.cursor()

    with open(METADATA_JSON, "r") as meta_file:
        for line_json in meta_file:
            accession = None
            try:
                parser = NCBIMetadataParser(line_json)

                accession = parser.accession

                if not parser.is_annotated():
                    raise Exception("Not annotated")

                row = (
                    parser.accession,
                    parser.find_collection_date(),
                    *parser.find_sgene_range(),
                    *parser.find_location(),
                    None
                )
                metadata.append(row)

            except Exception as exc:
                print(f"{accession} {exc}")

                if DEBUG:
                    print(line_json)
                    traceback.print_exc()

    print(f"Metadata size: {len(metadata)} entries. Writing into sqlite")

    cursor.executemany(
        """INSERT INTO 
        ncbi ('accession', 'collection_date', 'sgene_begin', 'sgene_end', 'location', 'region', 'sgene') 
        VALUES (?, ?, ?, ?, ?, ?, ?)""", metadata)

    sqlite_con.commit()


def parse_genome(sqlite_con):
    """ Parses genome data, cuts SGENE and writes it into sqlite database

        Args:
            sqlite_con - sqlite connection object 
    """
    cursor = sqlite_con.cursor()

    import time
    print("Parsing genome data")
    counter = 0
    with open("genomic.fna", "r") as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
            t0 = t1 = time.process_time()

            try:
                cursor.execute(
                    "SELECT sgene_begin, sgene_end FROM ncbi WHERE accession = ?", (seq.id,))

                row = cursor.fetchone()
                if row is None:
                    raise Exception(f"{seq.id} is not found")

                begin = row[0] - 1  # FASTA count starts from 1
                end = row[1]
                sgene = seq.sequence_as_string()[begin:end]

                cursor.execute(
                    "UPDATE ncbi SET sgene = ? WHERE accession = ?", (sgene, seq.id,))

                if counter % 100 == 0:
                    print(f"Parsed genome {counter} samples")
                    sqlite_con.commit()
                counter += 1

            except Exception as exc:
                # print(exc)
                if DEBUG:
                    traceback.print_exc()

    sqlite_con.commit()
    print(f"Time taken {(time.process_time() - t0)} sec")


def create_database(cursor):
    print("(Re-)creating database")

    cursor.executescript("""
    DROP TABLE IF EXISTS ncbi;

    CREATE TABLE "ncbi" (
        "accession"	TEXT UNIQUE,
        "location"	TEXT DEFAULT NULL,
        "region"	TEXT DEFAULT NULL,
        "sgene"	TEXT DEFAULT NULL,
        "collection_date"	TEXT DEFAULT NULL,
        "sgene_begin"	INTEGER DEFAULT NULL,
        "sgene_end"	INTEGER DEFAULT NULL
    );

    CREATE INDEX "accession_UK" ON "ncbi" ("accession");
    """)


def main():

    with sqlite3.connect('ncbi.sqlite') as conn:
        create_database(conn)
        parse_metadata(conn)
        parse_genome(conn)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        traceback.print_exc()
