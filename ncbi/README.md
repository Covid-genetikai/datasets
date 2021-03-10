Download **(SARS-CoV-2)** data from 
https://www.ncbi.nlm.nih.gov/datasets/coronavirus/genomes/    

Copy `genomic.fna, protein.faa, data_report.jsonl` to this directory

```
python3 -m pip install -r requirements.txt

# convert fasta to csv. Extract most important metadata from jsonl to csv
python3 convert.py

# build several datasets
python3 ncbi_dataset_builder.py
```

All information
- ncbi_full.csv - all available information (full genome, sgene, metadata)
- ncbi_sgene_full.fasta - all sgene nucleotides in fasta format

Only GOOD (sgenes without N, R, ... letters)
- ncbi_good.csv - only good available information (full genome, sgene, metadata)
- ncbi_sgene_good.fasta - only good sgene nucleotides in fasta format (do not contain N, R, ... letters)
- ncbi_sgene_good_unique.fasta - same as ncbi_sgene_good.fasta but contains only unique sgene_nucleotide

Columns:

- accesion (str) - ID 
- collection date (str) - when sample was collected. **Date format is NOT uniform**
- location (str) - country where sample was collected
- region (str) - continent where sample was collected
- genome (str) - full genome
- genome_desc (str) - full genome description from fast
- sgene_begin (int) - whene sgene begins in full genome. **FASTA format starts from 1**
- sgene_end (int) - where sgene ends in full genome
- sgene_nucleotide (str) - sgene cut from full genome
- protein_accession - protein ID
- sgene_protein (str) - sgene protein
- sgene_protein_desc (str) - sgene protein desctiption from fasta