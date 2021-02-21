Download **(SARS-CoV-2)** data from 
https://www.ncbi.nlm.nih.gov/datasets/coronavirus/genomes/    

```
python3 -m pip install -r requirements.txt

python3 ncbi_dataset_builder.py
```

Output:

Script parses metadata and genome files and creates output file containing following columns:
- accesion (str) - ID 
- sgene_begin (int) - whene sgene begins in full genome. **FASTA format starts from 1**
- sgene_end (int) - where sgene ends in full genome
- sgene_nucleotide (str) - sgene cut from full genome
- collection date (str) - when sample was collected. **Date format is NOT uniform**
- location (str) - country where sample was collected
- region (str) - continent where sample was collected
- protein_accession - protein ID
- sgene_protein (str) - sgene protein