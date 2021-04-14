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

Output used in next step: `ncbi_sgene_good_unique.fasta`

```
# align with small amount (50 sequences)
head -n 100 ncbi_sgene_good_unique.fasta > ncbi_sgene_good_unique_100.fasta
muscle -in ncbi_sgene_good_unique_100.fasta -clwout ncbi_sgene_good_unique_100.aln
```

```
# align all sequences
muscle -in ncbi_sgene_good_unique.fasta -clwout ncbi_sgene_good_unique.aln
```

Output used in next step: `ncbi_sgene_good_unique.aln`

```
# build tree (TO Do UPDATE)
https://github.com/stamatak/standard-RAxML
./standard-RAxML/raxmlHPC-PTHREADS-SSE3 -T 8 -f a -x 860647 -p 860647 -N 100 -m GTRGAMMA -O -n sgene_good_unique.tre -s $(PWD)/data/sgene_good_unique_modified.aln -w $(PWD)
```

Output used in next step: `RAxML_bestTree.sgene_good_unique.tre`

```
# Build binary tree 
https://github.com/Covid-genetikai/aligners/blob/main/src/binaryTreeGen_multiproc.py

python3 -m pip install anytree biopython
python3 binaryTreeGen_multiproc.py 
```

Output used in next step: `tree.dot`

```
https://www.graphviz.org/doc/info/output.html
sudo apt-get install graphviz 

# Convert tree.dot to tree.json
dot -Txdot_json -o tree.json tree.dot

# Convert tree.dot to svg (works also with png, jpg)
dot -Tsvg tree.dot -o tree.svg
```

All information
- ncbi_full.csv - all available information (full genome, sgene, metadata)
- ncbi_sgene_full.fasta - all sgene nucleotides in fasta format

Only GOOD (sgenes without N, R, ... letters)
- ncbi_good.csv - only good available information (full genome, sgene, metadata)
- ncbi_sgene_good.fasta - only good sgene nucleotides in fasta format (do not contain N, R, ... letters)
- ncbi_sgene_good_unique.fasta - same as ncbi_sgene_good.fasta but contains only unique sgene_nucleotide
- ncbi_sgene_good_unique.aln - aligned sgenes 
- tree.dot - binary tree in dot format


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