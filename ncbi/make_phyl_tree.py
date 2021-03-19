from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO

import sys
sys.setrecursionlimit(10000)

import time
from pathlib import Path

def main():

    path = Path("good_unique2")
    # path = Path("good_unique_small")

    # Read aligned sequences
    aln = AlignIO.read(path.with_suffix(".aln"), 'fasta')

    # Tree making magic
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    # Save tree in ascii
    Phylo.draw_ascii(tree, open(path.with_suffix(".ascii"), "w"))

    # Save tree in dnd
    Phylo.write(tree, path.with_suffix(".dnd"), format='newick')

    # Show tree image
    Phylo.draw(tree)

if __name__ == "__main__":
    
    t0 = time.process_time()

    try:
        main()
    except:
        raise
    finally:
        print(f"Time taken {int(time.process_time() - t0)} sec")