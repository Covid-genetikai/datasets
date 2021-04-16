import logging
from pathlib import Path
from collections import namedtuple

import numpy as np
import pandas as pd

from anytree import Node, RenderTree, AsciiStyle
from anytree.exporter import DotExporter

logging.basicConfig(format='%(asctime)s %(message)s',
                    level=logging.INFO,
                    datefmt="%Y-%m-%d %H:%M:%S")

def main(distances_path, reference_name):

    MinPair = namedtuple("MinPair", 'distance target leaf')
    path = Path(distances_path)

    distances = pd.read_csv(path, index_col=0)
    distances = distances[distances.gt(0)]  # pavercia nulius i NaN

    # distances = distances[names]
    # print(distances.shape)
    # distances = distances.loc[names]
    # print(distances.shape)

    # Drop ref sgene row to prevent child referencing back to ref sgene
    distances.drop(reference_name, inplace=True, axis=0)

    rooted_tree = Node(reference_name)
    leaves = [rooted_tree]

    logging.info("Building tree")
    while distances.shape[0] > 0:
        min_pair = None

        for leaf in leaves: 

            # Suranda maziausia atstuma didesni uz nuli
            distance = distances[leaf.name].min()
            # ir tos sekos ID. Indexas yra seku ID
            target = distances[leaf.name].idxmin()

            if min_pair is None or min_pair.distance > distance:
                min_pair = MinPair(distance, target, leaf)

        leaves.append(Node(min_pair.target, parent=min_pair.leaf))

        # distances.drop(min_pair.target, inplace=True, axis=1)
        distances.drop(min_pair.target, inplace=True, axis=0)

        # checking whether node has 2 leaves
        if len(min_pair.leaf.children) == 2:
            leaves.remove(min_pair.leaf)

        if distances.shape[0] % 10 == 0:
            # Tarpine informacija
            logging.info(f"Remaining leaves {distances.shape[0]}")

    # print(RenderTree(rooted_tree, style=AsciiStyle()).by_attr())
    DotExporter(rooted_tree).to_dotfile(path.with_name("tree.dot"))
    logging.info("Done")

if __name__ == "__main__":
    reference_name = "2246_6593532f926e48cc68421ef20a33018c"
    distances_path = "resources/2021-04-16/levenshtein_aligned/distances.csv"
    distances_path = "resources/2021-04-16/levenshtein_not_aligned/distances.csv"
    distances_path = "resources/2021-04-16/RAxML_aligned/distances.csv"

    main(distances_path, reference_name)
