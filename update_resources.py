"""This script moves resources and exports from the top level of the repo into
the corresponding folders within the Python package.

Resource files and exports are only checked into version control at the top
level and are copied into the package by setup.py upon installation. If a user
clones this repo with the intention of contributing to FamPlex then resources
and exports can be copied directly into the package using this script. Running
this script after manually updating any of these files will make the updates
available within the package."""

import os
import shutil


if __name__ == '__main__':
    print('Copying resource files from top level into FamPlex package.')
    HERE = os.path.dirname(os.path.abspath(__file__))
    RESOURCES_PATH = os.path.join(HERE, 'famplex', 'resources')
    EXPORT_PATH = os.path.join(HERE, 'famplex', 'export')
    for resource in ['entities.csv', 'relations.csv', 'equivalences.csv',
                     'grounding_map.csv', 'gene_prefixes.csv',
                     'descriptions.csv']:
        shutil.copy(os.path.join(HERE, resource), RESOURCES_PATH)
    print('Copying exports from top level into FamPlex package.')
    for export in ['famplex.belns', 'famplex.obo', 'hgnc_symbol_map.csv',
                   'famplex_groundings.tsv']:
        shutil.copy(os.path.join(HERE, 'export', export), EXPORT_PATH)
