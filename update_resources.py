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
import warnings

# A warning will be issued if anything from famplex is imported while the
# resource files are unavailable. We supress this warning here so that it
# appear the first time this script is used to move the resources into the
# package.
warnings.simplefilter('ignore')
import famplex.locations as loc

if __name__ == '__main__':
    print('Copying resource files from top level into FamPlex package.')
    HERE = os.path.dirname(os.path.abspath(__file__))
    for resource in ['entities.csv', 'relations.csv', 'equivalences.csv',
                     'grounding_map.csv', 'gene_prefixes.csv',
                     'descriptions.csv']:
        shutil.copy(os.path.join(HERE, resource), loc.RESOURCES_PATH)
    print('Copying exports from top level into FamPlex package.')
    for export in ['famplex.belns', 'famplex.obo', 'hgnc_symbol_map.csv',
                   'famplex_groundings.tsv']:
        shutil.copy(os.path.join(HERE, 'export', export), loc.EXPORT_PATH)
