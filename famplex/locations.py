"""This module contains paths to famplexes resource and export files."""
import os

FPLX_PATH = os.path.dirname(os.path.abspath(__file__))
RESOURCES_PATH = os.path.join(FPLX_PATH, 'resources')
EXPORT_PATH = os.path.join(FPLX_PATH, 'export')

# Paths to resources
ENTITIES_PATH = os.path.join(RESOURCES_PATH, 'entities.csv')
RELATIONS_PATH = os.path.join(RESOURCES_PATH, 'relations.csv')
EQUIVALENCES_PATH = os.path.join(RESOURCES_PATH, 'equivalences.csv')
GROUNDING_MAP_PATH = os.path.join(RESOURCES_PATH, 'grounding_map.csv')
GENE_PREFIXES_PATH = os.path.join(RESOURCES_PATH, 'gene_prefixes.csv')
DESCRIPTIONS_PATH = os.path.join(RESOURCES_PATH, 'descriptions.csv')

# Paths to exports
BELNS_PATH = os.path.join(EXPORT_PATH, 'famplex.belns')
OBO_PATH = os.path.join(EXPORT_PATH, 'famplex.obo')
HGNC_IDS_PATH = os.path.join(EXPORT_PATH, 'hgnc_symbol_map.csv')
GROUNDINGS_PATH = os.path.join(EXPORT_PATH, 'famplex_groundings.tsv')
