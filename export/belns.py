#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""Output FamPlex as a BEL namespace. Requires the `bel_resources` package."""

import os

from bel_resources import write_namespace
from bel_resources.constants import NAMESPACE_DOMAIN_GENE

path_this = os.path.dirname(os.path.abspath(__file__))
entities_file = os.path.join(path_this, os.pardir, 'entities.csv')
output_file = os.path.join(path_this, 'famplex.belns')


def _get_entities():
    with open(entities_file, 'r') as fh:
        return {l.strip(): 'GRPC' for l in fh.readlines()}


def _write_namespace(values):
    with open(output_file, 'w') as file:
        write_namespace(
            namespace_name='FamPlex',
            namespace_keyword='FPLX',
            namespace_domain=NAMESPACE_DOMAIN_GENE,
            author_name='John Bachman and Ben Gyori',
            citation_name='FamPlex',
            values=values,
            namespace_description='FamPlex is a collection of resources for grounding biological entities from text '
                                  'and describing their hierarchical relationships.',
            namespace_query_url='http://identifiers.org/fplx/',
            author_copyright='CC0 1.0 Universal',
            citation_url='https://github.com/sorgerlab/famplex',
            case_sensitive=True,
            cacheable=True,
            file=file,
        )


if __name__ == '__main__':
    entities = _get_entities()
    _write_namespace(entities)
