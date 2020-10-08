"""Implements functions for loading resource files into datastructures."""
import csv
from typing import Dict, List, Optional, Tuple
from famplex.locations import ENTITIES_PATH, EQUIVALENCES_PATH, \
    GROUNDING_MAP_PATH, RELATIONS_PATH, GENE_PREFIXES_PATH, DESCRIPTIONS_PATH


__all__ = ['load_grounding_map', 'load_equivalences', 'load_entities',
           'load_relations', 'load_gene_prefixes', 'load_descriptions']


def _load_csv(filename):
    """Load famplex csv file as list of rows

    Parameters
    ----------
    filename : str

    Returns
    -------
    rows : list
    """
    with open(filename) as f:
        csvreader = csv.reader(f, delimiter=str(u','),
                               lineterminator='\r\n',
                               quoting=csv.QUOTE_MINIMAL,
                               quotechar=str(u'"'))
        rows = [row for row in csvreader]
    return rows


def _construct_grounding_map(rows):
    """Construct grounding map from rows in a grounding_map csv file

    Parameters
    ----------
    rows : list
        List of rows from a grounding map csv file. File should contain seven
        columns, the first of which is an agent text. The remaining columns
        contain namespace, id pairs, each pair occupying two columns. Some
        columns may be blank but in this case the row must be padded out with
        commas.

    Returns
    -------
    gmap : dict
        Dictionary mapping agent texts to INDRA style db_refs dicts. Each
        db_refs dict maps namespaces to ids.
    """
    gmap = {}
    for row in rows:
        text = row[0]
        db_refs = {'TEXT': text}
        db_refs.update({ns: id_ for ns, id_ in zip(row[1::2], row[2::2]) if ns})
        gmap[text] = db_refs if len(db_refs) > 1 else None
    return gmap


def load_grounding_map() -> Dict[str, Optional[Dict[str, str]]]:
    """Returns the FamPlex grounding map in dictionary form

    Returns
    -------
    dict
        A dictionary mapping agent texts to INDRA style db_refs dictionaries.
    """
    rows = _load_csv(GROUNDING_MAP_PATH)
    return _construct_grounding_map(rows)


def load_equivalences() -> List[Tuple[str, str, str]]:
    """Returns FamPlex equivalences as a list of rows.

    Returns
    -------
    list
        List of lists corresponding to rows from equivalences.csv. Each row
        contains three entries. A namespace, an ID, and a FamPlex ID. For
        example ['BEL', 'AMP Activated Protein Kinase Complex', 'AMPK'].
    """
    return _load_csv(EQUIVALENCES_PATH)


def load_entities() -> List[str]:
    """Returns list of FamPlex entities

    Returns
    -------
    list
        A list of all FamPlex unique IDs sorted in Unix standard sorted order.
    """
    rows = _load_csv(ENTITIES_PATH)
    return [row[0] for row in rows]


def load_relations() -> List[Tuple[str, str, str, str, str]]:
    """Returns FamPlex relations as a list of rows

    Returns
    -------
    list
        List of lists corresponding to rows in relations.csv. Each row has
        five columns of the form [namespace1, id1, relation, namespace2, id2].
        For example ['FPLX', 'AMPK_alpha', 'partof', 'FPLX', 'AMPK'].
    """
    return _load_csv(RELATIONS_PATH)


def load_gene_prefixes() -> List[Tuple[str, str, str]]:
    """Returns FamPlex gene prefixes as a list of rows

    Returns
    -------
    list
        List of lists corresponding to rows in gene_prefixes.csv. Each row has
        three columns [Pattern, Category, Notes].
    """
    return _load_csv(GENE_PREFIXES_PATH)


def load_descriptions() -> List[Tuple[str, str, str]]:
    """Returns FamPlex descriptions as a list of rows"""
    return _load_csv(DESCRIPTIONS_PATH)
