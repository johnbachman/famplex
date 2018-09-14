from __future__ import print_function, unicode_literals
import csv
import sys
from os.path import abspath, join, dirname

data_path = join(dirname(abspath(__file__)), '../..')


def load_prefixes(filename):
    return {}


def load_entity_list(filename):
    rows = load_csv(filename)
    entities = [row[0] for row in rows]
    return entities


def load_relations(filename):
    relationships = []
    rows = load_csv(filename)
    for row in rows:
        relationships.append(((row[0], row[1]), row[2], (row[3], row[4])))
    return relationships


def load_equivalences(filename):
    equivalences = []
    rows = load_csv(filename)
    for row in rows:
        equivalences.append((row[0], row[1], row[2]))
    return equivalences


def load_grounding_map(filename):
    gm_rows = load_csv(filename)
    g_map = {}
    by_term = {}
    by_fplx = {}
    for row in gm_rows:
        gm_tuples.append(tuple(row))
        key = row[0]
        db_refs = {'TEXT': key}
        keys = [entry for entry in row[1::2] if entry != '']
        values = [entry for entry in row[2::2] if entry != '']
        if len(keys) != len(values):
            print('ERROR: Mismatched keys and values in row %s' % str(row))
            continue
        else:
            db_refs.update(dict(zip(keys, values)))
            if len(db_refs.keys()) > 1:
                g_map[key] = db_refs
            else:
                g_map[key] = None
    return g_map, tuple(gm_tuples)


def load_csv(filename):
    # Python 3 version
    if sys.version_info[0] >= 3:
        # Open the file in text mode with given encoding
        # Set newline arg to '' (see https://docs.python.org/3/library/csv.html)
        with open(filename, 'r', newline='', encoding='utf8') as f:
            csv_reader = csv.reader(f, delimiter=',', quotechar='"',
                                    lineterminator='\n')
            for row in csv_reader:
                yield row
    # Python 2 version
    else:
        # Open the file in binary mode
        with open(filename, 'rb') as f:
            # Next, get the csv reader, passing delimiter and quotechar as
            # bytestrings rather than unicode
            csv_reader = csv.reader(f, delimiter=','.encode('utf8'),
                                    quotechar='"'.encode('utf8'))
            # Iterate over the file and decode each string into unicode
            for row in csv_reader:
                yield [cell.decode(encoding) for cell in row]


entities = load_entity_list(join(data_path, 'entities.csv'))
grounding_map = load_grounding_map(join(data_path, 'grounding_map.csv'))
relations = load_relations(join(data_path, 'relations.csv'))
equivalences = load_equivalences(join(data_path, 'equivalences.csv'))
prefixes = load_prefixes(join(data_path, 'gene_prefixes.csv'))

