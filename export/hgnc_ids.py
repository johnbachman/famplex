"""Several files in this resource use HGNC gene symbols to identify individual
genes. However, the symbols assigned to HGNC IDs can change over time, and
therefore previously curated symbols can become invalid. This script
generates a mapping of current (i.e. at the time of running the script)
mappings of HGNC IDs to symbols so that the assumptions about the identity
of the genes in the various tables can be traced."""

import os
import csv
from indra.databases import hgnc_client

if __name__ == '__main__':
    path_this = os.path.dirname(os.path.abspath(__file__))
    hgnc_symbols = set()
    # Gather all HGNC symbols from relations.csv
    relations_file = os.path.join(path_this, os.pardir, 'relations.csv')
    with open(relations_file, 'r') as f:
        csvreader = csv.reader(f, delimiter=str(u','), lineterminator='\r\n',
                               quoting=csv.QUOTE_MINIMAL,
                               quotechar=str(u'"'))
        for row in csvreader:
            ns1, id1, rel, ns2, id2 = row
            if ns1 == 'HGNC':
                hgnc_symbols.add(id1)
            if ns2 == 'HGNC':
                hgnc_symbols.add(id2)

    # Gather all HGNC symbols from grounding_map.csv
    gm_file = os.path.join(path_this, os.pardir, 'grounding_map.csv')
    with open(gm_file, 'r') as f:
        csvreader = csv.reader(f, delimiter=str(u','), lineterminator='\r\n',
                               quoting=csv.QUOTE_MINIMAL,
                               quotechar=str(u'"'))
        for row in csvreader:
            namespaces = row[1::2]
            ids = row[2::2]
            for ns, id in zip(namespaces, ids):
                if ns == 'HGNC':
                    hgnc_symbols.add(id)

    # Create output file
    out_file = os.path.join(path_this, 'hgnc_symbol_map.csv')
    with open(out_file, 'w') as fh:
        for hgnc_symbol in sorted(list(hgnc_symbols)):
            hgnc_id = hgnc_client.get_hgnc_id(hgnc_symbol)
            fh.write('%s,%s\r\n' % (hgnc_symbol, hgnc_id))

