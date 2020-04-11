import os
import csv
from collections import Counter
from indra.databases import hgnc_client


path_this = os.path.dirname(os.path.abspath(__file__))
groundings_file = os.path.join(path_this, os.pardir, 'grounding_map.csv')


def get_groundings():
    groundings = []
    text_appearances = []
    with open(groundings_file, 'r') as f:
        csvreader = csv.reader(f, delimiter=str(u','),
                               lineterminator='\r\n',
                               quoting=csv.QUOTE_MINIMAL,
                               quotechar=str(u'"'))
        for row in csvreader:
            txt = row[0]
            text_appearances.append(txt)
            grounding_dict = {ns: id for ns, id in zip(row[1::2], row[2::2])}
            if 'FPLX' in grounding_dict:
                groundings.append((txt, grounding_dict['FPLX'], 'fplx',
                                   #'FamilyOrComplex'))
                                   'Family'))
            elif 'UP' in grounding_dict:
                groundings.append((txt, grounding_dict['UP'], 'uniprot',
                                   'Gene_or_gene_product'))
            elif 'HGNC' in grounding_dict:
                up_id = hgnc_client.get_uniprot_id(grounding_dict['HGNC'])
                if up_id:
                    groundings.append((txt, up_id, 'uniprot',
                                      'Gene_or_gene_product'))
                else:
                    groundings.append((txt, grounding_dict['HGNC'], 'hgnc',
                                       'Gene_or_gene_product'))
            elif 'IP' in grounding_dict:
                groundings.append((txt, grounding_dict['IP'],
                                   'interpro', 'FamilyOrComplex'))
            else:
                mappings = {'CHEBI': 'Simple_chemical',
                            'PUBCHEM': 'Simple_chemical',
                            'CHEMBL': 'Simple_chemical',
                            'HMDB': 'Simple_chemical',
                            'GO': 'BioProcess',
                            'MESH': 'BioProcess',
                            'NCIT': 'BioProcess'}
                for ns, type in mappings.items():
                    if ns in grounding_dict:
                        groundings.append((txt, grounding_dict[ns],
                                           ns.lower(), type))
                        break
                else:
                    print(txt, grounding_dict)
    cnt = Counter(text_appearances)
    ambiguous_txts = {t for t, c in cnt.items() if c >= 2}
    groundings = [g for g in sorted(groundings) if g[0] not in ambiguous_txts]
    return groundings


if __name__ == '__main__':
    groundings = get_groundings()
    grounding_export = os.path.join(path_this,
        'famplex_groundings.tsv')
    with open(grounding_export, 'w') as fh:
        fh.write('\n'.join(['\t'.join(entries) for entries in groundings]))