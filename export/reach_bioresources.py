import os
import csv
from indra.databases import hgnc_client


path_this = os.path.dirname(os.path.abspath(__file__))
entities_file = os.path.join(path_this, os.pardir, 'entities.csv')
groundings_file = os.path.join(path_this, os.pardir, 'grounding_map.csv')


def _get_entities():
    with open(entities_file, 'r') as fh:
        return {l.strip(): 'FamilyOrComplex' for l in fh.readlines()}


def _get_groundings():
    groundings = []
    with open(groundings_file, 'r') as f:
        csvreader = csv.reader(f, delimiter=str(u','),
                               lineterminator='\r\n',
                               quoting=csv.QUOTE_MINIMAL,
                               quotechar=str(u'"'))
        for row in csvreader:
            txt = row[0]
            grounding_dict = {ns: id for ns, id in zip(row[1::2], row[2::2])}
            if 'FPLX' in grounding_dict:
                groundings.append((txt, grounding_dict['FPLX'], 'FPLX',
                                   'FamilyOrComplex'))
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
                                           ns.lower(), 'Simple_chemical'))
                        break
                else:
                    print(txt, grounding_dict)
        groundings = sorted(groundings)
        return groundings


if __name__ == '__main__':
    entities = _get_entities()
    entities_export = os.path.join(path_this, 'famplex.tsv')
    with open(entities_export, 'w') as fh:
        fh.write('\n'.join([('%s\t%s' % (e, t))
                            for e, t in entities.items()]))

    groundings = _get_groundings()
    groundings_export = os.path.join(path_this, 'famplex_groundings.tsv')
    with open(groundings_export, 'w') as fh:
        fh.write('\n'.join(['\t'.join(entries) for entries in groundings]))
