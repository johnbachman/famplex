import io
import csv
import urllib
import requests
from collections import defaultdict
from indra.databases import hgnc_client

hgnc_fam_url = ('ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/csv/'
                'genefamily_db_tables/')
gene_fam_file = 'gene_has_family.csv'
family_file = 'family.csv'
hier_closure_file = 'hierarchy_closure.csv'
hier_file = 'hierarchy.csv'


def read_csv_from_ftp(fname):
    url = hgnc_fam_url + fname
    req = urllib.request.Request(url)
    res = urllib.request.urlopen(req)
    reader = csv.reader(io.TextIOWrapper(res))
    for row in reader:
        yield row


def _read_hgnc_family_genes():
    family_to_gene = defaultdict(list)
    gene_to_family = defaultdict(list)
    for gene_id, family_id in read_csv_from_ftp(gene_fam_file):
        family_to_gene[family_id].append(gene_id)
        gene_to_family[gene_id].append(family_id)
    return gene_to_family, family_to_gene


def _read_family_info():
    families = {}
    for idx, row in enumerate(read_csv_from_ftp(family_file)):
        if idx == 0:
            header = row
            continue
        families[row[0]] = {k: v for k, v in zip(header, row)}
    return families


def _read_hierarchy_info():
    children = defaultdict(list)
    for idx, (parent, child) in enumerate(read_csv_from_ftp(hier_file)):
        if idx == 0:
            continue
        children[parent].append(child)
    return children


families = _read_family_info()
children = _read_hierarchy_info()
gene_to_family, family_to_gene = _read_hgnc_family_genes()


def get_famplex_id(family):
    if family['abbreviation']:
        return family['abbreviation']
    else:
        return family['name'].replace(' ', '_').replace('-', '_')


def get_families_from_root(root_id):
    family_info = families[root_id]
    child_ids = children.get(root_id)
    famplex_id = get_famplex_id(family_info)
    if not child_ids:
        for gene in family_to_gene[root_id]:
            gene_name = hgnc_client.get_hgnc_name(gene)
            print('HGNC,%s,isa,FPLX,%s' % (gene_name, famplex_id))
    else:
        for child_id in child_ids:
            child_info = families[child_id]
            child_famplex_id = get_famplex_id(child_info)
            print('FPLX,%s,isa,FPLX,%s' % (child_famplex_id,
                                           famplex_id))
            get_families_from_root(child_id)


if __name__ == '__main__':
    get_families_from_root('177')
