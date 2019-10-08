import io
import os
import csv
import urllib
import itertools
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
        return family['abbreviation'].strip().replace(', ', '_')
    else:
        replaces = {' ': '_', '-': '_', ',': ''}
        name = family['name'].strip()
        for k, v in replaces.items():
            name = name.replace(k, v)
        return name


def get_relations_from_root(root_id, relations=None):
    if relations is None:
        relations = []
    family_info = families[root_id]
    child_ids = children.get(root_id)
    famplex_id = get_famplex_id(family_info)
    if not child_ids:
        for gene in family_to_gene[root_id]:
            gene_name = hgnc_client.get_hgnc_name(gene)
            rel = ('HGNC', gene_name, 'isa', 'FPLX', famplex_id, root_id)
            relations.append(rel)
    else:
        for child_id in child_ids:
            child_info = families[child_id]
            child_famplex_id = get_famplex_id(child_info)
            rel = ('FPLX', child_famplex_id, 'isa', 'FPLX', famplex_id,
                   root_id)
            relations.append(rel)
            get_relations_from_root(child_id, relations)
    return relations


def add_relations_to_famplex(relations):
    rel_file = os.path.join(os.path.dirname(__file__), os.pardir,
                            'relations.csv')
    with open(rel_file, 'a') as fh:
        for rel in relations:
            fh.write(','.join(rel[:-1]) + '\n')


def add_entities_to_famplex(entities):
    ents_file = os.path.join(os.path.dirname(__file__), os.pardir,
                             'entities.csv')
    with open(ents_file, 'a') as fh:
        for ent in entities:
            fh.write('%s\n' % ent)


def add_equivalences(relations):
    hgnc_fam_ids = sorted(list(set(int(r[5]) for r in relations)))
    equivs = []
    for fid in hgnc_fam_ids:
        equivs.append(('HGNC_GROUP', str(fid),
                       get_famplex_id(families[str(fid)])))
    equivs_file = os.path.join(os.path.dirname(__file__), os.pardir,
                               'equivalences.csv')
    with open(equivs_file, 'a') as fh:
        for eq in equivs:
            fh.write('%s\n' % ','.join(eq))


def find_overlaps(relations):
    all_gene_names = {r[1]: r[4] for r in relations if r[0] == 'HGNC'}

    rel_file = os.path.join(os.path.dirname(__file__), os.pardir,
                            'relations.csv')
    covered_genes = set()
    covered_families = set()
    fam_members = defaultdict(list)
    hgnc_families = set()
    with open(rel_file, 'r') as fh:
        for sns, sid, rel, tns, tid in csv.reader(fh):
            if sns == 'HGNC' and tns == 'FPLX':
                fam_members[tid].append(sid)
            if sns == 'HGNC' and sid in all_gene_names:
                covered_genes.add(sid)
                print('%s covered already' % sid)
                covered_families.add(tid)
                hgnc_families.add(all_gene_names[sid])

    fplx_fam_members = {}
    for famplex_fam in covered_families:
        fplx_fam_members[famplex_fam] = set(fam_members[famplex_fam])

    fplx_fam_members = sorted(fplx_fam_members.items(),
                              key=lambda x: list(x[1])[0])

    hgnc_fam_members = {}
    for hgnc_fam in hgnc_families:
        hgnc_fam_members[hgnc_fam] = set(g for g, f in all_gene_names.items()
                                         if f == hgnc_fam)
    hgnc_fam_members = sorted(hgnc_fam_members.items(),
                              key=lambda x: list(x[1])[0])

    totally_redundant = set()
    for ff, hf in zip(fplx_fam_members, hgnc_fam_members):
        if set(ff[1]) == set(hf[1]):
            totally_redundant.add(hf[0])
            print('FamPlex %s and HGNC-derived %s are exactly the same.' %
                  (ff[0], hf[0]))
        else:
            print('FamPlex %s and HGNC-derived %s are overlapping.' %
                  (ff[0], hf[0]))
        print('Members of %s are: %s' % (ff[0], ','.join(sorted(ff[1]))))
        print('Members of %s are: %s' % (hf[0], ','.join(sorted(hf[1]))))
    return totally_redundant


if __name__ == '__main__':
    relations = get_relations_from_root('292')
    relations += get_relations_from_root('294')
    relations = sorted(list(set(relations)), key= lambda x: (x[4], x[1]))
    totally_redundant = find_overlaps(relations)
    relations = [r for r in relations if r[4] not in totally_redundant]
    entities = sorted(list(set(r[4] for r in relations)))
    add_relations_to_famplex(relations)
    add_entities_to_famplex(entities)
    add_equivalences(relations)
