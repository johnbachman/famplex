import os
import sys
import csv
import json
import logging
import requests
from functools import lru_cache
from collections import defaultdict
from indra.statements import Agent
from indra.databases import hgnc_client
from indra.databases import uniprot_client
from indra.tools.expand_families import Expander
from indra.preassembler.hierarchy_manager import hierarchies


logger = logging.getLogger('reactome')


@lru_cache(10000)
def query_id(up_id):
    react_search_url = 'http://www.reactome.org/ContentService/search/query'
    params = {'query': up_id, 'cluster': 'true', 'species':'Homo sapiens'}
    headers = {'Accept': 'application/json'}
    res = requests.get(react_search_url, headers=headers, params=params)
    if not res.status_code == 200:
        return None
    json = res.json()
    results = json.get('results')
    if not results:
        print('No results for %s' % up_id)
        return None
    stable_ids = []
    for result in results:
        entries = result.get('entries')
        for entry in entries:
            stable_id = entry.get('stId')
            if not stable_id:
                continue
            name = entry.get('name')
            stable_ids.append(stable_id)
    return stable_ids


@lru_cache(100000)
def get_uniprot_id(reactome_id):
    react_url = 'http://www.reactome.org/ContentService/data/query/' \
                + reactome_id + '/referenceEntity'
    res = requests.get(react_url)
    if not res.status_code == 200:
        return None
    _, entry, entry_type = res.text.split('\t')
    if entry_type != 'ReferenceGeneProduct':
        return None
    id_entry = entry.split(' ')[0]
    db_ns, db_id = id_entry.split(':')
    if db_ns != 'UniProt':
        return None
    return db_id


@lru_cache(10000)
def get_participants(reactome_id):
    react_url = 'http://www.reactome.org/ContentService/data/event/' \
                + reactome_id + '/participatingReferenceEntities'
    headers = {'Accept': 'application/json'}
    res = requests.get(react_url, headers=headers)
    if not res.status_code == 200:
        return []
    json = res.json()
    up_ids = []
    for res in json:
        if not res.get('databaseName') == 'UniProt':
            continue
        up_id = res.get('identifier')
        if up_id is not None:
            up_ids.append(up_id)
    return up_ids


@lru_cache(10000)
def get_subunits(complex_id):
    react_url = 'http://www.reactome.org/ContentService/data/complex/' \
                + complex_id + '/subunits'
    headers = {'Accept': 'application/json'}
    res = requests.get(react_url, headers=headers)
    if not res.status_code == 200:
        return []
    json = res.json()
    up_ids = []
    for subunit in json:
        subunit_rx_id = subunit.get('stId')
        if not subunit_rx_id:
            continue
        up_id = get_uniprot_id(subunit_rx_id)
        if up_id:
            up_ids.append(up_id)
    return up_ids


@lru_cache(10000)
def get_parents(stable_id):
    react_data_url = 'http://www.reactome.org/ContentService/data/entity/' + \
                     stable_id + '/componentOf'
    headers = {'Accept': 'application/json'}
    res = requests.get(react_data_url, headers=headers)
    if not res.status_code == 200:
        return []
    json = res.json()
    names = []
    stable_ids = []
    schema_classes = []
    for parent_group in json:
        if not parent_group.get('type') in \
                            ['hasComponent', 'hasMember', 'hasCandidate']:
            continue
        names += parent_group.get('names')
        stable_ids += parent_group.get('stIds')
        schema_classes += parent_group.get('schemaClasses')
    parents_at_this_level = list(zip(names, stable_ids, schema_classes))
    parents_at_next_level_up = []
    for p_name, p_id, sc in parents_at_this_level:
        parents_at_next_level_up += get_parents(p_id)
    return parents_at_this_level + parents_at_next_level_up


@lru_cache(10000)
def get_all_parents(up_id):
    linked_stable_ids = query_id(up_id)
    if linked_stable_ids is None:
        return ([], [])
    parents = []
    for ls_id in linked_stable_ids:
        parents += get_parents(ls_id)
    sets = [tup for tup in parents if tup[2] != 'Complex']
    complexes = [tup for tup in parents if tup[2] == 'Complex']
    return sets, complexes


if __name__ == '__main__':

    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    """
    genes = []
    with open('../../bioentities/relations.csv', 'rt') as fh:
        csvreader = csv.reader(fh, delimiter=',', quotechar='"')
        for row in csvreader:
            if row[0] == 'HGNC':
                genes.append(row[1])
    with open('bioentities_genes.csv', 'wt') as f:
        for gene in genes:
            f.write('%s\n' % gene)

    gene_ids = []
    for hgnc_sym in genes:
        hgnc_id = hgnc_client.get_hgnc_id(hgnc_sym)
        up_id = hgnc_client.get_uniprot_id(hgnc_id)
        gene_ids.append((hgnc_sym, up_id))
    gene_ids = list(set(gene_ids))
    # Iterate over all genes and get all the sets (families)
    # For every bioentities family, get its child genes

    parents = []
    for hgnc_sym, up_id in gene_ids:
        parents = get_all_parents(up_id)
    """

    with open('entities.csv', 'rt') as fh:
        entities = [line.strip() for line in fh.readlines()]
    be_agents = [Agent(be_id, db_refs={'BE': be_id})
                 for be_id in entities]
    ex = Expander(hierarchies)
    child_map = {}
    rx_family_members = {}
    for be_agent in be_agents:
        print(be_agent.name)
        children = ex.get_children(be_agent)
        #child_map[be_agent.name] = children
        # For every child, get all parents
        rx_sets = set([])
        rx_complexes = set([])
        children_up_ids = set([])
        for child_ns, child_id in children:
            if child_ns == 'HGNC':
                hgnc_id = hgnc_client.get_hgnc_id(child_id)
                up_id = hgnc_client.get_uniprot_id(hgnc_id)
                children_up_ids.add(up_id)
            else:
                print("Unhandled NS: %s %s" % (child_ns, child_id))
                continue
            sets, complexes = get_all_parents(up_id)
            rx_sets = rx_sets.union(sets)
            rx_complexes = rx_complexes.union(complexes)
        # Next, for all of the sets, get their membership
        #matches = defaultdict(list)
        for rx_set in rx_sets.union(rx_complexes):
            name, rx_id, rx_type = rx_set
            # If we've already gotten info for this Reactome Set/Complex,
            # we can skip it
            if rx_id in rx_family_members:
                continue
            # Otherwise, get members
            if rx_type == 'DefinedSet' or rx_type == 'CandidateSet':
                rx_members = set(get_participants(rx_id))
            elif rx_type == 'Complex':
                rx_members = set(get_subunits(rx_id))
            else:
                print("Unrecognized type %s for %s, %s" %
                      (rx_type, name, rx_id))
                continue
            # Save info in dict
            rx_family_members[rx_id] = {'name': name, 'type': rx_type,
                                        'members': list(rx_members)}
            #print("Comparing %s with children %s (upids %s) to %s, rx_id %s "
            #      "with members %s" %
            #      (be_agent.name, children, children_up_ids, name, rx_id,
            #       rx_members))
            #if rx_members == children_up_ids:
            #    print("Found match for %s: %s" % (be_agent, name))
            #    matches[be_agent.name].append(rx_set)
    # Save to JSON
    with open('rx_family_members.json', 'wt') as f:
        json.dump(rx_family_members, f, indent=2)


