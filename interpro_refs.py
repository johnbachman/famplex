import json
from indra.statements import Agent
from indra.databases import hgnc_client
from indra.tools.expand_families import Expander
from indra.preassembler.hierarchy_manager import hierarchies
from indra.databases import uniprot_client

def get_be_up_ids():
    with open('entities.csv', 'rt') as fh:
        entities = [line.strip() for line in fh.readlines()]
    be_agents = [Agent(be_id, db_refs={'BE': be_id})
                 for be_id in entities]
    ex = Expander(hierarchies)
    child_map = {}
    all_up_ids = set([])
    for be_agent in be_agents:
        children = ex.get_children(be_agent)
        children_up_ids = set([])
        for child_ns, child_id in children:
            assert child_ns == 'HGNC'
            hgnc_id = hgnc_client.get_hgnc_id(child_id)
            up_id = hgnc_client.get_uniprot_id(hgnc_id)
            children_up_ids.add(up_id)
        child_map[be_agent.name] = children_up_ids
        all_up_ids = all_up_ids.union(children_up_ids)
    return child_map, all_up_ids

def get_ip_families_for_be(cache_file=None):
    if cache_file is not None:
        with open(cache_file, 'rt') as f:
            ipfs = [line.strip() for line in f.readlines()]
        return ipfs
    # If not cached, read from file
    with open('protein2ipr.dat', 'rt') as f:
        counter = 0
        while True:
            if counter % 1000000 == 0:
                print("At line %d" % counter)
            line = f.readline()
            if not line:
                continue
            entries = line.strip().split('\t')
            up_id, ip_id = entries[0:2]
            if up_id in all_up_ids:
                print(entries)
                ip_families_for_be.add(ip_id)
            counter += 1

if __name__ == '__main__':
    child_map, all_up_ids = get_be_up_ids()
    ip_families_for_be = get_ip_families_for_be(
                                cache_file='ip_families_for_be.txt')
    ip_family_members = {}

    with open('protein2ipr.dat', 'rt') as f:
        f.seek(80*2700000)
        f.readline()
        counter = 0
        while True:
            if counter % 100000 == 0:
                print("At line %d" % counter)
            if counter > 1000000:
                break;
            line = f.readline()
            if not line:
                continue
            entries = line.strip().split('\t')
            up_id, ip_id, ip_name = entries[0:3]
            if ip_id in ip_families_for_be and uniprot_client.is_human(up_id):
                print(entries)
                ip_fam_info = ip_family_members.get(ip_id)
                if ip_fam_info is None:
                    ip_family_members[ip_id] = {'name': ip_name,
                                                'members': [up_id]}
                else:
                    ip_fam_info['members'].append(up_id)
            counter += 1
    with open('ip_family_members.json', 'wt') as f:
        json.dump(ip_family_members, f)

