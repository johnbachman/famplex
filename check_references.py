from __future__ import print_function, unicode_literals
import csv
from collections import Counter

def load_csv(filename):
    with open(filename) as f:
        csvreader = csv.reader(f, delimiter=',', quotechar='"')
        rows = [row for row in csvreader]
    return rows

def load_grounding_map(filename):
    gm_rows = load_csv(filename)
    g_map = {}
    for row in gm_rows:
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
    return g_map

def load_entity_list(filename):
    with open(filename) as f:
        csvreader = csv.reader(f, delimiter=',', quotechar='"')
        entities = [row[0] for row in csvreader]
    return entities

def load_relationships(filename):
    relationships = []
    with open(filename) as f:
        csvreader = csv.reader(f, delimiter=',', quotechar='"')
        for row in csvreader:
            relationships.append(((row[0], row[1]), row[2], (row[3], row[4])))
    return relationships

def update_id_prefixes(filename):
    gm_rows = load_csv(filename)
    updated_rows = []
    for row in gm_rows:
        key = row[0]
        keys = [entry for entry in row[1::2]]
        values = [entry for entry in row[2::2]]
        if 'GO' in keys:
            go_ix = keys.index('GO')
            values[go_ix] = 'GO:%s' % values[go_ix]
        if 'CHEBI' in keys:
            chebi_ix = keys.index('CHEBI')
            values[chebi_ix] = 'CHEBI:%s' % values[chebi_ix]
        if 'CHEMBL' in keys:
            chembl_ix = keys.index('CHEMBL')
            values[chembl_ix] = 'CHEMBL%s' % values[chembl_ix]
        updated_row = [key]
        for pair in zip(keys, values):
            updated_row += pair
        updated_rows.append(updated_row)
    return updated_rows

def pubchem_and_chebi(db_refs):
    pubchem_id = db_refs.get('PUBCHEM')
    chebi_id = db_refs.get('CHEBI')
    if pubchem_id and not chebi_id:
        return 'chebi_missing'
    if chebi_id and not pubchem_id:
        return 'pubchem_missing'
    return None

if __name__ == '__main__':
    # Check the entity list for duplicates
    entities = load_entity_list('entities.csv')
    ent_counter = Counter(entities)
    print("-- Checking for duplicate entities --")
    found_duplicates = False
    for ent, freq in ent_counter.items():
        if freq > 1:
            print("ERROR: Duplicate entries for %s in entity list." % ent)
            found_duplicates = True
    if not found_duplicates:
        print("OK! No duplicates found.")

    print()
    print("-- Checking for undeclared Bioentities IDs in grounding map --")
    # Load the grounding map
    gm = load_grounding_map('grounding_map.csv')
    # Look through grounding map and find all instances with an 'BE' db
    entities_missing_gm = []
    for text, db_refs in gm.items():
        if db_refs is not None:
            for db_key, db_id in db_refs.items():
                if db_key == 'BE' and db_id not in entities:
                    entities_missing_gm.append(db_id)
                    print("ERROR: ID %s referenced in grounding map "
                          "is not in entities list." % db_id)

    print()
    print("-- Checking for CHEBI/PUBCHEM IDs--")
    chebi_id_missing = []
    pubchem_id_missing = []
    for text, db_refs in gm.items():
        if db_refs is not None:
            p_and_c = pubchem_and_chebi(db_refs)
            if p_and_c == 'chebi_missing':
                chebi_id_missing.append(db_refs['PUBCHEM'])
                print("WARNING: %s has PUBCHEM ID but no CHEBI ID." % text)
            if p_and_c == 'pubchem_missing':
                pubchem_id_missing.append(db_refs['CHEBI'])
                print("WARNING: %s has CHEBI ID but no PUBCHEM ID." % text)

    print()
    print("-- Checking for undeclared Bioentities IDs in relationships file --")
    # Load the relationships
    relationships = load_relationships('relations.csv')
    # Check the relationships for consistency with entities
    entities_missing_rel = []
    for subj, rel, obj in relationships:
        for term in (subj, obj):
            term_ns = term[0]
            term_id = term[1]
            if term_ns == 'BE' and term_id not in entities:
                entities_missing_rel.append(term_id)
                print("ERROR: ID %s referenced in relations "
                      "is not in entities list." % term_id)
    print()
    print("-- Checking for valid namespaces in relations --")
    for ix, (subj, rel, obj) in enumerate(relationships):
        for term in (subj, obj):
            term_ns = term[0]
            if term_ns not in ('BE', 'HGNC', 'UP'):
                print("ERROR: row %d: Invalid namespace in relations.csv: %s" %
                      (ix+1, term_ns))

    # This check requires the indra package
    try:
        from indra.databases import hgnc_client
        print()
        print("-- Checking for invalid HGNC IDs in relationships file --")
        for subj, rel, obj in relationships:
            for term in (subj, obj):
                term_ns = term[0]
                term_id = term[1]
                if term_ns == 'HGNC':
                    hgnc_id = hgnc_client.get_hgnc_id(term_id)
                    if not hgnc_id:
                        print("ERROR: ID %s referenced in relations is "
                              "not a valid HGNC ID." % term_id)
    except ImportError:
        pass

    # This check requires the indra package
    try:
        from indra.databases import hgnc_client
        print()
        print("-- Checking for invalid HGNC IDs in grounding map --")
        for text, db_refs in gm.items():
            if db_refs is not None:
                for db_key, db_id in db_refs.items():
                    if db_key == 'HGNC':
                        hgnc_id = hgnc_client.get_hgnc_id(db_id)
                        if not hgnc_id:
                            print("ERROR: ID %s in grounding map is "
                                   "not a valid HGNC ID." % db_id)
    except ImportError:
        pass

    # This check requires a ChEBI resource file to be available. You
    # can obtain it from here: ftp://ftp.ebi.ac.uk/pub/databases/
    #                          chebi/Flat_file_tab_delimited/compounds.tsv.gz
    try:
        with open('chebi_compounds.tsv', 'rt') as fh:
            chebi_ids = [lin.split('\t')[2] for lin in fh.readlines()]
        print()
        print("-- Checking for invalid ChEBI IDs in grounding map --")
        for text, db_refs in gm.items():
            if db_refs is not None:
                for db_key, db_id in db_refs.items():
                    if db_key == 'CHEBI':
                        if db_id not in chebi_ids:
                            print("ERROR: ID %s in grounding map is "
                                  "not a valid CHEBI ID." % db_id)
    except IOError:
        pass

    # This check requires the requests package to be installed
    try:
        import requests
        import logging
        logging.getLogger('requests').setLevel(logging.CRITICAL)
        logging.getLogger('urllib3').setLevel(logging.CRITICAL)
        print()
        print("-- Checking for invalid PUBCHEM CIDs in grounding map --")
        pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/' + \
                      'cid/%s/description/XML'
        for text, db_refs in gm.items():
            if db_refs is not None:
                for db_key, db_id in db_refs.items():
                    if db_key == 'PUBCHEM':
                        res = requests.get(pubchem_url % db_id)
                        if res.status_code != 200:
                            print("ERROR: ID %s in grounding map is "
                                  "not a valid PUBCHEM ID." % db_id)
    except ImportError:
        pass

    print()
    print("-- Checking for Bioentities whose relationships are undefined  --")
    # Check the relationships for consistency with entities
    rel_missing_entities = []
    for ent in entities:
        found = False
        for subj, rel, obj in relationships:
            subj_ns = subj[0]
            subj_id = subj[1]
            obj_ns = obj[0]
            obj_id = obj[1]
            if subj_ns == 'BE' and subj_id == ent:
                found = True
                break
            if obj_ns == 'BE' and obj_id == ent:
                found = True
                break
        if not found:
            rel_missing_entities.append(ent)
            print("ERROR: ID %s has no known relations." % ent)

