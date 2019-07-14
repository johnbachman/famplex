import os
import csv
import sys
import datetime
import collections


if sys.version_info.major < 3:
    raise Exception('This script should be run in Python 3.')


Reference = collections.namedtuple('Reference', ['ns', 'id'])
Synonym = collections.namedtuple('Synonym', ['name', 'status'])


class OboTerm(object):
    def __init__(self, term_id, name, rels, synonyms=None, xrefs=None):
        self.term_id = term_id
        self.name = name
        self.synonyms = synonyms
        if xrefs is not None:
            self.xrefs = xrefs
        else:
            self.xrefs = []
        self.rels = rels

    def to_obo(self):
        obo_str = '[Term]\n'
        obo_str += 'id: %s:%s\n' % (self.term_id.ns, self.term_id.id)
        obo_str += 'name: %s\n' % self.name
        for synonym in self.synonyms:
            obo_str += 'synonym: "%s" %s []\n' % (synonym.name, synonym.status)
        for xref in self.xrefs:
            if xref.ns == 'BEL':
                entry = 'BEL:"%s"' % xref.id
            elif xref.ns == 'NXP':
                entry = 'NEXTPROT-FAMILY:%s' % xref.id[3:]
            elif xref.ns == 'PF':
                entry = 'XFAM:%s' % xref.id
            elif xref.ns == 'GO':
                entry = xref.id
            else:
                entry = '%s:%s' % (xref.ns, xref.id)
            obo_str += 'xref: %s\n' % entry
        for rel_type, rel_entries in self.rels.items():
            for ref in rel_entries:
                obo_str += '%s: %s:%s\n' % (rel_type, ref.ns, ref.id)
        return obo_str

    def __str__(self):
        return self.to_obo()


def get_obo_terms():
    obo_terms = []
    path_this = os.path.dirname(os.path.abspath(__file__))
    entities_file = os.path.join(path_this, os.pardir, 'entities.csv')
    grounding_file = os.path.join(path_this, os.pardir, 'grounding_map.csv')
    equiv_file = os.path.join(path_this, os.pardir, 'equivalences.csv')
    rel_file = os.path.join(path_this, os.pardir, 'relations.csv')
    with open(entities_file, 'r') as fh:
        entities = [l.strip() for l in fh.readlines()]
    with open(equiv_file, 'r') as fh:
        csvreader = csv.reader(fh, delimiter=str(u','), lineterminator='\r\n',
                               quoting=csv.QUOTE_MINIMAL,
                               quotechar=str(u'"'))
        equivalences = collections.defaultdict(list)
        for row in csvreader:
            source_ns, source_id, fplx_id = row
            equivalences[fplx_id].append((source_ns, source_id))
    with open(grounding_file) as fh:
        csvreader = csv.reader(fh, delimiter=str(u','), lineterminator='\r\n',
                               quoting=csv.QUOTE_MINIMAL,
                               quotechar=str(u'"'))
        textrefs = collections.defaultdict(list)
        for row in csvreader:
            text_str = row[0]
            namespaces = row[1::2]
            ids = row[2::2]
            for ns, id in zip(namespaces, ids):
                if ns == 'FPLX':
                    textrefs[id].append(text_str)
    with open(rel_file) as fh:
        rels = {entity: collections.OrderedDict(is_a=[],
                                                part_of=[],
                                                inverse_is_a=[],
                                                has_part=[])
                for entity in entities}
        csvreader = csv.reader(fh, delimiter=str(u','), lineterminator='\r\n',
                               quoting=csv.QUOTE_MINIMAL,
                               quotechar=str(u'"'))
        for row in csvreader:
            ns1, id1, rel, ns2, id2 = row
            if ns1 == 'FPLX':
                if rel == 'isa':
                    rels[id1]['is_a'].append(Reference(ns2, id2))
                elif rel == 'partof':
                    rels[id1]['part_of'].append(Reference(ns2, id2))
            if ns2 == 'FPLX':
                if rel == 'isa':
                    rels[id2]['inverse_is_a'].append(Reference(ns1, id1))
                elif rel == 'partof':
                    rels[id2]['has_part'].append(Reference(ns1, id1))

    # For each entity in famplex
    for entity in entities:
        entity_id = Reference('FPLX', entity)
        # Construct string name
        name = entity.replace('_', '-')
        # Get synonyms
        refs = textrefs.get(entity)
        synonyms = []
        if refs:
            for textref in refs:
                synonyms.append(Synonym(textref, 'EXACT'))
        # Get xrefs
        equivs = equivalences.get(entity)
        xrefs = []
        if equivs:
            for equiv in equivs:
                xrefs.append(Reference(equiv[0], equiv[1]))
        # If the entity has no isa relations, connect it to the root
        if not rels[entity]['is_a'] and not rels[entity]['part_of']:
            rels[entity]['is_a'].append(Reference('FPLX', 'root'))
        term = OboTerm(entity_id, name, rels[entity], synonyms, xrefs)
        obo_terms.append(term)
    obo_terms.append(OboTerm(Reference('FPLX', 'root'),
                             'PROTEIN-FAMILY-OR-COMPLEX',
                             {}, [], {}))
    return obo_terms


def save_obo_terms(obo_terms, output_file=None):
    date = datetime.datetime.today()
    date_str = date.strftime('%d:%m:%Y %H:%M')
    path_this = os.path.dirname(os.path.abspath(__file__))
    if not output_file:
        output_file = os.path.join(path_this, os.pardir, 'famplex.obo')
    with open(output_file, 'wt') as fh:
        fh.write('format-version: 1.2\n')
        fh.write('date: %s\n' % date_str)
        fh.write('\n')
        for term in obo_terms:
            obo_str = term.to_obo()
            fh.write(obo_str)
            fh.write('\n')


if __name__ == '__main__':
    obo_terms = get_obo_terms()
    path_this = os.path.dirname(os.path.abspath(__file__))
    save_obo_terms(obo_terms, os.path.join(path_this, 'famplex.obo'))

