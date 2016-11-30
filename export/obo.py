from __future__ import print_function, unicode_literals
import os
import csv
import datetime
import collections

Reference = collections.namedtuple('Reference', ['ns', 'id'])
Synonym = collections.namedtuple('Synonym', ['name', 'status'])

class OboTerm(object):
    def __init__(self, term_id, name, synonyms=None, xrefs=None, isas=None):
        self.term_id = term_id
        self.name = name
        self.synonyms = synonyms
        if xrefs is not None:
            self.xrefs = xrefs
        else:
            self.xrefs = []
        if isas is not None:
            self.isas = isas
        else:
            self.isas = []

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
            else:
                entry = '%s:%s' % (xref.ns, xref.id)
            obo_str += 'xref: %s\n' % entry
        for isa in self.isas:
            obo_str += 'is_a: %s:%s\n' % (isa.ns, isa.id)
        return obo_str

    def __str__(self):
        return self.to_obo()

def get_obo_terms():
    obo_terms = []
    path_this = os.path.dirname(os.path.abspath(__file__))
    entities_file = os.path.join(path_this, os.pardir, 'entities.csv')
    grounding_file = os.path.join(path_this, os.pardir, 'grounding_map.csv')
    equiv_file = os.path.join(path_this, os.pardir, 'equivalences.csv')
    # For each entity in bioentities
    with open(entities_file, 'r') as fh:
        entities = [l.strip() for l in fh.readlines()]
    with open(equiv_file, 'r') as fh:
        csvreader = csv.reader(fh, delimiter=str(u','), lineterminator='\r\n',
                               quoting=csv.QUOTE_MINIMAL,
                               quotechar=str(u'"'))
        equivalences = collections.defaultdict(list)
        for row in csvreader:
            source_ns, source_id, be_id = row
            equivalences[be_id].append((source_ns, source_id))
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
                if ns == 'BE':
                    textrefs[id].append(text_str)
    for entity in entities:
        entity_id = Reference('BE', entity)
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
        # Get isa
        isas = []
        isas.append(Reference('ONT', 'PROTEIN-FAMILY'))
        term = OboTerm(entity_id, name, synonyms, xrefs, isas)
        obo_terms.append(term)
    return obo_terms

def save_obo_terms(obo_terms, output_file=None):
    date = datetime.datetime.today()
    date_str = date.strftime('%d:%m:%Y %H:%M')
    path_this = os.path.dirname(os.path.abspath(__file__))
    if not output_file:
        output_file = os.path.join(path_this, os.pardir, 'bioentities.obo')
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
    save_obo_terms(obo_terms, 'bioentities.obo')

