from __future__ import print_function, unicode_literals
import csv
import pygraphviz as pgv

rows = []

graph = pgv.AGraph(name='relations', directed=True, rankdir='LR')
hgnc_style = {'color': 'lightgray', 'style': 'filled',
               'fontname': 'arial'}
be_style = {'color': 'pink', 'style': 'filled',
               'fontname': 'arial'}
edge_style = {'fontname': 'arial'}

with open('relations.csv') as f:
    csvreader = csv.reader(f, delimiter=',', lineterminator='\r\n',
                           quoting=csv.QUOTE_MINIMAL,
                           quotechar='"')
    nodes = set([])
    isa_edges = set([])
    partof_edges = set([])
    for row in csvreader:
        ns1, id1, rel, ns2, id2 = row
        def node_label(ns, id):
            return '%s:%s' % (ns, id)
        for ns, id in ((ns1, id1), (ns2, id2)):
            if ns == 'HGNC':
                graph.add_node(node_label(ns, id), **hgnc_style)
            elif ns == 'BE':
                graph.add_node(node_label(ns, id), **be_style)
        edge = (node_label(ns1, id1), node_label(ns2, id2))
        if rel == 'isa':
            graph.add_edge(edge, label='isa', **edge_style)
        elif rel == 'partof':
            graph.add_edge(edge, label='partof', **edge_style)

graph.draw('relations.pdf', prog='dot')

