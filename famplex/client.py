from collections import defaultdict, deque

from famplex.util import load_csv, construct_grounding_map
from famplex.locations import ENTITIES_PATH, EQUIVALENCES_PATH, \
    GROUNDING_MAP_PATH, RELATIONS_PATH, GENE_PREFIXES_PATH, DESCRIPTIONS_PATH


class FamplexGraph(object):
    def __init__(self):
        graph = defaultdict(list)
        reverse_graph = defaultdict(list)
        relations = load_csv(RELATIONS_PATH)
        left_set, right_set = set(), set()
        for namespace1, id1, relation, namespace2, id2 in relations:
            graph[(namespace1, id1)].append((namespace2, id2, relation))
            reverse_graph[(namespace2, id2)].\
                append((namespace1, id1, relation))
            left_set.add((namespace1, id1))
            right_set.add((namespace2, id2))
        top_level_mapping = {}
        top_level = right_set - left_set
        for entry in top_level:
            for node in self._traverse(reverse_graph, entry,
                                       ['isa', 'partof']):
                top_level_mapping[node] = entry
        self._top_level = dict(top_level_mapping)
        self._graph = dict(graph)
        self._reverse_graph = dict(reverse_graph)

    def get_parents(self, namespace, id_):
        return [(ns, id2) for ns, id2, _ in self._graph[(namespace, id_)]]

    def get_children(self, namespace, id_):
        return [(ns, id2) for ns, id2, _ in
                self._reverse_graph[(namespace, id_)]]

    def get_progenitors(self, namespace, id_):
        return [node for node in self._top_level[(namespace, id_)]]

    def get_ancestors(self, namespace, id_):
        for ns, i in self._traverse(self._graph, (namespace, id_),
                                    ['isa', 'partof']):
            yield (ns, i)

    def get_descendants(self, namespace, id_):
        for ns, i in self._traverse(self._reverse_graph, (namespace, id_),
                                    ['isa', 'partof']):
            yield (ns, i)

    def get_leaves(self, namespace, id_):
        for ns, id_ in self.get_descendants(namespace, id_):
            if ns != 'FPLX':
                yield (ns, id_)

    def isa(self, namespace1, id1, namespace2, id2):
        return self._rel(namespace1, id1, namespace2, id2, ['isa'])

    def partof(self, namespace1, id1, namespace2, id2):
        return self._rel(namespace1, id1, namespace2, id2, ['partof'])

    def refinement_of(self, namespace, id1, namespace2, id2):
        return self._rel(namespace, id1, namespace2, id2, ['isa', 'partof'])

    def get_type(self, namespace, id_):
        edges = self._reverse_graph.get((namespace, id_))
        if edges is None:
            raise ValueError(f'{namespace}:{id_} is not in the'
                             ' FamPlex ontology')
        edge_types = set([edge[-1] for edge in edges])
        if not edge_types:
            output = 'gene/protein'
        elif 'partof' in edge_types:
            output = 'complex'
        elif edge_types == set(['isa']):
            output = 'family'
        else:
            output = None
        return output
    def _rel(self, namespace1, id1, namespace2, id2, relation_types):
        roots1 = self._top_level[(namespace1, id1)]
        roots2 = self._top_level[(namespace2, id2)]
        common_roots = roots1.keys() & roots2.keys()
        if not common_roots:
            return False
        root = common_roots.pop()
        node1, node2 = (namespace1, id1), (namespace2, id2)
        for node in self._traverse(self._graph, node1,
                                   relation_types):
            if node2 == node:
                return True
        return False

    def _traverse(self, graph, source, relation_types):
        visited = {source}
        queue = deque([source])
        while queue:
            node = queue.pop()
            try:
                children = graph[node]
            except KeyError:
                children = []
            for ns, id_, rel in children:
                if (ns, id_) not in visited and rel in relation_types:
                    queue.appendleft((ns, id_))
                    visited.add((ns, id_))
            yield node


def load_grounding_map():
    rows = load_csv(GROUNDING_MAP_PATH)
    return construct_grounding_map(rows)


def load_equivalences():
    return load_csv(EQUIVALENCES_PATH)


def load_entitites():
    return load_csv(ENTITIES_PATH)


def load_relations():
    return load_csv(RELATIONS_PATH)


def load_gene_prefixes():
    return load_csv(GENE_PREFIXES_PATH)


def load_descriptions():
    return load_csv(DESCRIPTIONS_PATH)

