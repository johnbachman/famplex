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
        root_class_mapping = {}
        root_classes = right_set - left_set
        for entry in root_classes:
            for node in self._traverse(reverse_graph, entry,
                                       ['isa', 'partof']):
                root_class_mapping[node] = entry
        self.root_classes = root_classes
        self._root_class_mapping = dict(root_class_mapping)
        self._graph = dict(graph)
        self._reverse_graph = dict(reverse_graph)

    def parent_terms(self, namespace, id_, relation_types=None):
        if relation_types is None:
            relation_types = ['isa', 'partof']
        return [(ns2, id2) for ns2, id2, rel in self._graph[(namespace, id_)]
                if rel in relation_types]

    def child_terms(self, namespace, id_, relation_types=None):
        if relation_types is None:
            relation_types = ['isa', 'partof']
        return [(ns2, id2) for ns2, id2, rel in
                self._reverse_graph[(namespace, id_)]
                if rel in relation_types]

    def root_terms(self, namespace, id_):
        return [node for node in self._root_class_mapping[(namespace, id_)]]

    def ancestral_terms(self, namespace, id_, relation_types=None):
        if relation_types is None:
            relation_types = ['isa', 'partof']
        output = []
        for ns2, id2 in self._traverse(self._graph, (namespace, id_),
                                       relation_types):
            output.append((ns2, id2))
        return output[1:]

    def descendant_terms(self, namespace, id_, relation_types=None):
        if relation_types is None:
            relation_types = ['isa', 'partof']
        output = []
        for ns2, id2 in self._traverse(self._reverse_graph, (namespace, id_),
                                       relation_types):
            output.append((ns2, id2))
        return output[1:]

    def individual_members(self, namespace, id_, relation_types=None):
        if relation_types is None:
            relation_types = ['isa', 'partof']
        output = []
        for ns2, id2 in self.descendant_terms(namespace, id_):
            if ns2 != 'FPLX':
                output.append[(ns2, id2)]
        return output

    def isa(self, namespace1, id1, namespace2, id2):
        return self._rel(namespace1, id1, namespace2, id2, ['isa'])

    def partof(self, namespace1, id1, namespace2, id2):
        return self._rel(namespace1, id1, namespace2, id2, ['partof'])

    def refinement_of(self, namespace, id1, namespace2, id2):
        return self._rel(namespace, id1, namespace2, id2, ['isa', 'partof'])

    def category(self, namespace, id_):
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

    def dict_representation(self, namespace, id_):
        out = {(namespace, id_): []}
        edges = self._reverse_graph.get((namespace, id_))
        if edges is None:
            raise ValueError(f'{namespace}:{id_} is not in the FamPlex'
                             ' ontology')
        if not edges:
            return out
        for namespace2, id2, relation in edges:
            out[(namespace, id_)].\
                append((self.dict_representation(namespace2, id2),
                        relation))
        return out

    def _rel(self, namespace1, id1, namespace2, id2, relation_types):
        roots1 = self._root_class_mapping.get((namespace1, id1))
        roots2 = self._root_class_mapping.get((namespace2, id2))
        if roots1 is None or roots2 is None:
            return False
        if roots1.keys() & roots2.keys():
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

