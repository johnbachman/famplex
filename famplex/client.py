from collections import defaultdict, deque

from famplex.util import load_csv, construct_grounding_map
from famplex.locations import ENTITIES_PATH, EQUIVALENCES_PATH, \
    GROUNDING_MAP_PATH, RELATIONS_PATH, GENE_PREFIXES_PATH, DESCRIPTIONS_PATH


class FamplexGraph(object):
    """Provides methods for working with graph of FamPlex entities and relations

    FamPlex is an ontology of protein families and complexes. Individual terms
    are genes/proteins. There are  higher level terms for families and
    complexes. Terms can be connected by isa or partof relationships.
    X isa Y expressing that X is a member of the Y family; Z partof W
    expressing that Z is a constituent of the W complex. The FamPlex graph is
    a Directed Acyclic Graph (DAG) with top level families and complexes as
    sinks and individual proteins/genes as sources.

    If X isa Y or X partof Y we say that X is a child of Y and Y is a parent of
    X. This is contrary to the terminology used within graph theory for trees
    where X is a child of Y if there is an edge from Y to X. However, this is
    consistent with how these terms are often used in the hierarchical
    relationships of ontologies.

    We say Y is above X in the FamPlex ontology if there is a path of isa and
    partof edges from X to Y. We also say that Y is an ancestor of X.
    X is then below Y in the FamPlex ontology and we also say X is a descendant
    of Y.

    Attributes
    ----------
    root_classes : set
        Set of top level families and complexes in the FamPlex ontology
    """
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
        """Returns terms immediately above a given term in the FamPlex ontology

        Parameters
        ----------
        namespace : str
            Namespace for a term. For the FamPlex ontology this should be one
            of HGNC, FPLX, or UP
        id_ : str
            Unique identifier for term in namespace. HGNC unique ID, FPLX ID
            or Uniprot ID.
        relation_types : Optional[list]
            Set of relation types that input term can have with returned
            parent terms. The valid relation types are 'isa' and 'partof'.
            If argument is None then there are no restrictions on relation
            type.
            Default: None

        Returns
        -------
        list
            List of tuples of the form (namespace, id) specifying parent terms
            of the input term.
        """
        if relation_types is None:
            relation_types = ['isa', 'partof']
        return [(ns2, id2) for ns2, id2, rel in self._graph[(namespace, id_)]
                if rel in relation_types]

    def child_terms(self, namespace, id_, relation_types=None):
        """Returns terms immediately below a given term in the FamPlex ontology

        Parameters
        ----------
        namespace : str
            Namespace for a term. For the FamPlex ontology this should be one
            of HGNC, FPLX, or UP
        id_ : str
            Unique identifier for term in namespace. HGNC unique ID, FPLX ID
            or Uniprot ID.
        relation_types : Optional[list]
            Restrict edges to relation types in this list. The valid relation
            types are the strings 'isa' and 'partof'.
            If argument is None then both isa and partof relations are
            included. Default: None

        Returns
        -------
        list
            List of tuples of the form (namespace, id) specifying child terms
            of the input term..
        """
        if relation_types is None:
            relation_types = ['isa', 'partof']
        return [(ns2, id2) for ns2, id2, rel in
                self._reverse_graph[(namespace, id_)]
                if rel in relation_types]

    def root_terms(self, namespace, id_):
        """Returns top level terms associated with input term

        Parameters
        ----------
        namespace : str
            Namespace for a term. For the FamPlex ontology this should be one
            of HGNC, FPLX, or UP
        id_ : str
            Unique identifier for term in namespace. HGNC unique ID, FPLX ID
            or Uniprot ID.

        Returns
        -------
        list
            List of terms above the input that are top level families and/or
            complexes within the FamPlex ontology.
        """
        return [node for node in self._root_class_mapping[(namespace, id_)]]

    def ancestral_terms(self, namespace, id_, relation_types=None):
        """
        Return list of all prior terms in the FamPlex Ontology

        Parameters
        ----------
        namespace : str
            Namespace for a term. For the FamPlex ontology this should be one
            of HGNC, FPLX, or UP
        id_ : str
            Unique identifier for term in namespace. HGNC unique ID, FPLX ID
            or Uniprot ID.
        relation_types : Optional[list]
            Restrict edges to relation types in this list. The valid relation
            types are the strings 'isa' and 'partof'.
            If argument is None then both isa and partof relations are
            included. Default: None

        Returns
        -------
        list
           List of terms are returned in breadth first order following
           relations upward from bottom to top in the ontology.
        """
        if relation_types is None:
            relation_types = ['isa', 'partof']
        output = []
        for ns2, id2 in self._traverse(self._graph, (namespace, id_),
                                       relation_types):
            output.append((ns2, id2))
        return output[1:]

    def descendant_terms(self, namespace, id_, relation_types=None):
        """
        Return list of all following terms in the FamPlex Ontology

        Parameters
        ----------
        namespace : str
            Namespace for a term. For the FamPlex ontology this should be one
            of HGNC, FPLX, or UP
        id_ : str
            Unique identifier for term in namespace. HGNC unique ID, FPLX ID
            or Uniprot ID.
        relation_types : Optional[list]
            Restrict edges to relation types in this list. The valid relation
            types are the strings 'isa' and 'partof'.
            If argument is None then both isa and partof relations are
            included. Default: None

        Returns
        -------
        list
           List of terms are returned in breadth first order following
           relations backwards from top to bottom in the ontology.
        """
        if relation_types is None:
            relation_types = ['isa', 'partof']
        output = []
        for ns2, id2 in self._traverse(self._reverse_graph, (namespace, id_),
                                       relation_types):
            output.append((ns2, id2))
        return output[1:]

    def individual_members(self, namespace, id_, relation_types=None):
        """Return terms beneath a given term that are not families or complexes

        Parameters
        ----------
        namespace : str
            Namespace for a term. For the FamPlex ontology this should be one
            of HGNC, FPLX, or UP
        id_ : str
            Unique identifier for term in namespace. HGNC unique ID, FPLX ID
            or Uniprot ID.
        relation_types : list
            Restrict edges to relation types in this list. The valid relation
            types are the strings 'isa' and 'partof'.
            If argument is None then both isa and partof relations are
            included. Default: None

        Returns
        -------
        list
            List of terms beneath the input term that have no children
            themselves. If relation_types includes only 'isa', then these will
            be the individual genes within a family. If only 'partof' relations
            are included then these will be the individual members of a
            complex.  There are some terms that are families of
            complexes. These have both partof and isa relationships. In these
            cases the returned list can contain families or complexes
            if partof or isa relationships are excluded respectively.
        """
        if relation_types is None:
            relation_types = ['isa', 'partof']
        output = []
        for ns2, id2 in self.descendant_terms(namespace, id_,
                                              relation_types):
            if not (ns2, id2) in self._reverse_graph:
                output.append[(ns2, id2)]
        return output

    def isa(self, namespace1, id1, namespace2, id2):
        """Return true if one term has an isa relationship with another

        Parameters
        ----------
        namespace1 : str
            Namespace of first term
        id1 : str
            Identifier of first term
        namespace2 : str
            Namespace of second term
        id2 : str
            Identifier of second term
        
        Returns
        -------
        bool
            True if the term given by (namespace1, id1) has an isa relationship
            with the term given by (namespace2, id2).
        """
        return self._rel(namespace1, id1, namespace2, id2, ['isa'])

    def partof(self, namespace1, id1, namespace2, id2):
        """Return true if one term has a partof relationship with another

        Parameters
        ----------
        namespace1 : str
            Namespace of first term
        id1 : str
            Identifier of first term
        namespace2 : str
            Namespace of second term
        id2 : str
            Identifier of second term
        
        Returns
        -------
        bool
            True if the term given by (namespace1, id1) has a partof
            relationship with the term given by (namespace2, id2).
        """
        return self._rel(namespace1, id1, namespace2, id2, ['partof'])

    def refinement_of(self, namespace, id1, namespace2, id2):
        """Return true if one term either isa or partof holds

        Parameters
        ----------
        namespace1 : str
            Namespace of first term
        id1 : str
            Identifier of first term
        namespace2 : str
            Namespace of second term
        id2 : str
            Identifier of second term
        
        Returns
        -------
        bool
            True if the term given by (namespace1, id1) has either an isa
            or partof relationship with the term given by (namespace2, id2).
        """
        return self._rel(namespace, id1, namespace2, id2, ['isa', 'partof'])

    def dict_representation(self, namespace, id_):
        """Return a nested dictionary representation of a FamPlex term

        Parameters
        ----------
        namespace : str
        id_ : str

        Returns
        -------
        dict
            Nested dictionary representing structure of a FamPlex term.
            Keys are tuples with namespace, id pairs. Values are lists of
            tuples of nested dictionary representations and relationships,
            as in the example below.

            {('FPLX': 'ESR'): [({('HGNC', 'ESR1'): []}, 'isa'),
                               ({('HGNC', 'ESR2'): []}, 'isa')]}
        """
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

