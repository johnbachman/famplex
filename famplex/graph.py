"""Work with the graph of FamPlex entities and relations."""
from typing import Container, Dict, Generator, List, Tuple

from collections import defaultdict, deque

from famplex.load import load_entities, load_equivalences, load_relations


class FamplexGraph(object):
    """Provides methods for working with graph of FamPlex entities and relations

    FamPlex is an ontology of protein families and complexes. Individual terms
    are genes/proteins. There are  higher level terms for families and
    complexes. Terms can be connected by isa or partof relationships.
    X isa Y expressing that X is a member of the Y family; Z partof W
    expressing that Z is a constituent of the W complex.

    Each term in the FamPlex ontology exists within a namespace and has an
    identifier which is unique within that namespace. Individual genes and
    proteins have either HGNC or Uniprot as a namespace. FamPlex has its own
    namespace for families and complexes and the unique identifiers are
    designed to be human readable. Identifiers for Uniprot are simply Uniprot
    IDs. For HGNC the HGNC Symbol is used instead of the HGNC unique ID.

    If X isa Y or X partof Y we say that X is a child of Y and Y is a parent of
    X. This is contrary to the terminology used within graph theory for trees
    where X is a child of Y if there is an edge from Y to X. However, this is
    consistent with how these terms are often used for the hierarchical
    relationships within ontologies.

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
        # Graphs are stored internally as a dictionary mapping tuples of
        # the form (namespace, id) to a list of tuples of the form
        # (namespace, id, relation_type). This is a variant of the adjacency
        # list representation of a graph but allowing for multiple edge types.

        # Contains forward isa and partof relationships between terms
        graph = defaultdict(list)
        # Contains reversed isa and partof relationships
        reverse_graph = defaultdict(list)
        relations = load_relations()
        left_set = set()
        right_set = set()
        # Loop through table populating edges of the above graphs.
        # By looking at the set of all terms that appear on the right in the
        # relations.csv table which do not appear on the left we can identify
        # the top level families and complexes within famplex.
        for namespace1, id1, relation, namespace2, id2 in relations:
            graph[(namespace1, id1)].append((namespace2, id2, relation))
            reverse_graph[(namespace2, id2)].\
                append((namespace1, id1, relation))
            left_set.add((namespace1, id1))
            right_set.add((namespace2, id2))
        graph = dict(graph)
        reverse_graph = dict(reverse_graph)
        # Sort edges in adjaceny lists by alphabetical order
        for node, edges in graph.items():
            graph[node] = sorted(edges, key=lambda x: (x[0].lower(),
                                                       x[1].lower()))
        for node, edges in reverse_graph.items():
            reverse_graph[node] = sorted(edges,
                                         key=lambda x: (x[0].lower(),
                                                        x[1].lower()))
        self._graph: Dict[Tuple[str, str], List[Tuple[str, str, str]]] = graph

        self._reverse_graph: Dict[Tuple[str, str],
                                  List[Tuple[str, str, str]]] = reverse_graph
        root_class_mapping = defaultdict(list)
        root_classes = sorted(right_set - left_set, key=lambda x: x[1].lower())
        # Build up an dictionary mapping terms to the top level families
        # or complexes to which they belong. Families and complexes can overlap
        # so there can be multiple top level terms above a given term.
        for entry in root_classes:
            for node in self.traverse(entry, ['isa', 'partof'],
                                      direction='down'):
                root_class_mapping[node].append(entry)
        root_class_mapping = dict(root_class_mapping)
        entities = load_entities()
        for entity in entities:
            entry = ('FPLX', entity)
            if entry not in root_class_mapping:
                root_class_mapping[entry] = [entry]
        for node, roots in root_class_mapping.items():
            root_class_mapping[node] = sorted(roots,
                                              key=lambda x: (x[0].lower(),
                                                             x[1].lower()))

        equivalences = defaultdict(list)
        reverse_equivalences = defaultdict(list)
        for ns, id_, fplx_id in load_equivalences():
            equivalences[fplx_id].append((ns, id_))
            reverse_equivalences[(ns, id_)].append(fplx_id)
        equivalences = dict(equivalences)
        reverse_equivalences = dict(reverse_equivalences)
        # Blank lines are to aid in reading of type hints
        self.root_classes: List[Tuple[str, str]] = root_classes

        self._root_class_mapping: Dict[Tuple[str, str],
                                       List[Tuple[str, str]]] = \
            root_class_mapping

        self._equivalences: Dict[str, List[Tuple[str, str]]] = equivalences
        self._reverse_equivalences: Dict[Tuple[str, str], List[str]] = \
            reverse_equivalences
        self.__error_message = 'Given input is not in the FamPlex ontology.'

    def in_famplex(self, namespace: str, id_: str) -> bool:
        """Returns True if input term is a member of the FamPlex ontology.

        Parameters
        ----------
        namespace : str
            Namespace for a term. This should be one of 'HGNC', 'FPLX' for
            FamPlex, or 'UP' for Uniprot.
        id_ : str
            Identifier for a term within namespace. See the FamplexGraph
            class Docstring for more info.

        Returns
        -------
        bool
        """
        return (namespace, id_) in self._root_class_mapping

    def raise_value_error_if_not_in_famplex(self, namespace: str,
                                            id_: str) -> None:
        """Raise a value error if input is not in FamPlex ontology

        This can be used in functions where we desire an exception to be
        raised when the input is not in the FamPlex ontology but the 
        most natural way of writing the function will not lead to any
        exceptions being raised.

        Parameters
        ----------
        namespace : str
            Namespace for a term. This should be one of 'HGNC', 'FPLX' for
            FamPlex, or 'UP' for Uniprot.
        id_ : str
            Identifier for a term within namespace. See the FamplexGraph
            class Docstring for more info.

        Raises
        ------
        ValueError
        If (namespace, id_) does not correspond to a term in FamPlex.
        """
        if not self.in_famplex(namespace, id_):
            raise ValueError(self.__error_message)

    def parent_edges(self, namespace: str,
                     id_: str) -> List[Tuple[str, str, str]]:
        """Returns node and relation type for all parents of input

        Parameters
        ----------
        namespace : str
            Namespace for a term. This should be one of 'HGNC', 'FPLX' for
            FamPlex, or 'UP' for Uniprot.
        id_ : str
            Identifier for a term within namespace. See the FamplexGraph
            class Docstring for more info.

        Returns
        -------
        list
            List of all tuples of the form (namespace, id, relation_type) where
            (namespace, id) is a parent of the input and relation_type is the
            type of relation connecting them.

        Raises
        ------
        ValueError
        If (namespace, id_) does not correspond to a term in FamPlex.
        """
        edges = self._graph.get((namespace, id_))
        if edges is None:
            self.raise_value_error_if_not_in_famplex(namespace, id_)
            return []
        return edges

    def child_edges(self, namespace: str,
                    id_: str) -> List[Tuple[str, str, str]]:
        """Returns node and relation type for all children of input

        Parameters
        ----------
        namespace : str
            Namespace for a term. This should be one of 'HGNC', 'FPLX' for
            FamPlex, or 'UP' for Uniprot.
        id_ : str
            Identifier for a term within namespace. See the FamplexGraph
            class Docstring for more info.

        Returns
        -------
        list
            List of all tuples of the form (namespace, id, relation_type) where
            (namespace, id) is a child of the input and relation_type is the
            type of relation connecting them.

        Raises
        ------
        ValueError
        If (namespace, id_) does not correspond to a term in FamPlex.
        """
        edges = self._reverse_graph.get((namespace, id_))
        if edges is None:
            self.raise_value_error_if_not_in_famplex(namespace, id_)
            return []
        return edges

    def root_terms(self, namespace: str, id_: str) -> List[Tuple[str, str]]:
        """Returns top level terms above the input term

        Parameters
        ----------
        namespace : str
            Namespace for a term. This should be one of 'HGNC', 'FPLX' for
            FamPlex, or 'UP' for Uniprot.
        id_ : str
            Identifier for a term within namespace. See the FamplexGraph
            class Docstring for more info.

        Returns
        -------
        list
            List of terms above the input that are top level families and/or
            complexes within the FamPlex ontology. Values are sorted in case
            insensitive alphabetical order, first by namespace and then by id.

        Raises
        ------
        ValueError
            If (namespace, id_) does not correspond to a term in FamPlex.
        """
        roots = self._root_class_mapping.get((namespace, id_))
        if roots is None:
            raise ValueError(self.__error_message)
        return roots

    def equivalences(self, fplx_id: str) -> List[Tuple[str, str]]:
        """Return list of equivalent terms from other namespaces.

        Parameters
        ----------
        fplx_id : str
            A valid Famplex ID

        Returns
        -------
        list
            List of tuples of the form (namespace, id) of equivalent terms
            from other namespaces.

        Raises
        ------
        ValueError
            If fplx_id an ID in the FamPlex ontology.
        """
        self.raise_value_error_if_not_in_famplex('FPLX', fplx_id)
        equiv = self._equivalences.get(fplx_id)
        if equiv is None:
            return []
        return equiv

    def reverse_equivalences(self, namespace: str, id_: str) -> List[str]:
        """Get equivalent FamPlex terms to a given term from another namespace

        Parameters
        ----------
        namespace : str
            Namespace of a term
        id_ : str
            id_ of a term

        Returns
        -------
        list
            List of FamPlex IDs for families or complexes equivalent to the
            term given by (namespace, id_)
        """
        equiv = self._reverse_equivalences.get((namespace, id_))
        equiv = [] if equiv is None else equiv
        return equiv

    def relation(self, namespace1: str, id1: str,
                 namespace2: str, id2: str,
                 relation_types: Container[str]) -> bool:
        """General function for determining if two entities are related

        Parameters
        ----------
        namespace1 : str
            Namespace of first term. This should be one of 'HGNC', 'FPLX' for
            FamPlex, or 'UP' for Uniprot.
        id1 : str
            Identifier of first term.
        namespace2 : str
            Namespace of second term. This should be one of 'HGNC', 'FPLX' for
            FamPlex, or 'UP' for Uniprot.
        id2 : str
            Identifier of second term.
        relation_types : container
            Function returns True if the first term is connected to the second
            by one of the relations in this container. Valid relations are
            'isa', and 'partof'.

        Returns
        -------
        bool
            True if the term given by (namespace1, id1) has one of the
            specified relations with the term given by (namespace2, id2). Will
            return False if either of (namespace1, id1) or (namespace2, id2) is
            not in the FamPlex ontology.
        """
        roots1 = self._root_class_mapping.get((namespace1, id1))
        roots2 = self._root_class_mapping.get((namespace2, id2))
        if roots1 is None or roots2 is None:
            return False
        if set(roots1) & set(roots2):
            node1, node2 = (namespace1, id1), (namespace2, id2)
            for node in self.traverse(node1, relation_types,
                                      direction='up'):
                if node2 == node:
                    return True
        return False

    def traverse(self, source: Tuple[str, str],
                 relation_types: Container[str],
                 direction: str) -> Generator[Tuple[str, str], None, None]:
        """Function for traversing FampPlex graph in breadth first order

        Parameters
        ----------
        source : tuple
            Tuple of the form (namespace, id) specifying where traversal is to
            begin.

        relation_types : container
            Traversal will follow edges from these specified relation_types.
            Valid relation types are isa and partof.

        direction : str
            One of 'up' or 'down'. If 'up' traversal will follow isa and partof
            edges in breadth first order to nodes above the source. If 'down'
            traversal will follow reversed edges to nodes below the source.

        Returns
        -------
        generator
            Generator iterating through nodes in the traversal. The source node
            is included in the traversal.
        """
        if direction == 'down':
            graph = self._reverse_graph
        elif direction == 'up':
            graph = self._graph
        else:
            raise ValueError
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
