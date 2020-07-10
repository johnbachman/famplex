"""Work with the graph of FamPlex entities and relations."""
from typing import Container, Dict, Generator, List, Optional, Tuple

from collections import defaultdict, deque

from famplex.load import load_equivalences, load_relations


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
        root_class_mapping = defaultdict(list)
        root_classes = sorted(right_set - left_set, key=lambda x: x[1].lower())
        # Build up an dictionary mapping terms to the top level families
        # or complexes to which they belong. Families and complexes can overlap
        # so there can be multiple top level terms above a given term.
        for entry in root_classes:
            for node in self._traverse(reverse_graph, entry,
                                       ['isa', 'partof']):
                root_class_mapping[node].append(entry)
        root_class_mapping = dict(root_class_mapping)
        for node, roots in root_class_mapping.items():
            root_class_mapping[node] = sorted(roots,
                                              key=lambda x: (x[0].lower(),
                                                             x[1].lower()))
        equivalences = defaultdict(list)
        for ns, id_, fplx_id in load_equivalences():
            equivalences[fplx_id].append((ns, id_))
        equivalences = dict(equivalences)
        # Blank lines are to aid in reading of type hints
        self.root_classes: List[Tuple[str, str]] = root_classes

        self._root_class_mapping: Dict[Tuple[str, str],
                                       List[Tuple[str, str]]] = \
            root_class_mapping

        self._graph: Dict[Tuple[str, str], List[Tuple[str, str, str]]] = graph

        self._reverse_graph: Dict[Tuple[str, str],
                                  List[Tuple[str, str, str]]] = reverse_graph

        self._equivalences: Dict[str, List[Tuple[str, str]]] = equivalences

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

    def parent_terms(self, namespace: str,
                     id_: str,
                     relation_types:
                     Optional[Container[str]] = None) -> \
            List[Tuple[str, str]]:
        """Returns terms immediately above a given term in the FamPlex ontology

        Parameters
        ----------
        namespace : str
            Namespace for a term. This should be one of 'HGNC', 'FPLX' for
            FamPlex, or 'UP' for Uniprot.
        id_ : str
            Identifier for a term within namespace. See the FamplexGraph
            class Docstring for more info.
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

        Raises
        ------
        ValueError
            If (namespace, id_) does not correspond to a term in FamPlex.
        """
        if relation_types is None:
            relation_types = ['isa', 'partof']
        edges = self._graph.get((namespace, id_))
        if edges is None:
            if (namespace, id_) in self.root_classes:
                return []
            raise ValueError(self.__error_message)
        return [(ns2, id2) for ns2, id2, rel in edges if rel in relation_types]

    def child_terms(self, namespace: str, id_: str,
                    relation_types:
                    Optional[Container[str]] = None) -> \
            List[Tuple[str, str]]:
        """Returns terms immediately below a given term in the FamPlex ontology

        Parameters
        ----------
        namespace : str
            Namespace for a term. This should be one of 'HGNC', 'FPLX' for
            FamPlex, or 'UP' for Uniprot.
        id_ : str
            Identifier for a term within namespace. See the FamplexGraph
            class Docstring for more info.
        relation_types : Optional[list]
            Restrict edges to relation types in this list. The valid relation
            types are the strings 'isa' and 'partof'.
            If argument is None then both isa and partof relations are
            included. Default: None

        Returns
        -------
        list
            List of tuples of the form (namespace, id) specifying child terms
            of the input term.

        Raises
        ------
        ValueError
            If (namespace, id_) does not correspond to a term in FamPlex.
        """
        if relation_types is None:
            relation_types = ['isa', 'partof']
        edges = self._reverse_graph.get((namespace, id_))
        if edges is None:
            if (namespace, id_) in self._graph:
                return []
            raise ValueError(self.__error_message)
        return [(ns2, id2) for ns2, id2, rel in edges if rel in relation_types]

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
            complexes within the FamPlex ontology.

        Raises
        ------
        ValueError
            If (namespace, id_) does not correspond to a term in FamPlex.
        """
        roots = self._root_class_mapping.get((namespace, id_))
        if roots is None:
            raise ValueError(self.__error_message)
        return roots

    def ancestral_terms(self, namespace: str, id_: str,
                        relation_types:
                        Optional[Container[str]] = None) -> \
            List[Tuple[str, str]]:
        """
        Return list of all terms above a given term in the FamPlex Ontology

        Parameters
        ----------
        namespace : str
            Namespace for a term. This should be one of 'HGNC', 'FPLX' for
            FamPlex, or 'UP' for Uniprot.
        id_ : str
            Identifier for a term within namespace. See the FamplexGraph
            class Docstring for more info.
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
        if not self.in_famplex(namespace, id_):
            raise ValueError(self.__error_message)
        if relation_types is None:
            relation_types = ['isa', 'partof']
        output = []
        for ns2, id2 in self._traverse(self._graph, (namespace, id_),
                                       relation_types):
            output.append((ns2, id2))
        return output[1:]

    def descendant_terms(self, namespace: str, id_: str,
                         relation_types:
                         Optional[Container[str]] = None) -> \
            List[Tuple[str, str]]:
        """
        Return list of all terms below a given term in the FamPlex Ontology

        Parameters
        ----------
        namespace : str
            Namespace for a term. This should be one of 'HGNC', 'FPLX' for
            FamPlex, or 'UP' for Uniprot.
        id_ : str
            Identifier for a term within namespace. See the FamplexGraph class
            Docstring for more info.
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
        if not self.in_famplex(namespace, id_):
            raise ValueError(self.__error_message)
        if relation_types is None:
            relation_types = ['isa', 'partof']
        output = []
        for ns2, id2 in self._traverse(self._reverse_graph, (namespace, id_),
                                       relation_types):
            output.append((ns2, id2))
        return output[1:]

    def individual_members(self, namespace: str, id_: str,
                           relation_types:
                           Optional[Container[str]] = None) -> \
            List[Tuple[str, str]]:
        """Return terms beneath a given term that are not families or complexes

        Parameters
        ----------
        namespace : str
            Namespace for a term. This should be one of 'HGNC', 'FPLX' for
            FamPlex, or 'UP' for Uniprot.
        id_ : str
            Identifier for a term within namespace. See the Famplexgraph class
            Docstring for more info.
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
            if not self.child_terms(ns2, id2, relation_types=relation_types):
                output.append((ns2, id2))
        return sorted(output, key=lambda x: (x[0].lower(), x[1].lower()))

    def isa(self, namespace1: str, id1: str,
            namespace2: str, id2: str) -> bool:
        """Return true if one term has an isa relationship with another

        See the FamplexGraph class Docstring for more info on namespaces and
        IDs.

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

        Returns
        -------
        bool
            True if the term given by (namespace1, id1) has an isa relationship
            with the term given by (namespace2, id2). Will return False if
            either of (namespace1, id1) or (namespace2, id2) is not in the
            FamPlex ontology.
        """
        return self._rel(namespace1, id1, namespace2, id2, ['isa'])

    def partof(self, namespace1: str, id1: str,
               namespace2: str, id2: str) -> bool:
        """Return true if one term has a partof relationship with another

        See the FamplexGraph class Docstring for more info on namespaces and
        IDs.

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

        Returns
        -------
        bool
            True if the term given by (namespace1, id1) has a partof
            relationship with the term given by (namespace2, id2). Will return
            False if either of (namespace1, id1) or (namespace2, id2) is not in
            the FamPlex ontology.
        """
        return self._rel(namespace1, id1, namespace2, id2, ['partof'])

    def refinement_of(self, namespace: str, id1: str,
                      namespace2: str, id2: str) -> bool:
        """Return true if one term either isa or partof holds

        See the FamplexGraph class Docstring for more info on namepaces and
        IDs.

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

        Returns
        -------
        bool
            True if the term given by (namespace1, id1) has either an isa or
            partof relationship with the term given by (namespace2, id2). Will
            return False if either of (namespace1, id1) or (namespace2, id2) is
            not in the FamPlex ontology.
        """
        return self._rel(namespace, id1, namespace2, id2, ['isa', 'partof'])

    def dict_representation(self, namespace: str,
                            id_: str) -> Dict[Tuple[str, str],
                                              List[Tuple[dict, str]]]:
        """Return a nested dictionary representation of a FamPlex term

        Parameters
        ----------
        namespace : str
            Namespace for a term. This should be one of 'HGNC', 'FPLX' for
            FamPlex, or 'UP' for Uniprot.
        id_ : str
            Identifier for a term within namespace. See the Famplexgraph class
            Docstring for more info.

        Returns
        -------
        dict
            Nested dictionary representing structure of a FamPlex term.
            Keys are tuples with namespace, id pairs. Values are lists of
            tuples of nested dictionary representations and relationships,
            as in the example below.

            {('FPLX', 'ESR'): [({('HGNC', 'ESR1'): []}, 'isa'),
                               ({('HGNC', 'ESR2'): []}, 'isa')]}

        Raises
        ------
        ValueError
            If (namespace, id_) does not correspond to a term in FamPlex.
        """
        out: Dict[Tuple[str, str], List[Tuple[dict, str]]] = \
            {(namespace, id_): []}
        edges = self._reverse_graph.get((namespace, id_))
        if not edges:
            if (namespace, id_) in self._graph:
                return out
            raise ValueError(self.__error_message)
        for namespace2, id2, relation in edges:
            out[(namespace, id_)].\
                append((self.dict_representation(namespace2, id2),
                        relation))
        return out

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
        equiv = self._equivalences.get(fplx_id)
        if equiv is None:
            raise ValueError('Input ID does not exist in FamPlex.')
        return equiv

    def _rel(self, namespace1: str, id1: str,
             namespace2: str, id2: str,
             relation_types: Container[str]) -> bool:
        roots1 = self._root_class_mapping.get((namespace1, id1))
        roots2 = self._root_class_mapping.get((namespace2, id2))
        if roots1 is None or roots2 is None:
            return False
        if set(roots1) & set(roots2):
            node1, node2 = (namespace1, id1), (namespace2, id2)
            for node in self._traverse(self._graph, node1,
                                       relation_types):
                if node2 == node:
                    return True
        return False

    def _traverse(self, graph: Dict[Tuple[str, str],
                                    List[Tuple[str, str, str]]],
                  source: Tuple[str, str],
                  relation_types:
                  Container[str]) -> Generator[Tuple[str, str], None, None]:
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
