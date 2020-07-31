"""Provides utilities for working with FampPlex entities and relations

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
X. We say Y is above X in the FamPlex ontology if there is a path of isa and
partof edges from X to Y. We also say that Y is an ancestor of X.
X is then below Y in the FamPlex ontology and we also say X is a descendant
of Y.
"""
import warnings
from typing import Container, Dict, List, Optional, Tuple

from famplex.graph import FamplexGraph

__all__ = ['in_famplex', 'parent_terms', 'child_terms', 'root_terms',
           'ancestral_terms', 'descendant_terms', 'individual_members', 'isa',
           'partof', 'refinement_of', 'dict_representation', 'equivalences',
           'reverse_equivalences', 'all_root_terms']


try:
    _famplex_graph = FamplexGraph()
except FileNotFoundError:
    warnings.warn(
        "Resource files are unavailable. If you've cloned this repository, "
        "run the script \"update_resources.py\" at the top level to move the "
        "resources into the package. See the README for more info.",
        Warning)


def in_famplex(namespace: str, id_: str) -> bool:
    """Returns True if input term is a member of the FamPlex ontology.

    Parameters
    ----------
    namespace : str
        Namespace for a term. This should be one of 'HGNC', 'FPLX' for
        FamPlex, or 'UP' for Uniprot.
    id_ : str
        Identifier for a term within namespace.

    Returns
    -------
    bool
    """
    return _famplex_graph.in_famplex(namespace, id_)


def parent_terms(namespace: str, id_: str,
                 relation_types: Optional[Container[str]] = None) \
                 -> List[Tuple[str, str]]:
    """Returns terms immediately above a given term in the FamPlex ontology

    Parameters
    ----------
    namespace : str
        Namespace for a term. This should be one of 'HGNC', 'FPLX' for
        FamPlex, or 'UP' for Uniprot.
    id_ : str
        Identifier for a term within namespace.
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
        of the input term. Values are sorted in case insensitive
        alphabetical order, first by namespace and then by id.

    Raises
    ------
    ValueError
        If (namespace, id_) does not correspond to a term in FamPlex.
    """
    if relation_types is None:
        relation_types = ['isa', 'partof']
    edges = _famplex_graph.parent_edges(namespace, id_)
    return [(ns2, id2) for ns2, id2, rel in edges if rel in relation_types]


def child_terms(namespace: str, id_: str,
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
        Identifier for a term within namespace.
    relation_types : Optional[list]
        Restrict edges to relation types in this list. The valid relation
        types are the strings 'isa' and 'partof'.
        If argument is None then both isa and partof relations are
        included. Default: None

    Returns
    -------
    list
        List of tuples of the form (namespace, id) specifying child terms
        of the input term. Values are sorted in case insensitive
        alphabetical order, first by namespace and then by id.

    Raises
    ------
    ValueError
        If (namespace, id_) does not correspond to a term in FamPlex.
    """
    if relation_types is None:
        relation_types = ['isa', 'partof']
    edges = _famplex_graph.child_edges(namespace, id_)
    return [(ns2, id2) for ns2, id2, rel in edges if rel in relation_types]


def root_terms(namespace: str, id_: str) -> List[Tuple[str, str]]:
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
    return _famplex_graph.root_terms(namespace, id_)


def ancestral_terms(namespace: str, id_: str,
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
       Edges from the same node are traversed in case insensitive
       alphabetical order, sorted first by namespace and then by id
       of the target node.

    Raises
    ------
    ValueError
        If (namespace, id_) does not correspond to a term in FamPlex.
    """
    _famplex_graph.raise_value_error_if_not_in_famplex(namespace, id_)
    if relation_types is None:
        relation_types = ['isa', 'partof']
    output = []
    for ns2, id2 in _famplex_graph.traverse((namespace, id_),
                                            relation_types, 'up'):
        output.append((ns2, id2))
    return output[1:]


def descendant_terms(namespace: str, id_: str,
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
       Edges from the same node are traversed in case insensitive
       alphabetical order, sorted first by namespace and then by id
       of the target node.

    Raises
    ------
    ValueError
        If (namespace, id_) does not correspond to a term in FamPlex.
    """
    _famplex_graph.raise_value_error_if_not_in_famplex(namespace, id_)
    if relation_types is None:
        relation_types = ['isa', 'partof']
    output = []
    for ns2, id2 in _famplex_graph.traverse((namespace, id_),
                                            relation_types, 'down'):
        output.append((ns2, id2))
    return output[1:]


def individual_members(namespace: str, id_: str,
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
        Values are sorted in case insensitive alphabetical order, first by
        namespace and then by id.

    Raises
    ------
    ValueError
        If (namespace, id_) does not correspond to a term in FamPlex.    Raises
    """
    if relation_types is None:
        relation_types = ['isa', 'partof']
    output = []
    for ns2, id2 in descendant_terms(namespace, id_, relation_types):
        if not child_terms(ns2, id2, relation_types=relation_types):
            output.append((ns2, id2))
    return sorted(output, key=lambda x: (x[0].lower(), x[1].lower()))


def isa(namespace1: str, id1: str, namespace2: str, id2: str) -> bool:
    """Return true if one term has an isa relationship with another

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
    return _famplex_graph.relation(namespace1, id1, namespace2, id2, ['isa'])


def partof(namespace1: str, id1: str, namespace2: str, id2: str) -> bool:
    """Return true if one term has a partof relationship with another

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
    return _famplex_graph.relation(namespace1, id1,
                                   namespace2, id2, ['partof'])


def refinement_of(namespace: str, id1: str, namespace2: str, id2: str) -> bool:
    """Return true if one term either isa or partof holds

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
    return _famplex_graph.relation(namespace, id1,
                                   namespace2, id2, ['isa', 'partof'])


def dict_representation(namespace: str,
                        id_: str) -> Dict[Tuple[str, str],
                                          List[Tuple[dict, str]]]:
    """Return a nested dictionary representation of a FamPlex term

    Parameters
    ----------
    namespace : str
        Namespace for a term. This should be one of 'HGNC', 'FPLX' for
        FamPlex, or 'UP' for Uniprot.
    id_ : str
        Identifier for a term within namespace.

    Returns
    -------
    dict
        Nested dictionary representing structure of a FamPlex term.
        Keys are tuples with namespace, id pairs. Values are lists of
        tuples of nested dictionary representations and relationships,
        as in the example below. Edges are sorted in case insensitive
        alphabetical order, first by namespace and then by id of the
        target node.

        {('FPLX', 'ESR'): [({('HGNC', 'ESR1'): []}, 'isa'),
                           ({('HGNC', 'ESR2'): []}, 'isa')]}

    Raises
    ------
    ValueError
        If (namespace, id_) does not correspond to a term in FamPlex.
    """
    out: Dict[Tuple[str, str], List[Tuple[dict, str]]] = \
        {(namespace, id_): []}
    edges = _famplex_graph.child_edges(namespace, id_)
    if not edges:
        return out
    for namespace2, id2, relation in edges:
        out[(namespace, id_)].\
            append((dict_representation(namespace2, id2), relation))
    return out


def equivalences(fplx_id: str,
                 namespaces: Optional[Container[str]] = None) -> \
                 List[Tuple[str, str]]:
    """Return list of equivalent terms from other namespaces.

    Parameters
    ----------
    fplx_id : str
        A valid Famplex ID

    namespaces : Optional[container]
        List of namespaces returned equivalences to which returned
        equivalences will be restricted. Can be used if one is interested
        in a particular type of equivalences.

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
    equivs = _famplex_graph.equivalences(fplx_id)
    if namespaces is not None:
        equivs = [(namespace, id_) for namespace, id_ in equivs
                  if namespace in namespaces]
    return equivs


def reverse_equivalences(namespace: str, id_: str) -> List[str]:
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
    return _famplex_graph.reverse_equivalences(namespace, id_)


def all_root_terms() -> List[Tuple[str, str]]:
    """Returns all top level families and complexes in FamPlex

    Returns
    -------
    list
        List of tuples of the form ('FPLX', id) where id runs over all
        top level families and complexes in FamPlex. List is in alphabetical
        order by id.
    """
    return _famplex_graph.root_classes
