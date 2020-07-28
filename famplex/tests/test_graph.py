import pytest


@pytest.fixture
def famplex_graph():
    from famplex.graph import FamplexGraph
    return FamplexGraph()


@pytest.mark.parametrize('test_input,expected',
                         [(('FPLX', 'AMPK'), True),
                          (('FPLX', 'AKT'), True),
                          (('HGNC', 'AKT1'), True),
                          (('FPLX',
                            'If_this_is_ever_a_real_family_then_we'
                            '_need_to_talk_to_the_biologists_who_named_it'),
                           False),
                          (('HGNC', 'GENE'), False)])
def test_in_famplex(famplex_graph, test_input, expected):
    assert famplex_graph.in_famplex(*test_input) == expected


@pytest.mark.parametrize('test_input,expected',
                         [(('HGNC', 'ESR1'), [('FPLX', 'ESR')]),
                          (('FPLX', 'ESR'), [('FPLX', 'ESR')]),
                          (('HGNC', 'SCN8A'),
                           [('FPLX', 'Cation_channels'),
                            ('FPLX', 'Voltage_gated_ion_channels')])])
def test_root_terms(famplex_graph, test_input, expected):
    assert famplex_graph.root_terms(*test_input) == expected


def test_root_terms_raises(famplex_graph):
    with pytest.raises(ValueError):
        famplex_graph.root_terms('HGNC', 'Gene')
