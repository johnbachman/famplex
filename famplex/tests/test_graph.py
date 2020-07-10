import pytest


@pytest.fixture
def famplex_graph():
    from famplex import FamplexGraph
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
                         [(('FPLX', 'Voltage_gated_ion_channels'), []),
                          (('FPLX', 'SCN'),
                           [('FPLX', 'Sodium_channels'),
                            ('FPLX', 'Voltage_gated_ion_channels')]),
                          (('FPLX', 'TRP'),
                           [('FPLX', 'Voltage_gated_ion_channels')]),
                          (('HGNC', 'TRPA1'), [('FPLX', 'TRP')])])
def test_parent_terms(famplex_graph, test_input, expected):
    assert famplex_graph.parent_terms(*test_input) == expected


@pytest.mark.parametrize('test_input,expected',
                         [(('HGNC', 'TRPA1'), []),
                          (('FPLX', 'ESR'),
                           [('HGNC', 'ESR1'), ('HGNC', 'ESR2')])])
def test_child_terms(famplex_graph, test_input, expected):
    assert famplex_graph.child_terms(*test_input) == expected


@pytest.mark.parametrize('test_input,expected',
                         [(('HGNC', 'ESR1'), [('FPLX', 'ESR')]),
                          (('FPLX', 'ESR'), [('FPLX', 'ESR')]),
                          (('HGNC', 'SCN8A'),
                           [('FPLX', 'Cation_channels'),
                            ('FPLX', 'Voltage_gated_ion_channels')])])
def test_root_terms(famplex_graph, test_input, expected):
    assert famplex_graph.root_terms(*test_input) == expected


@pytest.mark.parametrize('test_input,expected',
                         [(('FPLX', 'ESR'), []),
                          (('HGNC', 'ESR1'), [('FPLX', 'ESR')]),
                          (('HGNC', 'SCN8A'),
                          [('FPLX',
                            'Sodium_voltage_gated_channel_alpha_subunits'),
                           ('FPLX', 'SCN'), ('FPLX', 'Sodium_channels'),
                           ('FPLX', 'Voltage_gated_ion_channels'),
                           ('FPLX', 'Cation_channels')])])
def test_ancestral_terms(famplex_graph, test_input, expected):
    assert famplex_graph.ancestral_terms(*test_input) == expected


@pytest.mark.parametrize('test_input,expected',
                         [(('HGNC', 'ESR1'), [])])
def test_descendant_terms(famplex_graph, test_input, expected):
    assert famplex_graph.descendant_terms(*test_input) == expected


@pytest.mark.parametrize('test_input,expected',
                         [(('FPLX', 'ESR'),
                           [('HGNC', 'ESR1'), ('HGNC', 'ESR2')])])
def test_individual_members(famplex_graph, test_input, expected):
    assert famplex_graph.individual_members(*test_input) == expected


@pytest.mark.parametrize('test_input,expected',
                         [(('HGNC', 'ESR1', 'FPLX', 'ESR'), True),
                          (('FPLX', 'ESR', 'HGNC', 'ESR1'), False),
                          (('HGNC', 'SCN8A', 'FPLX', 'SCN'), True),
                          (('HGNC', 'PRKAB1', 'FPLX', 'AMPK_A2B1G1'), False),
                          (('FPLX', 'AMPK_A2B1G1', 'HGNC', 'PRKAB1'), False),
                          (('FPLX', 'AMPK_A2B1G1', 'FPLX', 'AMPK'), True),
                          (('HGNC', 'PRKAA1', 'FPLX', 'AMPK'), False),
                          (('HGNC', 'DAP3', 'FPLX',
                            'Mitochondrial_Ribosome'), False),
                          (('HGNC', 'SCN8A', 'FPLX', 'MEK'), False),
                          (('HGNC', 'GENE', 'FPLX', 'MEK'), False)])
def test_isa(famplex_graph, test_input, expected):
    assert famplex_graph.isa(*test_input) == expected


@pytest.mark.parametrize('test_input,expected',
                         [(('HGNC', 'ESR1', 'FPLX', 'ESR'), False),
                          (('FPLX', 'ESR', 'HGNC', 'ESR1'), False),
                          (('HGNC', 'SCN8A', 'FPLX', 'SCN'), False),
                          (('HGNC', 'PRKAB1', 'FPLX', 'AMPK_A2B1G1'), True),
                          (('FPLX', 'AMPK_A2B1G1', 'HGNC', 'PRKAB1'), False),
                          (('FPLX', 'AMPK_A2B1G1', 'FPLX', 'AMPK'), False),
                          (('HGNC', 'PRKAA1', 'FPLX', 'AMPK'), False),
                          (('HGNC', 'DAP3', 'FPLX',
                            'Mitochondrial_Ribosome'), True),
                          (('HGNC', 'SCN8A', 'FPLX', 'MEK'), False),
                          (('HGNC', 'GENE', 'FPLX', 'MEK'), False)])
def test_partof(famplex_graph, test_input, expected):
    assert famplex_graph.partof(*test_input) == expected


@pytest.mark.parametrize('test_input,expected',
                         [(('HGNC', 'ESR1', 'FPLX', 'ESR'), True),
                          (('FPLX', 'ESR', 'HGNC', 'ESR1'), False),
                          (('HGNC', 'SCN8A', 'FPLX', 'SCN'), True),
                          (('HGNC', 'PRKAB1', 'FPLX', 'AMPK_A2B1G1'), True),
                          (('FPLX', 'AMPK_A2B1G1', 'HGNC', 'PRKAB1'), False),
                          (('FPLX', 'AMPK_A2B1G1', 'FPLX', 'AMPK'), True),
                          (('HGNC', 'PRKAA1', 'FPLX', 'AMPK'), True),
                          (('HGNC', 'DAP3', 'FPLX',
                            'Mitochondrial_Ribosome'), True),
                          (('HGNC', 'SCN8A', 'FPLX', 'MEK'), False),
                          (('HGNC', 'GENE', 'FPLX', 'MEK'), False)])
def test_refinement_of(famplex_graph, test_input, expected):
    assert famplex_graph.refinement_of(*test_input) == expected


@pytest.mark.parametrize('test_input,expected',
                         # Estrogen Receptor Family
                         [(('FPLX', 'ESR'),
                           # expected
                           {('FPLX', 'ESR'):
                            [({('HGNC', 'ESR1'): []}, 'isa'),
                             ({('HGNC', 'ESR2'): []}, 'isa')]}),
                          # Estrogen Receptor 1
                          (('HGNC', 'ESR1'),
                           # expected
                           {('HGNC', 'ESR1'): []}),
                          # MAP2K
                          (('FPLX', 'MAP2K'),
                           # expected
                           {('FPLX', 'MAP2K'):
                            [({('FPLX', 'MEK'):
                               [({('HGNC', 'MAP2K1'): []}, 'isa'),
                                ({('HGNC', 'MAP2K2'): []}, 'isa')]},
                              'isa'),
                             ({('HGNC', 'MAP2K3'): []}, 'isa'),
                             ({('HGNC', 'MAP2K4'): []}, 'isa'),
                             ({('HGNC', 'MAP2K5'): []}, 'isa'),
                             ({('HGNC', 'MAP2K6'): []}, 'isa'),
                             ({('HGNC', 'MAP2K7'): []}, 'isa')]})])
def test_dict_representation(famplex_graph, test_input, expected):
    assert famplex_graph.dict_representation(*test_input) == expected


def test_equivalences(famplex_graph):
    assert set([('BEL', 'ESR Family')]) <= \
        set(famplex_graph.equivalences('ESR'))
