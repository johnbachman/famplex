import pytest

from famplex.api import child_terms, parent_terms, ancestral_terms, \
    descendant_terms, individual_members, isa, partof, refinement_of, \
    dict_representation, equivalences


@pytest.mark.parametrize('test_input,rel_types,expected',
                         [(('FPLX', 'Voltage_gated_ion_channels'), None, []),
                          (('FPLX', 'SCN'), None,
                           [('FPLX', 'Sodium_channels'),
                            ('FPLX', 'Voltage_gated_ion_channels')]),
                          (('FPLX', 'TRP'), None,
                           [('FPLX', 'Voltage_gated_ion_channels')]),
                          (('HGNC', 'TRPA1'), None, [('FPLX', 'TRP')]),
                          (('HGNC', 'PRKAA1'), None,
                           [('FPLX', 'AMPK_A1B1G1'),
                            ('FPLX', 'AMPK_A1B1G2'),
                            ('FPLX', 'AMPK_A1B1G3'),
                            ('FPLX', 'AMPK_A1B2G1'),
                            ('FPLX', 'AMPK_A1B2G2'),
                            ('FPLX', 'AMPK_A1B2G3'),
                            ('FPLX', 'AMPK_alpha')]),
                          (('HGNC', 'PRKAA1'), ['isa'],
                           [('FPLX', 'AMPK_alpha')]),
                          (('HGNC', 'PRKAA1'), ['partof'],
                           [('FPLX', 'AMPK_A1B1G1'),
                            ('FPLX', 'AMPK_A1B1G2'),
                            ('FPLX', 'AMPK_A1B1G3'),
                            ('FPLX', 'AMPK_A1B2G1'),
                            ('FPLX', 'AMPK_A1B2G2'),
                            ('FPLX', 'AMPK_A1B2G3')])])
def test_parent_terms(test_input, rel_types, expected):
    assert parent_terms(*test_input, relation_types=rel_types) == expected


def test_parent_terms_raises():
    with pytest.raises(ValueError):
        parent_terms('HGNC', 'GENE')


@pytest.mark.parametrize('test_input,rel_types,expected',
                         [(('HGNC', 'TRPA1'), None, []),
                          (('FPLX', 'ESR'), None,
                           [('HGNC', 'ESR1'), ('HGNC', 'ESR2')]),
                          (('FPLX', 'AMPK'), None,
                           [('FPLX', 'AMPK_A1B1G1'),
                            ('FPLX', 'AMPK_A1B1G2'),
                            ('FPLX', 'AMPK_A1B1G3'),
                            ('FPLX', 'AMPK_A1B2G1'),
                            ('FPLX', 'AMPK_A1B2G2'),
                            ('FPLX', 'AMPK_A1B2G3'),
                            ('FPLX', 'AMPK_A2B1G1'),
                            ('FPLX', 'AMPK_A2B1G2'),
                            ('FPLX', 'AMPK_A2B1G3'),
                            ('FPLX', 'AMPK_A2B2G1'),
                            ('FPLX', 'AMPK_A2B2G2'),
                            ('FPLX', 'AMPK_A2B2G3'),
                            ('FPLX', 'AMPK_alpha'),
                            ('FPLX', 'AMPK_beta'),
                            ('FPLX', 'AMPK_gamma')]),
                          (('FPLX', 'AMPK'), ['isa'],
                          [('FPLX', 'AMPK_A1B1G1'),
                           ('FPLX', 'AMPK_A1B1G2'),
                           ('FPLX', 'AMPK_A1B1G3'),
                           ('FPLX', 'AMPK_A1B2G1'),
                           ('FPLX', 'AMPK_A1B2G2'),
                           ('FPLX', 'AMPK_A1B2G3'),
                           ('FPLX', 'AMPK_A2B1G1'),
                           ('FPLX', 'AMPK_A2B1G2'),
                           ('FPLX', 'AMPK_A2B1G3'),
                           ('FPLX', 'AMPK_A2B2G1'),
                           ('FPLX', 'AMPK_A2B2G2'),
                           ('FPLX', 'AMPK_A2B2G3')]),
                          (('FPLX', 'AMPK'), ['partof'],
                          [('FPLX', 'AMPK_alpha'),
                           ('FPLX', 'AMPK_beta'),
                           ('FPLX', 'AMPK_gamma')])])
def test_child_terms(test_input, rel_types, expected):
    assert child_terms(*test_input, relation_types=rel_types) == expected


def test_child_terms_raises():
    with pytest.raises(ValueError):
        child_terms('FPLX', 'Complex')


@pytest.mark.parametrize('test_input,rel_types,expected',
                         [(('FPLX', 'ESR'), None, []),
                          (('HGNC', 'ESR1'), None, [('FPLX', 'ESR')]),
                          (('HGNC', 'SCN8A'), None,
                          [('FPLX',
                            'Sodium_voltage_gated_channel_alpha_subunits'),
                           ('FPLX', 'SCN'), ('FPLX', 'Sodium_channels'),
                           ('FPLX', 'Voltage_gated_ion_channels'),
                           ('FPLX', 'Cation_channels')]),
                          (('HGNC', 'PRKAA1'), None,
                           [('FPLX', 'AMPK_A1B1G1'),
                            ('FPLX', 'AMPK_A1B1G2'),
                            ('FPLX', 'AMPK_A1B1G3'),
                            ('FPLX', 'AMPK_A1B2G1'),
                            ('FPLX', 'AMPK_A1B2G2'),
                            ('FPLX', 'AMPK_A1B2G3'),
                            ('FPLX', 'AMPK_alpha'),
                            ('FPLX', 'AMPK')]),
                          (('HGNC', 'PRKAA1'), ['isa'],
                           [('FPLX', 'AMPK_alpha')]),
                          (('HGNC', 'PRKAA1'), ['partof'],
                           [('FPLX', 'AMPK_A1B1G1'),
                            ('FPLX', 'AMPK_A1B1G2'),
                            ('FPLX', 'AMPK_A1B1G3'),
                            ('FPLX', 'AMPK_A1B2G1'),
                            ('FPLX', 'AMPK_A1B2G2'),
                            ('FPLX', 'AMPK_A1B2G3')])])
def test_ancestral_terms(test_input, rel_types, expected):
    assert ancestral_terms(*test_input,
                           relation_types=rel_types) == expected


def test_ancestral_terms_raises():
    with pytest.raises(ValueError):
        ancestral_terms('HGNC', 'GENE')


@pytest.mark.parametrize('test_input,rel_types,expected',
                         [(('HGNC', 'ESR1'), None, []),
                          (('FPLX', 'AMPK'), None,
                           [('FPLX', 'AMPK_A1B1G1'),
                            ('FPLX', 'AMPK_A1B1G2'),
                            ('FPLX', 'AMPK_A1B1G3'),
                            ('FPLX', 'AMPK_A1B2G1'),
                            ('FPLX', 'AMPK_A1B2G2'),
                            ('FPLX', 'AMPK_A1B2G3'),
                            ('FPLX', 'AMPK_A2B1G1'),
                            ('FPLX', 'AMPK_A2B1G2'),
                            ('FPLX', 'AMPK_A2B1G3'),
                            ('FPLX', 'AMPK_A2B2G1'),
                            ('FPLX', 'AMPK_A2B2G2'),
                            ('FPLX', 'AMPK_A2B2G3'),
                            ('FPLX', 'AMPK_alpha'),
                            ('FPLX', 'AMPK_beta'),
                            ('FPLX', 'AMPK_gamma'),
                            ('HGNC', 'PRKAA1'),
                            ('HGNC', 'PRKAB1'),
                            ('HGNC', 'PRKAG1'),
                            ('HGNC', 'PRKAG2'),
                            ('HGNC', 'PRKAG3'),
                            ('HGNC', 'PRKAB2'),
                            ('HGNC', 'PRKAA2')]),
                          (('FPLX', 'AMPK'), ['isa'],
                           [('FPLX', 'AMPK_A1B1G1'),
                            ('FPLX', 'AMPK_A1B1G2'),
                            ('FPLX', 'AMPK_A1B1G3'),
                            ('FPLX', 'AMPK_A1B2G1'),
                            ('FPLX', 'AMPK_A1B2G2'),
                            ('FPLX', 'AMPK_A1B2G3'),
                            ('FPLX', 'AMPK_A2B1G1'),
                            ('FPLX', 'AMPK_A2B1G2'),
                            ('FPLX', 'AMPK_A2B1G3'),
                            ('FPLX', 'AMPK_A2B2G1'),
                            ('FPLX', 'AMPK_A2B2G2'),
                            ('FPLX', 'AMPK_A2B2G3')]),
                          (('FPLX', 'AMPK'), ['partof'],
                           [('FPLX', 'AMPK_alpha'),
                            ('FPLX', 'AMPK_beta'),
                            ('FPLX', 'AMPK_gamma')])])
def test_descendant_terms(test_input, rel_types, expected):
    assert descendant_terms(*test_input,
                            relation_types=rel_types) == expected


def test_descendant_terms_raises():
    with pytest.raises(ValueError):
        descendant_terms('FPLX', 'Complex')


@pytest.mark.parametrize('test_input,rel_types,expected',
                         [(('FPLX', 'ESR'), None,
                           [('HGNC', 'ESR1'), ('HGNC', 'ESR2')]),
                          (('FPLX', 'AMPK'), None,
                           [('HGNC', 'PRKAA1'),
                            ('HGNC', 'PRKAA2'),
                            ('HGNC', 'PRKAB1'),
                            ('HGNC', 'PRKAB2'),
                            ('HGNC', 'PRKAG1'),
                            ('HGNC', 'PRKAG2'),
                            ('HGNC', 'PRKAG3')]),
                          (('FPLX', 'AMPK'), ['isa'],
                           [('FPLX', 'AMPK_A1B1G1'),
                            ('FPLX', 'AMPK_A1B1G2'),
                            ('FPLX', 'AMPK_A1B1G3'),
                            ('FPLX', 'AMPK_A1B2G1'),
                            ('FPLX', 'AMPK_A1B2G2'),
                            ('FPLX', 'AMPK_A1B2G3'),
                            ('FPLX', 'AMPK_A2B1G1'),
                            ('FPLX', 'AMPK_A2B1G2'),
                            ('FPLX', 'AMPK_A2B1G3'),
                            ('FPLX', 'AMPK_A2B2G1'),
                            ('FPLX', 'AMPK_A2B2G2'),
                            ('FPLX', 'AMPK_A2B2G3')]),
                          (('FPLX', 'AMPK'), ['partof'],
                           [('FPLX', 'AMPK_alpha'),
                            ('FPLX', 'AMPK_beta'),
                            ('FPLX', 'AMPK_gamma')])])
def test_individual_members(test_input, rel_types, expected):
    assert individual_members(*test_input,
                              relation_types=rel_types) == expected


def test_individual_members_raises():
    with pytest.raises(ValueError):
        individual_members('FPLX', 'Complex')


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
def test_isa(test_input, expected):
    assert isa(*test_input) == expected


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
def test_partof(test_input, expected):
    assert partof(*test_input) == expected


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
def test_refinement_of(test_input, expected):
    assert refinement_of(*test_input) == expected


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
def test_dict_representation(test_input, expected):
    assert dict_representation(*test_input) == expected


def test_dict_representation_raises():
    with pytest.raises(ValueError):
        dict_representation('FPLX', 'Complex')


def test_equivalences():
    assert set([('BEL', 'ESR Family')]) <= set(equivalences('ESR'))


def test_equivalences_raises():
    with pytest.raises(ValueError):
        equivalences('Complex')
