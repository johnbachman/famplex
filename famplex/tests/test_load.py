from famplex.load import load_grounding_map


def test_load_grounding_map():
    gm = load_grounding_map()
    for text, db_refs in gm.items():
        assert db_refs['TEXT'] == text
        assert '' not in db_refs
