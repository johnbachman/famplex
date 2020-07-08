import csv


def load_csv(filename):
    """Load famplex csv file as list of rows

    Parameters
    ----------
    filename : str

    Returns
    -------
    rows : list
    """
    with open(filename) as f:
        csvreader = csv.reader(f, delimiter=str(u','),
                               lineterminator='\r\n',
                               quoting=csv.QUOTE_MINIMAL,
                               quotechar=str(u'"'))
        rows = [row for row in csvreader]
    return rows


def construct_grounding_map(rows):
    """Construct grounding map from rows in a grounding_map csv file

    Parameters
    ----------
    rows : list
        List of rows from a grounding map csv file. File should contain seven
        columns, the first of which is an agent text. The remaining columns
        contain namespace, id pairs, each pair occupying two columns. Some
        columns may be blank but in this case the row must be padded out with
        commas.

    Returns
    -------
    gmap : dict
        Dictionary mapping agent texts to INDRA style db_refs dicts. Each
        db_refs dict maps namespaces to ids.
    """
    gmap = {}
    for row in rows:
        text = row[0]
        db_refs = {'TEXT': text}
        db_refs.update({ns: id_ for ns, id_ in zip(row[1::2], row[2::2])})
        gmap[text] = db_refs if len(db_refs) > 1 else None
    return gmap


def update_id_prefixes(filename):
    """Return list of rows in grounding map with IDs corrected

    Parameters
    ----------
    filename : str
        Location of a grounding map csv file

    Returns
    -------
    list
        List of rows from grounding_map with GO, CHEBI, and CHEMBL IDs
        correctly prefixed with these respective namespaces. Leaves already
        correct IDs unchanged.
    """
    gm_rows = load_csv(filename)
    updated_rows = []
    for row in gm_rows:
        key = row[0]
        keys = [entry for entry in row[1::2]]
        values = [entry for entry in row[2::2]]
        if 'GO' in keys:
            go_ix = keys.index('GO')
            id_ = values[go_ix]
            if not id_.startswith('GO:'):
                values[go_ix] = 'GO:%s' % id_
        if 'CHEBI' in keys:
            chebi_ix = keys.index('CHEBI')
            id_ = values[chebi_ix]
            if not id_.startswith('CHEBI:'):
                values[chebi_ix] = 'CHEBI:%s' % id_
        if 'CHEMBL' in keys:
            chembl_ix = keys.index('CHEMBL')
            id_ = values[chembl_ix]
            if not id_.startswith('CHEMBL'):
                values[chembl_ix] = 'CHEMBL%s' % id_
        updated_row = [key]
        for pair in zip(keys, values):
            updated_row += pair
        updated_rows.append(updated_row)
    return updated_rows
