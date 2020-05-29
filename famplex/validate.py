def find_invalid_rows(rows, expected_length):
    """Find all rows with incorrect length from input list

    Parameters
    ----------
    rows : list

    expected_length : int

    Returns
    -------
    list
        List of tuples of the form (index, observed_length) for each row in the
        input list with observed length different from the expected length.
    """
    bad_rows = []
    for index, row in enumerate(rows):
        if len(row) != expected_length:
            bad_rows.append((index, len(row)))
    return bad_rows


def find_grounding_map_mismatches(gmap_rows):
    """Find rows in grounding map with mismatching namespace/id amounts

    Parameters
    ----------
    gmap_rows : list
        List of rows from a grounding map csv file

    Returns
    -------
    list
        List of all rows in input where the number of namespace columns
        does not match the number of id columns
    """
    output = []
    for row in gmap_rows:
        keys = [entry for entry in row[1::2] if entry != '']
        values = [entry for entry in row[2::2] if entry != '']
        if len(keys) != len(values):
            output.append(row)
    return output
    
