import csv
import sys


def _read_csv(fh, delimiter, quotechar):
    if sys.version_info.major < 3:
        csvreader = csv.reader(fh, delimiter=bytes(delimiter),
                               quotechar=bytes(quotechar))
        rows = [[cell.decode('utf-8') for cell in row] for row in csvreader]
    else:
        csvreader = csv.reader(fh, delimiter=delimiter, quotechar=quotechar)
        rows = [row for row in csvreader]
    return rows


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
