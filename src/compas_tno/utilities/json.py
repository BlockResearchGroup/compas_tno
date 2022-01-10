import json


__all__ = [
    'update_json'
]


def update_json(infile):
    """Update JSON for new compas edge data information

    Parameters
    ----------
    infile : str
        The ``.json`` address of the file in the old format

    Returns
    -------
    outfile
        The address of the new, saved file.
    """

    outfile = ''.join(infile.split('.')[:-1]) + '_review.' + infile.split('.')[-1]

    file = open(infile)
    data = json.load(file)

    edge_data = {}
    for key in data['data']['edgedata']:
        a, b = key.split('-')
        newkey = str((int(a), int(b)))
        edge_data[newkey] = data['data']['edgedata'][key]

    data['data']['edgedata'] = edge_data

    with open(outfile, 'w') as savefile:
        json.dump(data, savefile)

    return outfile
