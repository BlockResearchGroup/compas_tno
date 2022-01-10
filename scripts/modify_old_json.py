import json


infile = '/Users/mricardo/compas_dev/me/min_thk/dome/PAPER_CAS/dome_minthk.json'


def update_json(infile):
    outfile = ''.join(infile.split('.')[:-1]) + '_review.' + infile.split('.')[-1]

    file = open(infile)
    data = json.load(file)

    for key in data['data']:
        print(key)

    edge_data = {}
    for key in data['data']['edgedata']:
        print(key)
        a, b = key.split('-')
        newkey = str((int(a), int(b)))
        print(newkey)
        edge_data[newkey] = data['data']['edgedata'][key]

    data['data']['edgedata'] = edge_data

    with open(outfile, 'w') as savefile:
        json.dump(data, savefile)

    return outfile


