
import json
import compas_tno


def callback_save_json(xopt, *args, **kwargs):
    """Save the iteration in a json file"""

    DATA_FILENAME = compas_tno.get('output.json')

    with open(DATA_FILENAME, mode='r', encoding='utf-8') as f:
        data = json.load(f)

    i = len(data['iterations'])
    data['iterations'][i] = list(xopt.flatten())

    with open(DATA_FILENAME, mode='w', encoding='utf-8') as f:
        json.dump(data, f)

    return


def callback_create_json():
    """Create a json to store the iterations of an optimisation"""

    data = {'iterations': {}}

    DATA_FILENAME = compas_tno.get('output.json')

    with open(DATA_FILENAME, mode='w', encoding='utf-8') as f:
        json.dump(data, f)

    return
