
import json
import compas_tno


__all__ = [
    'callback_save_json',
    'callback_create_json',
]


def callback_save_json(xopt):
    """Save the iteration in a json file"""

    DATA_FILENAME = compas_tno.get('output_fixed.json')

    with open(DATA_FILENAME, mode='r', encoding='utf-8') as f:
        data = json.load(f)

    i = len(data['iterations'])
    data['iterations'][i] = list(xopt.flatten())

    with open(DATA_FILENAME, mode='w', encoding='utf-8') as f:
        json.dump(data, f)

    return


def callback_create_json():

    data = {'iterations': {}}

    DATA_FILENAME = compas_tno.get('output_fixed.json')

    with open(DATA_FILENAME, mode='w', encoding='utf-8') as f:
        json.dump(data, f)

    return
