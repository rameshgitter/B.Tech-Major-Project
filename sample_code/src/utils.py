import json
import os

def save_json(obj, path):
    with open(path, 'w') as f:
        json.dump(obj, f, indent=2)

def ensure_dir(path):
    d = os.path.dirname(path)
    if d and not os.path.exists(d):
        os.makedirs(d)
