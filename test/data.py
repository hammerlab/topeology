from os import path

def data_path(name):
    return path.join(path.dirname(__file__), "data", name)
