import matplotlib as mpl
import numpy as np

def get_matplotlib_version_tuple():
    return list(map(int, mpl.__version__.split('.')))

def get_numpy_version_tuple():
    return list(map(int, np.__version__.split('.')))
