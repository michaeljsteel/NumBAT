import re

import matplotlib as mpl
import numpy as np

def _version_tuple(version_str):
    # Keep the leading digits of each component so distro suffixes such as '3.10.7+dfsg1' parse.
    parts = []
    for part in version_str.split('.'):
        match = re.match(r'\d+', part)
        if match is None:
            break
        parts.append(int(match.group()))
    return parts

def get_matplotlib_version_tuple():
    return _version_tuple(mpl.__version__)

def get_numpy_version_tuple():
    return _version_tuple(np.__version__)
