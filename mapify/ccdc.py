"""
Reading CCDC raw results, typically JSON or pickle formats for test data.
"""
import os
import json
import pickle
from typing import Tuple


def jsonpaths(root: str) -> list:
    """
    Create a list of file paths to files that end in .json in the given directory.

    Args:
        root: directory path

    Returns:
        list of JSON file paths
    """
    return [os.path.join(root, f)
            for f in os.listdir(root)
            if f[-5:] == '.json']


def picklepaths(root: str) -> list:
    """
    Create a list of file paths to files that end in .p in the given directory.

    Args:
        root: directory path

    Returns:
        list of pickle file paths
    """
    return [os.path.join(root, f)
            for f in os.listdir(root)
            if f[-5:] == '.p']


def pathcoords(path: str) -> Tuple[int, int]:
    """
    Pull the Chip X and Chip Y coords from the file path.

    Args:
        path: file path

    Returns:
        chip upper left x/y based on the file name
    """
    parts = os.path.split(path)[-1].split('_')
    return int(parts[1]), int(parts[2][:-5])


def loadjfile(path: str) -> dict:
    """
    Load a JSON formatted file into a dictionary.

    Args:
        path: file path

    Returns:
        dictionary representation of the JSON
    """
    return json.load(open(path, 'r'))


def loadjstr(string: str) -> dict:
    """
    Load a JSON formatted string into a dictionary.

    Args:
        string: JSON formatted string

    Returns:
        dictionary representation of the JSON
    """
    return json.loads(string)


def loadpfile(path: str) -> object:
    """
    Loads whatever object is contained in the pickle file.

    Args:
        path: file path

    Returns:
        some object
    """
    return pickle.load(open(path, 'rb'))
