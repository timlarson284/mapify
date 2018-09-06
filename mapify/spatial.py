"""
Functions handling spatial coordinates, and reading/writing rasters.
"""

import os
from typing import Tuple

from osgeo import gdal
import numpy as np


def createtif(path: str, rows: int, cols: int, affine: tuple,
              datatype: int, proj: str, bands: int) -> gdal.Dataset:
    """
    Create a GeoTif data set to work with.
    """
    if os.path.exists(path):
        os.remove(path)

    ds = (gdal
          .GetDriverByName('GTiff')
          .Create(path, cols, rows, bands, datatype))

    ds.SetGeoTransform(affine)
    ds.SetProjection(proj)

    return ds


def writedata(ds: gdal.Dataset, data: np.ndarray, offx: int, offy: int) -> None:
    """
    Write a chip of data to the given data set.
    """
    ds.GetRasterBand(1).WriteArray(data.reshape(100, 100), offx, offy)
    return


def transform_geo(x: float, y: float, affine: tuple) -> Tuple(int, int):
    """
    Perform the affine transformation from a geospatial coordinate to row/col
    space.
    """
    # Spelled out for clarity
    col = (x - affine[0] - affine[3] * affine[2]) / affine[1]
    row = (y - affine[3] - affine[0] * affine[4]) / affine[5]

    return int(row), int(col)


def transform_rc(row: int, col: int, affine: tuple) -> Tuple(int, int):
    """
    Perform the affine transformation from a row/col coordinate to a geospatial
    space.
    """
    # Spelled out for clarity
    x = affine[0] + col * affine[1] + row * affine[2]
    y = affine[3] + col * affine[4] + row * affine[5]

    return x, y


def determine_hv(x: float, y: float, aff: tuple) -> Tuple[int, int]:
    """
    Determine the ARD tile H/V that contains the given coordinate.
    """
    return transform_geo(x, y, aff)[::-1]
