"""
Functions handling spatial coordinates, and reading/writing rasters.
"""

import os
from typing import Tuple

from osgeo import gdal
import numpy as np

from mapify.config import cu_tileaff


def createtif(path: str, rows: int, cols: int, affine: tuple,
              datatype: int, proj: str, bands: int) -> gdal.Dataset:
    """
    Create a GeoTif and return the data set to work with.
    If the file exists at the given path, this will attempt to remove it.

    Args:
        path: file path to create
        rows: number of rows
        cols: number of columns
        affine: gdal GeoTransform tuple
        datatype: gdal data type for the file
        proj: projection well known text
        bands: number of bands to create

    Returns:
        gdal data set for the file
    """
    if os.path.exists(path):
        os.remove(path)

    ds = (gdal
          .GetDriverByName('GTiff')
          .Create(path, cols, rows, bands, datatype))

    ds.SetGeoTransform(affine)
    ds.SetProjection(proj)

    return ds


def writedata(ds: gdal.Dataset, data: np.ndarray,
              col_off: int=0, row_off: int=0, band: int=1) -> None:
    """
    Write a chip of data to the given data set and band.

    Args:
        ds: gdal data set to write to
        data: data to write
        col_off: column offset to start writing data
        row_off: row offset to start writing data
        band: which band if it is a tiff-stack

    Returns:
        None
    """
    ds.GetRasterBand(band).WriteArray(data, col_off, row_off)
    return


def transform_geo(x: float, y: float, affine: tuple) -> Tuple(int, int):
    """
    Perform the affine transformation from a x/y coordinate to row/col
    space.

    Args:
        x: projected geo-spatial x coord
        y: projected geo-spatial y coord
        affine: gdal GeoTransform tuple

    Returns:
        containing pixel row/col
    """
    # Spelled out for clarity
    col = (x - affine[0] - affine[3] * affine[2]) / affine[1]
    row = (y - affine[3] - affine[0] * affine[4]) / affine[5]

    return int(row), int(col)


def transform_rc(row: int, col: int, affine: tuple) -> Tuple(int, int):
    """
    Perform the affine transformation from a row/col coordinate to projected x/y
    space.

    Args:
        row: pixel/array row number
        col: pixel/array column number
        affine: gdal GeoTransform tuple

    Returns:
        x/y coordinate
    """
    # Spelled out for clarity
    x = affine[0] + col * affine[1] + row * affine[2]
    y = affine[3] + col * affine[4] + row * affine[5]

    return x, y


def determine_hv(x: float, y: float, affine: tuple=cu_tileaff) -> Tuple[int, int]:
    """
    Determine the ARD tile H/V that contains the given coordinate.

    Args:
        x: projected geo-spatial x coord
        y: projected geo-spatial y coord
        affine: gdal GeoTransform tuple

    Returns:
        ARD tile h/v
    """
    return transform_geo(x, y, affine)[::-1]
