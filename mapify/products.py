"""
Functions for producing product values from CCDC results.
"""

import datetime as dt
from typing import Union, NamedTuple, Tuple, Sequence
from operator import attrgetter

import numpy as np
from osgeo import gdal

from mapify.config import dfc, chg_begining, chg_magbands, lc_map, nlcdxwalk


__ordbegin = dt.datetime.strptime(chg_begining, '%Y-%m-%d').toordinal()


class BandModel(NamedTuple):
    """
    Container for change detection spectral models.
    """
    name: str
    magnitude: float
    rmse: float
    intercept: float
    coeficients: Tuple


class CCDCModel(NamedTuple):
    """
    Container for the unified CCDC model.
    """
    start_day: int
    end_day: int
    break_day: Union(int, None)
    obs_count: Union(int, None)
    change_prob: Union(float, None)
    curve_qa: Union(int, None)
    bands: Union(Tuple, None)
    class_split: Union(str, None)
    class_probs1: Union(Tuple, None)
    class_probs2: Union(Tuple, None)
    class_vals: Union(Tuple, None)


def sortmodels(models: Sequence, key: str='start_day') -> list:
    """
    Sort a sequence of CCDC models based upon given key.

    Args:
        models: sequence of CCDC of namedtuples
        key: attribute to sort on

    Returns:
        sorted sequence
    """
    return sorted(models, key=attrgetter(key))


def classprobs(model: CCDCModel, ordinal: int) -> tuple:
    """
    Simple function to extract the class probabilities that go with the associated
    ordinal date. This function makes no assumptions on whether the date is
    actually contained within the segment, it simply does a date comparison
    against the class_split attribute if it exists.

    Args:
        model: classified CCDCmodel namedtuple
        ordinal: ordinal date

    Returns:
        class probabilities
    """
    if model.class_split and model.class_split >= ordinal:
        return model.class_probs2
    else:
        return model.class_probs1


def growth(model: CCDCModel, lc_mapping: dict=lc_map) -> bool:
    """
    Determine if the CCDC model represents a growth segment:
    tree -> grass/shrub, indicates decline
    or
    grass/shrub -> tree, indicates growth

    Args:
        model: classified CCDCmodel namedtuple
        lc_mapping: mapping of land cover names to values

    Returns:
        True if it is growth
    """
    if model.class_split:
        if np.argmax(model.class_probs1) == lc_mapping['grass']:
            return True

    return False


def decline(model: CCDCModel, lc_mapping: dict=lc_map) -> bool:
    """
    Determine if the CCDC model represents a growth segment:
    tree -> grass/shrub, indicates decline
    or
    grass/shrub -> tree, indicates growth

    Args:
        model: classified CCDCmodel namedtuple
        lc_mapping: mapping of land cover names to values

    Returns:
        True if it is decline
    """
    if model.class_split:
        if np.argmax(model.class_probs1) == lc_mapping['tree']:
            return True

    return False


def landcover(models: Sequence, ordinal: int, rank: int, dfcmap: dict=dfc,
              fill_begin: bool=True, fill_end: bool=True, fill_samelc: bool=True,
              fill_difflc: bool=True) -> int:
    """
    Given a sequence of CCDC models representing pixel history and an ordinal date,
    what is the Primary Land Cover value?

    Args:
        models: sorted sequence of CCDC namedtuples that represent the pixel history
        ordinal: date where we want cover for
        rank: which numeric rank to pull,
            1 - primary, 2- secondary, 3 - tertiary ...
        dfcmap: data format mapping, determines what values to assign for the
            various conditionals that could occur
        fill_begin: if the date falls before a known segment,
            use the first segment's value
        fill_end: if the date falls after all segments have ended,
            use the last segment's value
        fill_samelc: if the date falls between two segments,
            and they have the same class, use that class
        fill_difflc: if the date falls between two segments,
            if the date falls before the break date of the first segment, then use the first,
            if the date falls after the break date, then use the second

    Returns:
        primary land cover class value

    """
    # No model filling is provided by another function.
    if ordinal <= 0 or not models:
        return dfcmap['lc_insuff']

    # ord date before time series models -> cover back
    if fill_begin and ordinal < models[0].start_day:
        return models[0].class_vals[np.argsort(classprobs(models[0], ordinal))[rank]]

    # ord date after time series models -> cover forward
    if fill_end and ordinal > models[-1].end_day:
        return models[-1].class_vals[np.argsort(classprobs(models[-1], ordinal))[rank]]

    prev_end = 0
    prev_br = 0
    prev_class = 0
    for m in models:
        curr_class = m.class_vals[np.argsort(classprobs(m, ordinal))[rank]]
        # Date is contained within the model
        if m.start_day <= ordinal <= m.end_day:
            return curr_class
        # Same land cover fill
        elif fill_samelc and curr_class == prev_class and prev_end < ordinal < m.start_day:
            return curr_class
        # Different land cover fill, previous break -> current model
        elif fill_difflc and prev_br <= ordinal < m.start_day:
            return curr_class
        # Different land cover fill, model end -> break
        elif fill_difflc and m.end_day < ordinal < m.break_day:
            return curr_class

        prev_end = m.end_day
        prev_br = m.break_day
        prev_class = curr_class

    return dfcmap['lc_insuff']


def landcover_conf(models: Sequence, ordinal: int, rank: int, dfcmap: dict=dfc,
                   fill_begin: bool=True, fill_end: bool=True, fill_samelc: bool=True,
                   fill_difflc: bool=True) -> int:
    """
    Given a sequence of CCDC models representing pixel history and an ordinal date,
    what is the Primary Land Cover Confidence value?

    Args:
        models: sorted sequence of CCDC namedtuples that represent the pixel history
        ordinal: date where we want cover for
        rank: which numeric rank to pull,
            1 - primary, 2- secondary, 3 - tertiary ...
        dfcmap: data format mapping, determines what values to assign for the
            various conditionals that could occur
        fill_begin: if the date falls before a known segment,
            use the first segment's value
        fill_end: if the date falls after all segments have ended,
            use the last segment's value
        fill_samelc: if the date falls between two segments,
            and they have the same class, use that class
        fill_difflc: if the date falls between two segments,
            if the date falls before the break date of the first segment, then use the first,
            if the date falls after the break date, then use the second

    Returns:
        primary land cover confidence value
    """
    # No model filling is provided by another function.
    if ordinal <= 0 or not models:
        return 0

    # ord date before time series models -> cover back
    if fill_begin and ordinal < models[0].start_day:
        return dfcmap['lccf_back']

    # ord date after time series models -> cover forward
    if fill_end and ordinal > models[-1].end_day:
        if models[-1].change_prob == 1:
            return dfcmap['lccf_afterbr']

        return dfcmap['lccf_forwards']

    prev_end = 0
    prev_class = 0
    for m in models:
        probs = np.argsort(classprobs(models[-1], ordinal))
        curr_class = m.class_vals[probs[rank]]
        # Date is contained within the model
        if m.start_day <= ordinal <= m.end_day:
            # Annualized classification mucking jazz
            if growth(m):
                return dfcmap['lcc_growth']
            elif decline(m):
                return dfcmap['lcc_decline']
            return int(probs[rank] * 100)
        # Same land cover fill
        elif fill_samelc and curr_class == prev_class and prev_end < ordinal < m.start_day:
            return dfcmap['lccf_samelc']
        # Different land cover fill, prev model end -> current model start
        elif fill_difflc and prev_end <= ordinal < m.start_day:
            return dfcmap['lccf_difflc']

        prev_end = m.end_day
        prev_class = curr_class

    return 0


def crosswalk(inarr, xwalkmap: dict=nlcdxwalk) -> np.ndarray:
    """
    Cross-walks values in a data set to another set of values.

    Args:
        inarr: values to crosswalk
        xwalkmap: mapping of how to crosswalk

    Returns:
        np array of cross-walked values
    """

    outarr = np.copy(inarr)

    for old, new in xwalkmap.items():
        outarr[inarr == old] = new

    return outarr


def lc_nodatafill(lc_arr: np.ndarray, lcc_arr: np.ndarray, nlcd: np.ndarray,
                  xwalk: bool=True, xwalkmap: dict=nlcdxwalk,
                  dfcmap: dict=dfc) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fill areas without a model with NLCD data.

    Args:
        lc_arr: array of land cover values
        lcc_arr: array of land cover confidence values
        nlcd: array of NLCD values for same region
        xwalk: boolean if the NLCD values need to be cross-walked
        xwalkmap: mapping to on how to cross-walk the NLCD values
        dfcmap: data format mapping, determines what values to assign for the
            various conditionals that could occur

    Returns:
        filled land cover and land cover confidence arrays
    """
    if xwalk:
        nlcd = crosswalk(nlcd, xwalkmap)

    outlc = np.copy(lc_arr)
    outlcc = np.copy(lcc_arr)

    mask = outlc == dfcmap['lc_insuff']
    outlc[mask] = nlcd[mask]
    outlcc[mask] = dfcmap['lccf_nomodel']

    return outlc, outlcc


def lc_primary(models: Sequence, ordinal: int, dfcmap: dict=dfc,
               fill_begin: bool=True, fill_end: bool=True, fill_samelc: bool=True,
               fill_difflc: bool=True) -> int:
    return landcover(models, ordinal, 1, dfcmap,
                     fill_begin, fill_end, fill_samelc, fill_difflc)


def lc_secondary(models: Sequence, ordinal: int, dfcmap: dict=dfc,
                 fill_begin: bool=True, fill_end: bool=True, fill_samelc: bool=True,
                 fill_difflc: bool=True) -> int:
    return landcover(models, ordinal, 2, dfcmap,
                     fill_begin, fill_end, fill_samelc, fill_difflc)


def lc_primaryconf(models: Sequence, ordinal: int, dfcmap: dict=dfc,
                   fill_begin: bool=True, fill_end: bool=True,
                   fill_samelc: bool=True, fill_difflc: bool=True) -> int:
    return landcover_conf(models, ordinal, 1, dfcmap,
                          fill_begin, fill_end, fill_samelc, fill_difflc)


def lc_secondaryconf(models: Sequence, ordinal: int, dfcmap: dict=dfc,
                     fill_begin: bool=True, fill_end: bool=True,
                     fill_samelc: bool=True, fill_difflc: bool=True) -> int:
    return landcover_conf(models, ordinal, 2, dfcmap,
                          fill_begin, fill_end, fill_samelc, fill_difflc)


def lc_fromto(models: Sequence, ordinal: int) -> int:
    """
    Traditional from-to for the primary land cover.
    Assumes all Land Cover filling strategies are used.

    Args:
        models: sorted sequence of CCDC namedtuples that represent the pixel history
        ordinal: date where we want cover for

    Returns:
        fromto value

    """
    prev_yr = dt.date.fromordinal(ordinal)
    prev_yr = dt.date(year=prev_yr.year - 1, month=prev_yr.month, day=prev_yr.day)

    curr = landcover(models, ordinal, 1)
    prev = landcover(models, prev_yr.toordinal(), 1)

    if prev == curr:
        return curr

    return prev * 10 + curr


def chg_doy(models: Sequence, ordinal: int) -> int:
    """
    The day of year that a change happened, if a change happened in the
    same year as the ordinal given.

    Args:
        models: sorted sequence of CCDC namedtuples that represent the pixel history
        ordinal: date for a year that we want to know if change occurs in

    Returns:
        day of year or 0

    """
    if ordinal <= 0:
        return 0

    query_date = dt.date.fromordinal(ordinal)

    for m in models:
        if m.break_day <= 0:
            continue

        break_date = dt.date.fromordinal(m.break_day)

        if query_date.year == break_date.year and m.change_prob == 1:
            return break_date.timetuple().tm_yday

    return 0


def chg_mag(models: Sequence, ordinal: int, bands: Sequence=chg_magbands) -> float:
    """
    The spectral magnitude of the change (if one occurred) in the same
    year as the given ordinal.

    Args:
        models: sorted sequence of CCDC namedtuples that represent the pixel history
        ordinal: magnitude for a change, if it occurred in the year
        bands: spectral band names to perform the calculation over

    Returns:
        magnitude or 0

    """
    if ordinal <= 0:
        return 0

    query_date = dt.date.fromordinal(ordinal)

    for m in models:
        if m.break_day <= 0:
            continue

        break_date = dt.date.fromordinal(m.break_day)

        if query_date.year == break_date.year and m.change_prob == 1:
            mags = [b.magnitude for b in m.bands if b.name in bands]
            return np.linalg.norm(mags)

    return 0


def chg_modelqa(models: Sequence, ordinal: int) -> int:
    """
    Information on the quality of the curve fit the intercept the ordinal date.

    Args:
        models: sorted sequence of CCDC namedtuples that represent the pixel history
        ordinal: magnitude for a change, if it occurred in the year

    Returns:
        curve_qa or 0

    """
    if ordinal <= 0:
        return 0

    for m in models:
        if m.start_day <= ordinal <= m.end_day:
            return m.curve_qa

    return 0


def chg_seglength(models: Sequence, ordinal: int, ordbegin: int=__ordbegin) -> int:
    """
    How long, in days, has the current model been underway. This includes
    between or outside of CCD segments as well.

    Args:
        models: sorted sequence of CCDC namedtuples that represent the pixel history
        ordinal: magnitude for a change, if it occurred in the year
        ordbegin: when to start counting from

    Returns:
        number of days
    """
    if ordinal <= 0:
        return 0

    diff = [ordinal - ordbegin]
    for m in models:
        if ordinal > m.end_day:
            diff.append(ordinal - m.end_day)
        else:
            diff.append(ordinal - m.start_day)

    return min(filter(lambda x: x >= 0, diff), default=0)


def chg_lastbrk(models: Sequence, ordinal: int) -> int:
    """
    How long ago, in days, was the last spectral break.
    0 if before the first break.

    Args:
        models: sorted sequence of CCDC namedtuples that represent the pixel history
        ordinal: magnitude for a change, if it occurred in the year

    Returns:
        number of days
    """
    if ordinal <= 0:
        return 0

    diff = [(ordinal - m.break_day) for m in models if m.change_prob == 1]

    return min(filter(lambda x: x >= 0, diff), default=0)


def prodmap() -> dict:
    """
    Container function to hold the product mapping, showing which names map to
    which functions and other such information ...

    {'product name': [function, gdal data type],
     ...}

    Returns:
        product mapping
    """
    return {'Change': [chg_doy, gdal.GDT_UInt16],
            'LastChange': [chg_lastbrk, gdal.GDT_UInt16],
            'SegLength': [chg_seglength, gdal.GDT_UInt16],
            'ChangeMag': [chg_mag, gdal.GDT_Float32],
            'Quality': [chg_modelqa, gdal.GDT_Byte],
            'CoverPrim': [lc_primary, gdal.GDT_Byte],
            'CoverSec': [lc_secondary, gdal.GDT_Byte],
            'CoverConfPrim': [lc_primaryconf, gdal.GDT_Byte],
            'CoverConfSec': [lc_secondaryconf, gdal.GDT_Byte],
            'FromTo': [lc_fromto, gdal.GDT_Byte]}
