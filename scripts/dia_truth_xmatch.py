#!/usr/bin/env python3.6
# -*- coding: UTF-8 -*-

"""This performs cross matching between DIA pipeline outputs and the run1.2
truth catalog.
"""

import sys
from pathlib import Path

# User installed packages
site_pckgs_path = '/global/u1/d/djp81/.local/lib/python3.6/site-packages'
gcr_catalogs_path = '/global/u1/d/djp81/public/gcr-catalogs'
sys.path.extend((site_pckgs_path, gcr_catalogs_path))

from astropy.table import Table, hstack, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u
from tqdm import tqdm

import GCRCatalogs
import lsst.daf.persistence as dafPersist

from .postage_stamps import get_valid_dataids


def get_diasrc_for_id(butler, dataid, bbox=None):
    """For a given data Id, return a table of all DIA sources

    Args:
        butler (Butler): A data access butler
        dataid   (dict): A dictionary of values uniquely identifying an image
        bbox     (BBox): Optionally return only targets within a bbox

    Returns:
       An astropy table
    """

    diasrc_cat = butler.get('deepDiff_diaSrc', dataId=dataid).asAstropy()

    diasrc_table = diasrc_cat['id', 'coord_ra', 'coord_dec']
    diasrc_table['ra'] = diasrc_table['coord_ra'].to('deg')
    diasrc_table['dec'] = diasrc_table['coord_dec'].to('deg')
    diasrc_table['visit'] = dataid['visit']
    diasrc_table['filter'] = dataid['filter']
    diasrc_table['detector'] = dataid['detector']

    diasrc_table.remove_columns(['coord_ra', 'coord_dec'])
    return diasrc_table


def get_truth_catalog(cat_name, bbox=None):
    """Load Id, RA, and Dec values from a truth catalog

    Args:
        cat_name (str): Name of the catalog to load
            (e.g. "dc2_truth_run1.2_variable_summary")
        bbox    (BBox): Optionally return only targets within a bbox

    Returns:
       An astropy table object with id, ra, and dec fields
    """

    # Load truth data
    truth_gcr = GCRCatalogs.load_catalog(cat_name)
    truth_data = truth_gcr.get_quantities(['uniqueId', 'ra', 'dec'])
    truth_table = Table(truth_data)

    truth_table.rename_column('uniqueId', 'id')
    truth_table['ra'].unit = u.degree
    truth_table['dec'].unit = u.degree
    return truth_table


def match_dataid(butler, dataid, truth_cat, radius):
    """Cross match DIA sources from a given data Id with the truth catalog

    Args:
        butler   (Butler): A data access butler
        dataid     (dict): A dictionary of values uniquely identifying an image
        truth_cat (Table): GCR truth catalog
        radius    (float): Match radius in arc-seconds (Default: 1)

    Returns:
        Cross match results as an astropy table
    """

    source_cat = get_diasrc_for_id(butler, dataid)
    source_skycoord = SkyCoord(ra=source_cat['ra'], dec=source_cat['dec'])
    truth_skycoord = SkyCoord(ra=truth_cat['ra'], dec=truth_cat['dec'])

    idx, d2d, d3d = truth_skycoord.match_to_catalog_sky(source_skycoord)
    matched_data = source_cat[idx]

    out_data = hstack([truth_cat, matched_data], table_names=['truth', 'src'])
    out_data['d2d'] = d2d

    # Only keep matches within radius
    sep_constraint = d2d < radius * u.arcsec
    return out_data[sep_constraint]


def match_dataid_list(butler, dataid_list, truth_cat, radius=1):
    """Cross match DIA sources from a list of data Ids with the truth catalog

    Args:
        butler          (Butler): A data access butler
        dataid_list (list[dict]): A list of data identifiers
        truth_cat        (Table): GCR truth catalog
        radius           (float): Match radius in arc-seconds (Default: 1)

    Returns:
        Cross match results as an astropy table
    """

    tables = []
    for dataid in dataid_list:
        tables.append(match_dataid(butler, dataid, truth_cat, radius=radius))

    return vstack(tables)


def run(diff_im_dir, out_dir):
    """Create postage stamps for all DIA sources

    Args:
        diff_im_dir   (str): Output directory of DIA pipeline
        out_dir       (str): Directory to save cross match results into
    """

    out_dir = Path(out_dir).resolve()

    tqdm.write('Initializing butler...')
    butler = dafPersist.Butler(diff_im_dir)

    tqdm.write('Checking for valid dataId values...')
    dataid_list = get_valid_dataids(butler, 'deepDiff_diaSrc')

    tqdm.write('Cross matching sources with truth catalog...')
    truth_cat = get_truth_catalog('dc2_truth_run1.2_variable_summary')
    xm_results = match_dataid_list(butler, dataid_list, truth_cat)

    out_path = out_dir / 'xmatch.csv'
    tqdm.write(f'Writing to {out_path}')
    xm_results.write(out_path, overwrite=True)
