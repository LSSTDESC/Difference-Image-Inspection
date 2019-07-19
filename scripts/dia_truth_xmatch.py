#!/usr/bin/env python3.6
# -*- coding: UTF-8 -*-

"""This performs cross matching between DIA pipeline outputs and the run1.2
truth catalog.
"""

import argparse
import sys
from pathlib import Path

# User installed packages
site_pckgs_path = '/global/u1/d/djp81/.local/lib/python3.6/site-packages'
gcr_catalogs_path = '/global/u1/d/djp81/public/gcr-catalogs'
sys.path.extend((site_pckgs_path, gcr_catalogs_path))

import numpy as np
from astropy.table import Table, hstack, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u
from tqdm import tqdm

import GCRCatalogs
import lsst.daf.persistence as dafPersist
import lsst.afw.geom as afwGeom
from lsst.geom import SpherePoint

from .postage_stamps import get_valid_dataids


def get_diasrc_catalog(butler, dataid):
    """For a given data Id, return a table of all DIA sources

    Args:
        butler (Butler): A data access butler
        dataid   (dict): A dictionary of values uniquely identifying an image

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


def get_truth_catalog(cat_name):
    """Load Id, RA, and Dec values from a truth catalog

    Args:
        cat_name (str): Name of the catalog to load
            (e.g. "dc2_truth_run1.2_variable_summary")

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


def match_catalogs(catalog1, catalog2, radius):
    """positional cross match to find sources from catalog1 in catalog2

    Args:
        catalog1 (Table): A dictionary of values uniquely identifying an image
        catalog2 (Table): GCR truth catalog
        radius   (float): Match radius in arc-seconds (Default: 1)

    Returns:
        Cross match results as an astropy table
    """

    source_skycoord = SkyCoord(ra=catalog1['ra'], dec=catalog1['dec'])
    truth_skycoord = SkyCoord(ra=catalog2['ra'], dec=catalog2['dec'])

    idx, d2d, d3d = truth_skycoord.match_to_catalog_sky(source_skycoord)
    matched_data = catalog1[idx]

    out_data = hstack([catalog2, matched_data], table_names=['truth', 'src'])
    out_data['d2d'] = d2d

    # Only keep matches within radius
    sep_constraint = d2d < radius * u.arcsec
    return out_data[sep_constraint]


def get_bbox_for_ids(butler, dataid_list, dataset_type):
    """Return bbox objects for a lost of all data Id values

    Args:
        butler          (Butler): A data access butler
        dataid_list (list[dict]): A list of dataID values
        dataset_type       (str): The name of the data set to load

    Returns:
        A list with one BBox for each object Id
        A list with the WCS for each Id
    """

    bbox_list, wcs_list = [], []
    for dataid in dataid_list:
        image = butler.get(dataset_type, dataID=dataid)
        bbox_list.append(image.getBBox())
        wcs_list.append(image.getWcs())

    return bbox_list, wcs_list


def get_catalog_bbox_indices(ra, dec, bbox_list, wcs_list):
    """Crop a table of targets to only include targets in a list of boundaries

    Args:
        ra       (list[float]): Target RA coordinates
        dec      (list[float]): Target Dec coordinates
        bbox_list (list[BBox]): List of boundaries

    Returns:
        A list of booleans for whether each coordinate is in the boundaries
    """

    if not len(ra) == len(dec):
        raise ValueError(
            f'Unequal number of ra / dec coordinates: {len(ra)} != {len(dec)}')

    indices = []
    for ra_coord, dec_coord in zip(ra, dec):
        is_in = []
        for bbox, wcs in zip(bbox_list, wcs_list):
            radec = SpherePoint(ra_coord, dec_coord, afwGeom.degrees)
            xy = afwGeom.PointI(wcs.skyToPixel(radec))
            is_in.append(xy in bbox)

        indices.append(all(is_in))

    return np.ones(len(ra))


def match_dataid_list(butler, dataid_list, truth_cat, radius=1,
                      bbox_list=None, wcs_list=None):
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
    for dataid, bbox_this, wcs_this in zip(dataid_list, bbox_list, wcs_list):

        indices = get_catalog_bbox_indices(
            truth_cat['ra'], truth_cat['dec'], bbox_this, wcs_this)

        truth_cat_this = truth_cat[indices]
        source_cat = get_diasrc_catalog(butler, dataid)
        tables.append(
            match_catalogs(source_cat, truth_cat_this, radius=radius))

    # Todo: Add targets to table that were not matched and mask unmatched data

    return vstack(tables)


def run(diff_im_dir, out_path):
    """Create postage stamps for all DIA sources

    Args:
        diff_im_dir (str): Output directory of DIA pipeline
        out_path     (str): Directory to save cross match results into
    """

    out_path = Path(out_path).resolve()

    tqdm.write('Initializing butler...')
    butler = dafPersist.Butler(diff_im_dir)

    tqdm.write('Checking for valid dataId values...')
    dataid_list = get_valid_dataids(butler, 'deepDiff_diaSrc')

    tqdm.write('Getting tract / patch boundaries')
    sky_bbox, sky_wcs = get_bbox_for_ids(butler, dataid_list, 'deepCoadd_skymap')
    visit_bbox, visit_wcs = get_bbox_for_ids(butler, dataid_list, 'deepDiff_differenceExp')
    all_bboxes = np.array([sky_bbox, visit_bbox]).T
    all_wcs = np.array([sky_wcs, visit_wcs]).T

    tqdm.write('Cross matching sources with truth catalog...')
    truth_cat = get_truth_catalog('dc2_truth_run1.2_variable_summary')
    xm_results = match_dataid_list(
        butler, dataid_list, truth_cat, bbox_list=all_bboxes, wcs_list=all_wcs)

    out_path = out_path / 'xmatch.csv'
    tqdm.write(f'Writing to {out_path}')
    xm_results.write(out_path, overwrite=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create postage stamps from DIA pipeline results')

    parser.add_argument(
        '-r', '--repo',
        type=str,
        required=True,
        help='Path of DIA pipeline output directory')

    parser.add_argument(
        '-o', '--out_path',
        type=str,
        required=True,
        help='Where to write postage stamps')

    args = parser.parse_args()
    out = Path(args.out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    run(args.repo, out)
