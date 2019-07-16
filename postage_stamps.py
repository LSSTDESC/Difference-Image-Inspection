#!/usr/bin/env python3.6
# -*- coding: UTF-8 -*-

import argparse
import sys
from pathlib import Path

# Todo: Find a way to do this that doesn't involve my user directory
# User installed packages
site_pckgs_path = '/global/u1/d/djp81/.local/lib/python3.6/site-packages'
gcr_catalogs_path = '/global/u1/d/djp81/gcr-catalogs'
sys.path.extend((site_pckgs_path, gcr_catalogs_path))

import numpy as np
from astropy.table import Table, hstack, vstack
from tqdm import tqdm
from astropy.coordinates import SkyCoord
from astropy import units as u

import GCRCatalogs
import lsst.afw.geom as afw_geom
import lsst.daf.persistence as dafPersist
from lsst.pex.exceptions import LengthError


def get_valid_dataids(butler, dataset_type):
    """Return all data Id values that exist for a given data set

    Args:
        butler    (Butler): A data access butler
        dataset_type (str): The name of the data set to load

    Returns:
        A list of dataID values
    """

    subset = butler.subset(dataset_type)
    id_list = [dr.dataId for dr in tqdm(subset) if dr.datasetExists()]
    if not id_list:
        raise RuntimeError('No valid dataId found.')

    return id_list


def get_truth_catalog(catalog_name):
    """Load Id, RA, and Dec values from a truth catalog

    Args:
        catalog_name (str): Name of the catalog to load
            (e.g. "dc2_truth_run1.2_variable_summary")

    Returns:
       An astropy table object with id, ra, and dec fields
    """

    # Load truth data
    truth_gcr = GCRCatalogs.load_catalog(catalog_name)
    truth_data = truth_gcr.get_quantities(['uniqueId', 'ra', 'dec'])
    truth_table = Table(truth_data)

    truth_table.rename_column('uniqueId', 'id')
    truth_table['ra'].unit = u.degree
    truth_table['dec'].unit = u.degree
    return truth_table


def get_diasrc_for_id(butler, dataid):
    """For a given data Id, return a table of all DIA sources

    Args:
        butler (Butler): A data access butler
        dataid   (dict): A dictionary of values uniquely identifying an image

    Returns:
       An astropy table
    """

    diasrc_cat = butler.get('deepDiff_diaSrc', dataId=dataid).asAstropy()
    print(diasrc_cat.colnames)
    diasrc_table = diasrc_cat['id', 'coord_ra', 'coord_dec']
    diasrc_table['ra'] = diasrc_table['coord_ra'].to('deg')
    diasrc_table['dec'] = diasrc_table['coord_dec'].to('deg')
    diasrc_table.remove_columns(['coord_ra', 'coord_dec'])
    
    diasrc_table['visit'] = dataid['visit']
    diasrc_table['filter'] = dataid['filter']
    diasrc_table['detector'] = dataid['detector']
    return diasrc_table


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


def create_postage_stamp(butler, out_path, dataid, xpix, ypix, side_length,
                         dataset_type='deepDiff_differenceExp'):
    """Create a singe postage stamp and save it to file

    Args:
        butler     (Butler): A data access butler
        out_path      (str): Output path of fits file
        dataid       (dict): A valid data identifier
        xpix        (float): x pixel coordinate of cutout in degrees
        ypix        (float): y pixel coordinate of cutout in degrees
        side_length (float): Side length of cutout in pixels
        dataset_type  (str): Name of data set to create postage stamp for
    """

    cutout_size = afw_geom.ExtentI(side_length, side_length)
    xy = afw_geom.PointI(xpix, ypix)
    bbox = afw_geom.BoxI(xy - cutout_size // 2, cutout_size)

    try:
        cutout_image = butler.get(
            datasetType=f'{dataset_type}_sub',
            bbox=bbox,
            immediate=True,
            dataId=dataid)

    except LengthError:
        pass

    else:
        cutout_image.writeFits(str(out_path))


def save_stamps(butler, out_dir, dataid_list, cutout_size):
    """Create postage stamps for a list of data ids

    Creates postage stamps for all sources found for all data ids.

    Args:
        butler          (Butler): A data access butler
        out_dir           (Path): Output directory for fits files
        dataid_list (list[dict]): A list of dataID values
        cutout_size      (float): Side length of cutout in pixels
    """

    for dataid in tqdm(dataid_list, position=0, desc='Images'):
        sub_dir_name = '{visit:08d}-{filter}-{detector:03d}'.format(**dataid)
        out_sub_dir = out_dir / sub_dir_name
        out_sub_dir.mkdir(exist_ok=True)

        sources = butler.get('deepDiff_diaSrc', dataId=dataid)
        for source in tqdm(sources, position=1, desc='Sources'):
            out_path = out_sub_dir / f'{source["id"]}.fits'

            x_pix = source['base_NaiveCentroid_x']
            y_pix = source['base_NaiveCentroid_y']

            # noinspection PyTypeChecker
            create_postage_stamp(
                butler, out_path, dataid, x_pix, y_pix, cutout_size)


def run(diff_im_dir, out_dir, cutout_size):
    """Create postage stamps for all DIA sources

    Args:
        diff_im_dir   (str): Output directory of DIA pipeline
        out_dir       (str): Directory to save postave stamps into
        cutout_size (float): Side length of postage stamps in pixels
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

    tqdm.write('Creating postage stamps...')
    save_stamps(butler, out_dir, dataid_list, cutout_size)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create postage stamps from DIA pipeline results')

    parser.add_argument(
        '-r', '--repo',
        type=str,
        required=True,
        help='Path of DIA pipeline output directory')

    parser.add_argument(
        '-o', '--out_dir',
        type=str,
        required=True,
        help='Where to write postage stamps')

    parser.add_argument(
        '-s', '--cutout_size',
        type=int,
        required=True,
        help='Side length of cutout images')

    args = parser.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    run(args.repo, out, args.cutout_size)
