#!/usr/bin/env python3.6
# -*- coding: UTF-8 -*-

import argparse
import sys
from copy import copy
from pathlib import Path

# Todo: Find a way to do this that doesn't involve my user directory
# User installed packages
site_pckgs_path = '/global/u1/d/djp81/.local/lib/python3.6/site-packages'
gcr_catalogs_path = '/global/u1/d/djp81/gcr-catalogs'
sys.path.extend((site_pckgs_path, gcr_catalogs_path))

import numpy as np
from astropy.table import Table, vstack, hstack
from tqdm import tqdm
from astropy.coordinates import SkyCoord
from astropy import units as u

import GCRCatalogs
import lsst.afw.geom as afw_geom
import lsst.daf.persistence as dafPersist


def get_valid_dataids(butler, **kwargs):
    """Return all dataId values that exist

    Args:
        butler (Butler): A data access butler

    Returns:
        A list of dataID values
    """

    subset = butler.subset('deepDiff_diaSrc')
    id_list = [dr.dataId for dr in tqdm(subset) if dr.datasetExists()]
    if not id_list:
        raise RuntimeError('No valid dataId found.')
    return id_list


def get_truth_catalog():
    """Load catalog data from "dc2_truth_run1.2_variable_summary"

    Returns:
       A SourceCatalog object with id, coord_ra, and coord_dec fields
    """

    # Load truth data
    truth_gcr = GCRCatalogs.load_catalog('dc2_truth_run1.2_variable_summary')
    truth_data = truth_gcr.get_quantities(['uniqueId', 'ra', 'dec'])

    truth_table = Table(
        data=[truth_data['uniqueId'], truth_data['ra'], truth_data['dec']],
        names=['truth_id', 'truth_ra', 'truth_dec'],
        dtype=[np.int64, float, float]
    )

    truth_table['truth_ra'].unit = u.degree
    truth_table['truth_dec'].unit = u.degree
    return truth_table


def get_diasrc_catalog(butler, data_id):
    """For a collection of data Ids, return a table of all DIA sources

    Args:
        butler (Butler): A data access butler
        data_id  (dict): A dictionary of values uniquely identifying an image

    Returns:
       An astropy table
    """

    diasrc_table = Table(
        names=['src_id', 'src_ra', 'src_dec', 'visit', 'filter',
               'detector'],
        dtype=[np.int64, float, float, np.int64, str, 'U10']
    )

    visit = data_id['visit']
    filter_ = data_id['filter']
    detector = data_id['detector']
    for source in butler.get('deepDiff_diaSrc', dataId=data_id):
        id_ = source['id']
        ra = source['coord_ra']
        dec = source['coord_dec']
        diasrc_table.add_row([id_, ra, dec, visit, filter_, detector])

    diasrc_table['src_ra'].unit = u.rad
    diasrc_table['src_dec'].unit = u.rad
    return diasrc_table


def match_dataid(butler, data_id, truth_ctlg, radius):
    """Cross match DIA sources from a given data Id with the truth catalog

    Args:
        butler     (Butler): A data access butler
        data_id      (dict): A dictionary of values uniquely identifying an image
        truth_ctlg  (Table): GCR truth catalog
        radius      (float): Match radius in arc-seconds (Default: 1)

    Returns:
        Cross match results as an astropy table
    """

    source_ctlg = get_diasrc_catalog(butler, data_id)
    source_skyc = SkyCoord(
        ra=source_ctlg['src_ra'], dec=source_ctlg['src_dec'])

    truth_skyc = SkyCoord(
        ra=truth_ctlg['truth_ra'], dec=truth_ctlg['truth_dec'])

    idx, d2d, d3d = truth_skyc.match_to_catalog_sky(source_skyc)

    out_table = copy(truth_ctlg)
    out_table.add_columns(list(source_ctlg[idx].itercols()))
    out_table['d2d'] = d2d

    # Only keep matches within radius
    sep_constraint = d2d < radius * u.arcsec
    return out_table[sep_constraint]


def match_dataid_list(butler, data_id_list, truth_ctlg, radius=1):
    """Cross match DIA sources from a list of data Ids with the truth catalog

    Args:
        butler           (Butler): A data access butler
        data_id_list (list[dict]): A list of data identifiers
        truth_ctlg        (Table): GCR truth catalog
        radius            (float): Match radius in arc-seconds (Default: 1)

    Returns:
        Cross match results as an astropy table
    """

    data_id_list = copy(data_id_list)
    first_id = data_id_list.pop()
    combined_table = match_dataid(butler, first_id, truth_ctlg, radius=radius)
    for data_id in data_id_list:
        match_table = match_dataid(butler, data_id, truth_ctlg, radius=radius)
        combined_table = vstack(combined_table, match_table)

    return combined_table


def create_postage_stamp(butler, out_path, data_id, xpix, ypix, side_length,
                         dataset_type='deepDiff_differenceExp'):
    """Create a singe postage stamp and save it to file

    Args:
        butler     (Butler): A data access butler
        out_path      (str): Output path of fits file
        data_id      (dict): A valid data identifier
        xpix        (float): x pixel coordinate of cutout in degrees
        ypix        (float): y pixel coordinate of cutout in degrees
        side_length (float): Side length of cutout in pixels
        dataset_type  (str): Name of data set to create postage stamp for
    """

    cutout_size = afw_geom.ExtentI(side_length, side_length)
    xy = afw_geom.PointI(xpix, ypix)
    bbox = afw_geom.BoxI(xy - cutout_size // 2, cutout_size)

    cutout_image = butler.get(
        datasetType=f'{dataset_type}_sub',
        bbox=bbox,
        immediate=True,
        dataId=data_id)

    cutout_image.writeFits(str(out_path))


def save_stamps(butler, out_dir, data_id_list, cutout_size):
    """Create postage stamps for a list of data ids

    Creates postage stamps for all sources found for all data ids.

    Args:
        butler           (Butler): A data access butler
        out_dir            (Path): Output directory for fits files
        data_id_list (list[dict]): A list of dataID values
        cutout_size       (float): Side length of cutout in pixels
    """

    for data_id in tqdm(data_id_list, position=0, desc='Images'):
        sub_dir_name = '{visit:08d}-{filter}-{detector:03d}'.format(**data_id)
        out_sub_dir = out_dir / sub_dir_name
        out_sub_dir.mkdir(exist_ok=True)

        sources = butler.get('deepDiff_diaSrc', dataId=data_id)
        for source in tqdm(sources, position=1, desc='Sources'):
            out_path = out_sub_dir / f'{source["id"]}.fits'

            x_pix = source['base_NaiveCentroid_x']
            y_pix = source['base_NaiveCentroid_y']

            # noinspection PyTypeChecker
            create_postage_stamp(
                butler,
                out_path,
                data_id,
                x_pix,
                y_pix,
                cutout_size
            )


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
    data_id_list = get_valid_dataids(butler)

    tqdm.write('Cross matching sources with truth catalog...')
    truth_ctlg = get_truth_catalog()
    xm_results = match_dataid_list(butler, data_id_list, truth_ctlg)

    out_path = out_dir / 'xmatch.csv'
    tqdm.write(f'Writing to {out_path}')
    xm_results.write(out_path, overwrite=True)

    tqdm.write('Creating postage stamps...')
    save_stamps(butler, out_dir, data_id_list, cutout_size)


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
