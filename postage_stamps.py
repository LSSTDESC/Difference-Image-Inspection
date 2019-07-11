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
from astropy.table import Table
from tqdm import tqdm
from astropy.coordinates import SkyCoord
from astropy import units as u

import GCRCatalogs
import lsst.afw.geom as afw_geom
import lsst.daf.persistence as dafPersist


def get_valid_ids(butler):
    """Return all dataId values that exist

    Args:
        butler (Butler): A data access butler

    Returns:
        A list of dataID values
    """

    subset = butler.subset('deepDiff_diaSrc')
    return [dr.dataId for dr in tqdm(subset) if dr.datasetExists()]


def get_diasrc_catalog(butler, data_id):
    """For a collection of data Ids, return a table of all DIA sources

    Args:
        butler (Butler): A data access butler
        data_id  (dict): A dictionary of values uniquely identifying an image

    Returns:
       An astropy table
    """

    diasrc_table = Table(
        names=['src_id', 'coord_ra', 'coord_dec', 'visit', 'filter',
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

    diasrc_table['coord_ra'].unit = u.rad
    diasrc_table['coord_dec'].unit = u.rad
    return diasrc_table


def get_truth_catalog():
    """Load catalog data from "dc2_truth_run1.2_variable_summary"

    Returns:
       A SourceCatalog object with id, coord_ra, and coord_dec fields
    """

    # Load truth data
    truth_gcr = GCRCatalogs.load_catalog('dc2_truth_run1.2_variable_summary')
    truth_data = truth_gcr.get_quantities([])

    truth_table = Table(
        data=[truth_data['uniqueId'], truth_data['ra'], truth_data['dec']],
        names=['src_id', 'coord_ra', 'coord_dec'],
        dtype=[np.int64, float, float]
    )

    truth_table['coord_ra'].unit = u.degree
    truth_table['coord_dec'].unit = u.degree
    return truth_table


def match_truth_catalog(source_ctlg, truth_ctlg, radius=1):
    """Cross match DIA sources with the truth catalog

    Args:
        source_ctlg (SourceCatalog): A list of dataID values
        truth_ctlg  (SourceCatalog): GCR truth catalog
        radius            (float): Match radius in arcseconds (Default: 1)

    Returns:
        Cross match results as an astropy table
    """

    source_skyc = SkyCoord(
        ra=source_ctlg['coord_ra'], dec=source_ctlg['coord_dec'])

    truth_skyc = SkyCoord(
        ra=truth_ctlg['coord_ra'], dec=truth_ctlg['coord_dec'])

    idx, d2d, d3d = truth_skyc.match_to_catalog_sky(source_skyc)

    sep_constraint = d2d < radius * u.arcsec
    out_table = source_ctlg[idx][sep_constraint]
    for col_name in ('uniqueId', 'ra', 'dec'):
        out_table[col_name] = truth_ctlg[col_name][idx]

    return out_table


def create_postage_stamp(butler, out_path, data_id, xpix, ypix, side_length,
                         dataset_type='deepDiff_differenceExp'):
    """Create a singe postage stamp and save it to file

    Args:
        butler     (Butler): A data access butler
        out_path      (str): Output path of fits file
        data         (dict): A valid data identifier
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
        file_pattern = '{visit:08d}-{filter}-{detector:03d}_{source_id}.fits'

        sources = butler.get('deepDiff_diaSrc', dataId=data_id)
        for source in tqdm(sources, position=1, desc='Sources'):
            out_path = out_dir / file_pattern.format(
                source_id=source['id'], **data_id)

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
    data_id_list = get_valid_ids(butler)

    tqdm.write('Building source catalog...')
    src_ctlg = get_diasrc_catalog(butler, data_id_list)

    tqdm.write('Building truth catalog...')
    truth_ctlg = get_truth_catalog()

    tqdm.write('Cross matching with truth catalog...')
    xm_results = match_truth_catalog(src_ctlg, truth_ctlg)
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
    if not out.exists():
        out.mkdir(parents=True, exist_ok=True)

    run(args.diff_dir, out, args.cutout_size)
