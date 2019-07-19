#!/usr/bin/env python3.6
# -*- coding: UTF-8 -*-

"""This script creates postage stamps for all sources identified by
the DIA pipeline.
"""

import argparse
import sys
from pathlib import Path

# User installed packages
site_pckgs_path = '/global/u1/d/djp81/.local/lib/python3.6/site-packages'
gcr_catalogs_path = '/global/u1/d/djp81/public/gcr-catalogs'
sys.path.extend((site_pckgs_path, gcr_catalogs_path))

from tqdm import tqdm

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

    # Clear console line
    tqdm.write('\n')


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
