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

import GCRCatalogs
import lsst.afw.geom as afw_geom
import lsst.afw.table as afw_table
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


def get_diasrc_catalog(butler, data_id_list):
    """For a collection of data Ids, return a catalog of all DIA sources

    Args:
        butler           (Butler): A data access butler
        data_id_list (list[dict]): A list of dataID values

    Returns:
       A SourceCatalog object with id, coord_ra, and coord_dec fields
    """

    dia_schema = afw_table.SourceTable.makeMinimalSchema()
    dia_schema.addField('visit', type=np.int64, doc='Visit id')
    dia_schema.addField('filter', type=str, size=1)
    dia_schema.addField('raftName', type=str, size=3)
    dia_schema.addField('detectorName', type=str, size=6)
    dia_schema.addField('detector', type=np.int32, doc='detector')
    src_catalog = afw_table.SourceCatalog(dia_schema)

    for data_id in tqdm(data_id_list):
        visit = data_id['visit']
        filter_ = data_id['filter']
        raft_name = data_id['raftName']
        detector_name = data_id['detectorName']
        detector = data_id['detector']

        for source in butler.get('deepDiff_diaSrc', dataId=data_id):
            id_ = source['id']
            ra = afw_geom.Angle(source['coord_ra'], afw_geom.radians)
            dec = afw_geom.Angle(source['coord_dec'], afw_geom.radians)

            record = src_catalog.addNew()
            record.set('id', id_)
            record.set('coord_ra', ra)
            record.set('coord_dec', dec)
            record.set('visit', visit)
            record.set('filter', filter_)
            record.set('raftName', raft_name)
            record.set('detectorName', detector_name)
            record.set('detector', detector)

    return src_catalog


def get_truth_catalog():
    """Load catalog data from "dc2_truth_run1.2_variable_summary"

    Returns:
       A SourceCatalog object with id, coord_ra, and coord_dec fields
    """

    # Load truth data
    truth_gcr = GCRCatalogs.load_catalog('dc2_truth_run1.2_variable_summary')
    truth_data = truth_gcr.get_quantities(['uniqueId', 'ra', 'dec'])
    truth_data_zip = zip(truth_data['uniqueId'], truth_data['ra'],
                         truth_data['dec'])

    # Create an empty catalog and populate
    # minimal schema only contains the `id`, `coord_ra`, and `coord_dec` fields
    truth_schema = afw_table.SourceTable.makeMinimalSchema()
    truth_catalog = afw_table.SourceCatalog(truth_schema)

    # Populate catalog
    data_iter = tqdm(truth_data_zip, total=len(truth_data['uniqueId']))
    for id_, ra, dec in data_iter:
        record = truth_catalog.addNew()
        record.set('id', id_)
        record.set('coord_ra', afw_geom.Angle(ra, afw_geom.degrees))
        record.set('coord_dec', afw_geom.Angle(dec, afw_geom.degrees))

    return truth_catalog


def match_truth_catalog(source_ctlg, truth_ctlg, radius=1):
    """Cross match DIA sources with the truth catalog

    Args:
        source_ctlg (SourceCatalog): A list of dataID values
        truth_ctlg  (SourceCatalog): GCR truth catalog
        radius            (float): Match radius in arcseconds (Default: 1)

    Returns:
        Cross match results as an astropy table
    """

    radius = afw_geom.Angle(radius, afw_geom.arcseconds)
    matches = afw_table.matchRaDec(source_ctlg, truth_ctlg, radius)

    out_table = Table(
        names=['src_id', 'src_ra', 'src_dec', 'truth_id', 'truth_ra',
               'truth_dec', 'sep', 'visit', 'filter', 'raftName',
               'detectorName', 'detector'],
        dtype=[int, float, float, int, float, float, float, int, 'U100',
               'U100', 'U100', int]
    )

    for match in matches:
        src_id = match.first['id']
        src_ra = match.first['coord_ra']
        src_dec = match.first['coord_dec']
        truth_id = match.second['id']
        truth_ra = match.second['coord_ra']
        truth_dec = match.second['coord_dec']
        sep = np.degrees(match.distance) * 3600 * 1000
        visit = match.first['visit']
        filter_ = match.first['filter']
        raft_name = match.first['raftName']
        detector_name = match.first['detectorName']
        detector = match.first['detector']

        out_table.add_row(
            [src_id, src_ra, src_dec, truth_id,
             truth_ra, truth_dec, sep, visit, filter_, raft_name,
             detector_name, detector])

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
        dataset_type=f'{dataset_type}_sub',
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
            out_path = out_dir / file_pattern.format(source_id=source['id'], **data_id)
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
