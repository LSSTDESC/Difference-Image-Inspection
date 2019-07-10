#!/usr/bin/env python3.6
# -*- coding: UTF-8 -*-

import sys
from pathlib import Path

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
    tqdm.write(str((type(source_ctlg), type(truth_ctlg))))
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


# Todo
def create_postage_stamp(butler, out_path, data_id, ra, dec, cutout_size):
    """Create a singe postage stamp and save it to file

    Args:
        butler     (Butler): A data access butler
        out_path     (Path): Output path of fits file
        data         (dict): A valid data identifier
        ra          (float): RA coordinate of cutout in degrees
        dec         (float): Dec coordinate of cutout in degrees
        cutout_size (float): Side length of cutout in pixels
    """

    pass


def save_stamps(butler, out_dir, data_id_list, cutout_size):
    """Create postage stamps for a list of data ids
    
    Creates postage stamps for all sources found for all data ids.
    
    Args:
        butler           (Butler): A data access butler
        out_dir            (Path): Output directory for fits files
        data_id_list (list[dict]): A list of dataID values
        cutout_size       (float): Side length of cutout in pixels
    """
    
    for data_id in data_id_list:
        file_pattern = (
            f"{data_id['visit']}-{data_id['filter']}-{data_id['raftName']}-"
            f"{data_id['detectorName']}-{data_id['detector']}-{{}}.fits"
        )

        for source in butler.get('deepDiff_diaSrc', dataId=data_id):
            out_path = out_dir / file_pattern.format(source['id'])
            create_postage_stamp(
                butler,
                out_path,
                data_id,
                source['coord_ra'],  # Todo: Convert to correct units
                source['coord_dec'],
                cutout_size
            )


def main(diff_im_dir, postage_output_dir):
    """Create postage stamps for all DIA sources
    
    Args:
        diff_im_dir        (str): Output directory of DIA pipeline
        postage_output_dir (str): Directory to save postave stamps into
    """

    postage_output_dir = Path(postage_output_dir).resolve()

    tqdm.write('Initializing butler...')
    butler = dafPersist.Butler(diff_im_dir)

    tqdm.write('Checking for valid dataId values...')
    # Todo: For development we skip the next line because it's slow
    # data_id_list = get_valid_ids(butler)
    data_id_list = [
        {'visit': 431306,
         'filter': 'u',
         'raftName': 'R10',
         'detectorName': 'S01',
         'detector': 28}
    ]

    tqdm.write('Building source catalog...')
    src_ctlg = get_diasrc_catalog(butler, data_id_list)

    tqdm.write('Building truth catalog...')
    truth_ctlg = get_truth_catalog()

    tqdm.write('Crossmatching with truth catalog...')
    xm_results = match_truth_catalog(src_ctlg, truth_ctlg)
    out_path = postage_output_dir / 'xmatch.csv'
    tqdm.write(f'Writing to {out_path}')
    xm_results.write(out_path, overwrite=True)

    tqdm.write('Creating postage stamps...')
    save_stamps(butler, data_id_list)


# Todo: Switch to CLI parsing for arguments
if __name__ == '__main__':
    main('/global/u1/d/djp81/test_imdiff/', './stamps')
