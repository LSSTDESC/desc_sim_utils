#!/usr/bin/env python
"""
Script to extract overlaps of sensor-visits with skymap tracts and patches
using a DM data repo registry.sqlite3 file and the DESC version of the
opsim db.
"""
import pickle
import sqlite3
import multiprocessing
import argparse
import numpy as np
from desc.sim_utils import process_registry_file, SkyMapPolygons

def process_rows(registry_file, opsim_db, sky_map_polygons, constraint,
                 row_bounds=None):
    """
    Process selected rows from a registry file.
    """
    df = process_registry_file(registry_file, opsim_db, sky_map_polygons,
                               constraint=constraint, row_bounds=row_bounds)
    if row_bounds is None:
        outfile = 'sensor_visit_info.pkl'
    else:
        outfile = 'sensor_visit_info_%06d_%06d.pkl' % row_bounds
    df.to_pickle(outfile)

parser = argparse.ArgumentParser()
parser.add_argument('registry_file', type=str, help='registry.sqlite3 file')
parser.add_argument('skyMap_file', type=str, help='skyMap.pickle file')
parser.add_argument('opsim_db_file', type=str, help='Opsim db file.')
parser.add_argument('--tract_info_file', type=str, default='tract_info.pkl',
                    help='Name of file to hold tract info extracted from the '
                    'skyMap file')
parser.add_argument('--constraint', type=str, default=None,
                    help='Constraint to apply on query of raw table in '
                    'the registry.sqlite3 file')
parser.add_argument('--processes', type=int, default=1,
                    help='number of concurrent processes to use')
args = parser.parse_args()

with open(args.skyMap_file, 'rb') as fd:
    sky_map = pickle.load(fd)
sky_map_polygons = SkyMapPolygons(sky_map, tract_info_file=args.tract_info_file)


if args.processes == 1:
    # Just run the processing function directly.
    process_rows(args.registry_file, args.opsim_db_file, sky_map_polygons,
                 args.constraint, row_bounds=None)
else:
    # Divide the rows to process among multiprocessing pool workers.
    #
    # Determine the total number of rows from the registry.sqlite file.
    with sqlite3.connect(args.registry_file) as conn:
        query = 'select count(*) from raw'
        if args.constraint is not None:
            query += ' where {}'.format(args.constraint)
        curs = conn.execute(query)
        nrows = [entry for entry in curs][0][0]

    row_array = [int(_) for _ in np.linspace(0, nrows, args.processes + 1)]

    with multiprocessing.Pool(processes=args.processes) as pool:
        workers = []
        for rmin, rmax in zip(row_array[:-1], row_array[1:]):
            args = (args.registry_file, args.opsim_db_file, sky_map_polygons,
                    args.constraint)
            kwargs = dict(row_bounds=(rmin, rmax))
            workers.append(pool.apply_async(process_rows, args, kwargs))
        pool.close()
        pool.join()
        results = [_.get() for _ in workers]
