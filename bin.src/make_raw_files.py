#!/usr/bin/env python
"""
Script to convert imSim eimage files to raw files using the
multiprocessing module.
"""
import os
import argparse
import warnings
import multiprocessing
from astropy._erfa import ErfaWarning
import desc.imsim

class WriteAmpFile:
    """Functor class to use in multiprocessing.apply_async callback."""
    def __init__(self, opsim_db=None):
        self.opsim_db = opsim_db

    def __call__(self, eimage_file, outdir='.'):
        global logger
        logger.info(os.path.basename(eimage_file))
        image_source \
            = desc.imsim.ImageSource.create_from_eimage(eimage_file,
                                                        opsim_db=self.opsim_db)
        image_source.write_fits_file(self.outfile(eimage_file, outdir=outdir))

    @staticmethod
    def outfile(eimage_file, outdir='.'):
        """Return the raw file name based on the eimage filename."""
        return os.path.join(outdir, os.path.basename(eimage_file).replace('lsst_e', 'lsst_a'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Script to convert imsim "
                                     "eimages to raw files using the "
                                     "multiprocessing module.")
    parser.add_argument('eimage_files', type=str, nargs="+",
                        help="eimage files to process")
    parser.add_argument('--opsim_db', type=str, default=None,
                        help="opsim db file.  This is only needed for older "
                        "eimage files that do not have pointing info in "
                        "the PHDU.")
    parser.add_argument('--processes', type=int, default=1,
                        help="Number of parallel processes to use.")
    parser.add_argument('--outdir', type=str, default='.',
                        help='output directory for raw files')
    parser.add_argument('--log_level', type=str, default='INFO',
                        choices='DEBUG INFO WARN ERROR CRITICAL'.split())
    args = parser.parse_args()

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)

    logger = desc.imsim.get_logger(args.log_level, name='make_raw_files')

    write_amp_file = WriteAmpFile(opsim_db=args.opsim_db)
    results = []
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', 'ERFA function', ErfaWarning)
        with multiprocessing.Pool(processes=args.processes, maxtasksperchild=1)\
             as pool:
            for item in args.eimage_files:
                outfile = write_amp_file.outfile(item, args.outdir)
                if os.path.isfile(outfile):
                    continue
                results.append(pool.apply_async(write_amp_file, (item,),
                                                dict(outdir=args.outdir)))

            pool.close()
            pool.join()
            for res in results:
                res.get()
