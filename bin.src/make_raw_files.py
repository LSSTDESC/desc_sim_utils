#!/usr/bin/env python
"""
Script to convert imSim eimage files to raw files using the
multiprocessing module.
"""
import os
import sys
import argparse
import multiprocessing
import desc.imsim

class WriteAmpFile:
    """Functor class to use in multiprocessing.apply_async callback."""
    def __init__(self, opsim_db=None):
        self.opsim_db = opsim_db

    def __call__(self, eimage_file, outdir='.'):
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
                        help="opsim db file.  If None, then the NERSC location "
                        "of minion_1016_desc_dithered_v4.db will be used.")
    parser.add_argument('--processes', type=int, default=1,
                        help="Number of parallel processes to use.")
    parser.add_argument('--outdir', type=str, default='.',
                        help='output directory for raw files')
    args = parser.parse_args()

    opsim_db = args.opsim_db if args.opsim_db is not None else \
               '/global/projecta/projectdirs/lsst/groups/SSim/DC2/minion_1016_desc_dithered_v4.db'

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    write_amp_file = WriteAmpFile(opsim_db=opsim_db)
    results = []
    with multiprocessing.Pool(processes=args.processes, maxtasksperchild=1) \
         as pool:
        for item in args.eimage_files:
            outfile = write_amp_file.outfile(item, args.outdir)
            if os.path.isfile(outfile):
                continue
            print("processing", os.path.basename(item))
            sys.stdout.flush()
            results.append(pool.apply_async(write_amp_file, (item,),
                                            dict(outdir=args.outdir)))

        pool.close()
        pool.join()
        for res in results:
            res.get()
