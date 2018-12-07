#!/usr/bin/env python
"""
Script to make bias or dark raw files for imSim.
"""
import os
import glob
import multiprocessing
import warnings
import numpy as np
from astropy.io import fits
from astropy._erfa import ErfaWarning
import astropy.time
import lsst.utils as lsstUtils
import desc.imsim

desc.imsim.get_config()

class DarkFrames:
    """Class to make dark frame exposures for the full LSST focalplane."""

    def __init__(self, cp_header_dir, visit0, exptime, mjd=59580,
                 ny_nx=(4072, 4000)):
        self.visit0 = visit0
        self.exptime = exptime
        self.mjd = mjd
        self.date_obs = astropy.time.Time(mjd, format='mjd')
        self.date_end \
            = self.date_obs + astropy.time.TimeDelta(exptime, format='sec')
        self.ny_nx = ny_nx
        self.header_files \
            = sorted(glob.glob(os.path.join(cp_header_dir, 'lsst_e*')))

    def _process_sensor(self, visit, template_file):
        global logger
        with fits.open(template_file) as eimage:
            chip_id = eimage[0].header['CHIPID']
            logger.info("  %s", chip_id)
            eimage[0].data = np.zeros(self.ny_nx, dtype=np.float)
            eimage[0].header['EXPTIME'] = self.exptime
            eimage[0].header['MJD-OBS'] = self.mjd
            eimage[0].header['DATE-OBS'] = self.date_obs.isot
            eimage[0].header['DATE-END'] = self.date_end.isot
            eimage[0].header['AIRMASS'] = 0
            eimage[0].header['OBSID'] = visit
            seed = self.random_seed(visit, chip_id)
            eimage[0].header['CR_SEED'] = seed
            eimage[0].data \
                = self.add_cosmic_rays(eimage[0].data, self.exptime, seed)
            raw_image \
                = desc.imsim.ImageSource(eimage[0].data, self.exptime, chip_id,
                                         visit=visit, logger=logger)
            raw_image.eimage = eimage
            raw_image.eimage_data = eimage[0].data
            raw_image._read_pointing_info(None)
            raw_file = 'lsst_a_{}_{}.fits'.format(visit, chip_id)
            raw_image.write_fits_file(raw_file)
        # Adjust header keywords for calibration product type.
        with fits.open(raw_file) as raw_frame:
            hdr = raw_frame[0].header
            hdr['OBSTYPE'] = 'BIAS' if self.exptime == 0 else 'DARK'
            hdr['IMGTYPE'] = hdr['OBSTYPE']
            raw_frame.writeto(raw_file, overwrite=True)

    def make(self, ivisit, processes=1):
        """Method to make a dark frame."""
        global logger
        visit = self.visit0 + ivisit

        logger.info("processing %s", visit)

        if processes == 1:
            # Processes each sensor serially.
            for template_file in self.header_files:
                self._process_sensor(visit, template_file)
        else:
            # Use multiprocessing.Pool.
            with multiprocessing.Pool(processes=processes, maxtasksperchild=1) \
                 as pool:
                results = []
                for template_file in self.header_files:
                    results.append(pool.apply_async(self._process_sensor,
                                                    (visit, template_file)))
                pool.close()
                pool.join()
                for res in results:
                    res.get()

    @staticmethod
    def add_cosmic_rays(imarr, exptime, seed):
        """
        Add cosmic rays to an array of pixels
        """
        config = desc.imsim.read_config()
        ccd_rate = config['cosmic_rays']['ccd_rate']
        if exptime == 0 or ccd_rate == 0:
            return imarr
        catalog = os.path.join(lsstUtils.getPackageDir('imsim'),
                               'data', 'cosmic_ray_catalog.fits.gz')
        crs = desc.imsim.CosmicRays.read_catalog(catalog, ccd_rate=ccd_rate)
        crs.set_seed(seed)
        return crs.paint(imarr, exptime=exptime)

    @staticmethod
    def random_seed(visit, chip_id):
        """
        Return a random seed generated from a hash of the visit number
        and chip id.
        """
        # Convert to obs_lsstSim CCD name so that seed is consistent
        # with imSim.
        detname = "R:{},{} S:{},{}".format(*[x for x in chip_id if x.isdigit()])
        return desc.imsim.CosmicRays.generate_seed(visit, detname)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('visit0', type=int, help='initial visit number')
    parser.add_argument('--nvisits', type=int, default=1,
                        help='number of visits. Default: 1')
    parser.add_argument('--exptime', type=float, default=0,
                        help='exposure time (s). Default: 0')
    parser.add_argument('--header_dir', type=str, default=None,
                        help='directory of template headers for calibration '
                        'products. If None, then the default location '
                        'will be used.')
    parser.add_argument('--log_level', type=str, default='INFO',
                        choices='DEBUG INFO WARN ERROR CRITICAL'.split())
    parser.add_argument('--processes', type=int, default=1,
                        help='number of subprocesses to use with '
                        'the multiprocessing module. Default: 1')
    args = parser.parse_args()

    header_dir = args.header_dir
    if header_dir is None:
        header_dir = os.path.join(os.environ['DESC_SIM_UTILS_DIR'],
                                  'data', 'calib_headers')

    logger = desc.imsim.get_logger(args.log_level)

    dark_frames = DarkFrames(header_dir, args.visit0, args.exptime)

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', 'ERFA function', ErfaWarning)

        for iv in range(args.nvisits):
            dark_frames.make(iv, processes=args.processes)
