#!/usr/bin/env python
"""
Use multiprocessing module to run multiple sensor-flats on KNL.
"""
import os
import glob
import warnings
import multiprocessing
from astropy.io import fits
from astropy._erfa import ErfaWarning
import astropy.time
import galsim
from lsst.afw import cameraGeom
from lsst.sims.GalSimInterface import LSSTCameraWrapper, make_galsim_detector
import desc.imsim

class ImSimWcs(galsim.fitswcs.GSFitsWCS):
    """
    Subclass of GSFitsWCS that enables the extra keywords needed by
    eimages to be set when writing a FITS file.
    """
    def __init__(self, *args, **kwds):
        super(ImSimWcs, self).__init__(*args, **kwds)
        self.header = dict()

    def set_keywords(self, eimage_file, det_name, obs_md, phot_params):
        """Set the eimage keywords."""
        self.header['WCSAXES'] = 2
        self.header['MJD-OBS'] = obs_md.mjd.TAI
        exptime = phot_params.nexp*phot_params.exptime
        self.header['EXPTIME'] = exptime
        date_obs = astropy.time.Time(obs_md.mjd.TAI, format='mjd')
        self.header['DATE-OBS'] = date_obs.isot
        date_end = date_obs + astropy.time.TimeDelta(exptime, format='sec')
        self.header['DATE-END'] = date_end.isot
        self.header['AIRMASS'] = 0
        self.header['RADESYS'] = 'ICRS'
        self.header['EXTTYPE'] = 'IMAGE'
        self.header['FILTER'] = obs_md.bandpass
        self.header['RATEL'] = obs_md.pointingRA
        self.header['DECTEL'] = obs_md.pointingDec
        self.header['ROTANGLE'] = obs_md.rotSkyPos
        self.header['CHIPID'] \
            = "R{}_S{}".format(det_name[2:5:2], det_name[8:11:2])
        self.header['OBSID'] = obs_md.OpsimMetaData['obshistID']
        self.header['OUTFILE'] = eimage_file

    def _writeHeader(self, header, bounds):
        """Re-implementation of the galsim.FitsWCS._writeHeader method."""
        header.update(super(ImSimWcs, self)._writeHeader(header, bounds))
        header.update(self.header)
        return header


class FitsWCS(dict):
    """Container class for WCS objects indexed by det_name."""
    def __init__(self, eimage_dir):
        super(FitsWCS, self).__init__()
        eimage_files \
            = sorted(glob.glob(os.path.join(eimage_dir, 'lsst_e*fits*')))
        for eimage_file in eimage_files:
            with fits.open(eimage_file) as hdulist:
                chipid = hdulist[0].header['CHIPID']
            det_name = "R:{},{} S:{},{}".format(*tuple(x for x in chipid
                                                       if x.isdigit()))
            self[det_name] = ImSimWcs(eimage_file)
            self[det_name].eimage_file = eimage_file

def make_sensor_flat(det_name, wcs, counts_per_iter, niter, rng,
                     overwrite=True):
    """Pickleable function to use with multiprocessing module."""
    global camera_wrapper, phot_params, obs_md
    logger = desc.imsim.get_logger('INFO', name=det_name)

    visit = obs_md.OpsimMetaData['obshistID']
    ccd_id = "R{}_S{}".format(det_name[2:5:2], det_name[8:11:2])
    outfile = 'lsst_a_{}_{}_{}.fits'.format(visit, ccd_id, obs_md.bandpass)
    if not overwrite and os.path.isfile(outfile):
        logger.info("%s exists. Skipping.", outfile)
        return

    logger.info("running %s", det_name)
    gs_det = make_galsim_detector(camera_wrapper, det_name, phot_params, obs_md)

    desc.imsim.add_treering_info([gs_det])

    wcs.set_keywords(outfile, det_name, obs_md, phot_params)
    my_flat = desc.imsim.make_flat(gs_det, counts_per_iter, niter, rng,
                                   logger=logger, wcs=wcs)

    exptime = wcs.header['EXPTIME']
    with fits.open(wcs.eimage_file) as eimage:
        eimage[0].header.update(wcs.header)
        eimage[0].data = my_flat.array
        raw_image = desc.imsim.ImageSource(eimage[0].data, exptime, ccd_id,
                                           visit=visit)
        raw_image.eimage = eimage
        raw_image.eimage_data = eimage[0].data
        raw_image._read_pointing_info(None)
        raw_image.write_fits_file(outfile, image_type='FLAT')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('visit', type=int, help="visit number")
    parser.add_argument('--processes', type=int, default=1,
                        help='number of subprocesses to use with '
                        'the multiprocessing module. Default: 1')
    parser.add_argument('--counts_per_iter', type=int, default=4000,
                        help='counts per iteration. Default: 4000')
    parser.add_argument('--niter', type=int, default=20,
                        help='number of iterations. Default: 20')
    parser.add_argument('--header_dir', type=str, default=None,
                        help='directory of template headers for calibration '
                        'products. If None, then the default location '
                        'will be used.')
    parser.add_argument('--exptime', type=float, default=10.,
                        help='flat exposure time in seconds. Default: 10')
    args = parser.parse_args()

    header_dir = args.header_dir
    if header_dir is None:
        header_dir = os.path.join(os.environ['DESC_SIM_UTILS_DIR'],
                                  'data', 'calib_headers')
    fits_wcs = FitsWCS(header_dir)

    desc.imsim.read_config()

    instcat = os.path.join(os.environ['DESC_SIM_UTILS_DIR'],
                           'data', 'flat_instcat.txt')
    obs_md, phot_params, _ \
        = desc.imsim.parsePhoSimInstanceFile(instcat, (), numRows=30)
    phot_params._exptime = args.exptime
    obs_md.OpsimMetaData['obshistID'] = args.visit

    camera_wrapper = LSSTCameraWrapper()
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', 'ERFA function', ErfaWarning)

        with multiprocessing.Pool(processes=args.processes, maxtasksperchild=1)\
             as pool:
            results = []
            for i, det in enumerate(camera_wrapper.camera):
                det_name = det.getName()
                if det.getType() != cameraGeom.SCIENCE:
                    continue
                wcs = fits_wcs[det_name]
                rng = galsim.UniformDeviate(args.visit + i)
                fargs = (det_name, wcs, args.counts_per_iter, args.niter, rng)
                results.append(pool.apply_async(make_sensor_flat, fargs))
            pool.close()
            pool.join()
            for res in results:
                res.get()
