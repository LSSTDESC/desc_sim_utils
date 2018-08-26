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
        self._header = dict()

    def set_keywords(self, eimage_file, det_name, obs_md, phot_params):
        """Set the eimage keywords."""
        self._header['WCSAXES'] = 2
        self._header['MJD-OBS'] = obs_md.mjd.TAI
        self._header['DATE-OBS'] \
            = astropy.time.Time(obs_md.mjd.TAI, format='mjd').isot
        self._header['EXPTIME'] = phot_params.nexp*phot_params.exptime
        self._header['RADESYS'] = 'ICRS'
        self._header['EXTTYPE'] = 'IMAGE'
        self._header['FILTER'] = obs_md.bandpass
        self._header['RATEL'] = obs_md.pointingRA
        self._header['DECTEL'] = obs_md.pointingDec
        self._header['ROTANGLE'] = obs_md.rotSkyPos
        self._header['CHIPID'] = det_name
        self._header['OBSID'] = obs_md.OpsimMetaData['obshistID']
        self._header['OUTFILE'] = eimage_file

    def _writeHeader(self, header, bounds):
        """Re-implementation of the galsim.FitsWCS._writeHeader method."""
        header.update(super(ImSimWcs, self)._writeHeader(header, bounds))
        header.update(self._header)
        return header


class FitsWCS(dict):
    """Container class for WCS objects indexed by det_name."""
    def __init__(self, eimage_dir):
        super(FitsWCS, self).__init__()
        eimages = sorted(glob.glob(os.path.join(eimage_dir, 'lsst_e*fits*')))
        for eimage in eimages:
            chipid = fits.open(eimage)[0].header['CHIPID']
            det_name = "R:{},{} S:{},{}".format(*tuple(x for x in chipid
                                                       if x.isdigit()))
            self[det_name] = ImSimWcs(eimage)


def make_sensor_flat(det_name, wcs, counts_per_iter, niter, rng):
    """Pickleable function to use with multiprocessing module."""
    global camera_wrapper, phot_params, obs_md

    visit = obs_md.OpsimMetaData['obshistID']
    ccd_id = "R{}_S{}".format(det_name[2:5:2], det_name[8:11:2])
    outfile = 'lsst_e_{}_{}_{}.fits'.format(visit, ccd_id, obs_md.bandpass)
    if os.path.isfile(outfile + '.gz'):
        return

    print("running", det_name)
    gs_det = make_galsim_detector(camera_wrapper, det_name, phot_params, obs_md)

    desc.imsim.add_treering_info([gs_det])
    logger = desc.imsim.get_logger('INFO', name=det_name)

    wcs.set_keywords(outfile, det_name, obs_md, phot_params)
    my_flat = desc.imsim.make_flat(gs_det, counts_per_iter, niter, rng,
                                   logger=logger, wcs=wcs)
    my_flat.write(outfile)

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
        = desc.imsim.parsePhoSimInstanceFile(instcat, numRows=30)
    phot_params._exptime = 80000.
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
                args = (det_name, wcs, args.counts_per_iter, args.niter, rng)
                results.append(pool.apply_async(make_sensor_flat, args))
            pool.close()
            pool.join()
            for res in results:
                res.get()
