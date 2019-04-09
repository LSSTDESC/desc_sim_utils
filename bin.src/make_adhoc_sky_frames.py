#!/usr/bin/env python
"""
Script to make perfect sky frames.  These are intended to be used with
imSim data that has a flat response across the focalplane.
"""
import os
import glob
import argparse
import numpy as np
from astropy.io import fits
import lsst.utils

parser = argparse.ArgumentParser()
parser.add_argument('--outdir', type=str, default='sky_frames',
                    help="Output folder for sky frames for each band.")
parser.add_argument('--template_dir', type=str, default=None,
                    help="Folder containing sky frame template files. If None, "
                    "then use $DESC_SIM_UTILS_DIR/data/sky_frame_templates.")
args = parser.parse_args()

if args.template_dir is None:
    template_dir = os.path.join(lsst.utils.getPackageDir('desc_sim_utils'),
                                'data', 'sky_frame_templates')
else:
    template_dir = args.template_dir

template_files = sorted(glob.glob(os.path.join(template_dir, 'SKY*fits')))

for template in template_files:
    with fits.open(template) as hdus:
        calib_id = hdus[0].header['CALIB_ID']
        hdus[1].data = np.zeros(hdus[1].data.shape, dtype=np.float32)
        hdus[3].data = np.zeros(hdus[1].data.shape, dtype=np.float32)
        for band in 'ugrizy':
            hdus[0].header['CALIB_ID'] \
                = calib_id.replace('filter=u', 'filter={}'.format(band))
            outdir = os.path.join(args.outdir, band)
            os.makedirs(outdir, exist_ok=True)
            basename \
                = os.path.basename(template).replace('-u-', '-{}-'.format(band))
            outfile = os.path.join(outdir, basename)
            hdus.writeto(outfile, overwrite=True)
