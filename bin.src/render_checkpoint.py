#!/usr/bin/env python
"""
Script to render imSim checkpoint files as raw files.
"""
import argparse
from lsst.sims.GalSimInterface import (make_gs_interpreter,
                                       make_galsim_detector,
                                       LSSTCameraWrapper)
import desc.imsim


def render_checkpoint(instcat, checkpoint_file, create_eimage=False):
    """
    Render a checkpoint file as a raw file.
    """
    config = desc.imsim.read_config()
    RS_digits = [x for x in checkpoint_file.split('-')[-1] if x.isdigit()]
    det_name = 'R:{},{} S:{},{}'.format(*RS_digits)
    obs_md, phot_params, _ \
        = desc.imsim.parsePhoSimInstanceFile(instcat, [det_name], numRows=50,
                                             log_level='DEBUG')

    camera_wrapper = LSSTCameraWrapper()
    det = make_galsim_detector(camera_wrapper, det_name, phot_params, obs_md)
    gs_interpreter = make_gs_interpreter(obs_md, [det], None, None)
    gs_interpreter.checkpoint_file = checkpoint_file
    gs_interpreter.restore_checkpoint(camera_wrapper, phot_params, obs_md)

    desc.imsim.add_cosmic_rays(gs_interpreter, phot_params)
    full_well = int(config['ccd']['full_well'])
    desc.imsim.apply_channel_bleeding(gs_interpreter, full_well)

    visit = obs_md.OpsimMetaData['obshistID']
    name_root = f'lsst_a_{visit}'
    for name, gs_image in gs_interpreter.detectorImages.items():
        raw = desc.imsim.ImageSource.create_from_galsim_image(gs_image)
        outfile = '_'.join((name_root, name))
        raw.write_fits_file(outfile, compress=True)

    if create_eimage:
        name_root = f'lsst_e_{visit}'
        gs_interpreter.writeImages(nameRoot=name_root)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('instcat', type=str, help='instance catalog')
    parser.add_argument('checkpoint_file', type=str, help='checkpoint file')
    parser.add_argument('--create_eimage', default=False, action='store_true',
                        help='create an eimage file')
    args = parser.parse_args()

    render_checkpoint(args.instcat, args.checkpoint_file,
                      create_eimage=args.create_eimage)
