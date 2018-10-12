#!/usr/bin/env python
"""
Script to repackage phosim amplifier files into single sensor
multi-extension FITS files, with an HDU per amplifier.
"""
import argparse
import desc.sim_utils

parser = argparse.ArgumentParser(description="Repackager for phosim amp files")
parser.add_argument('visit_dir', type=str, help="visit directory")
parser.add_argument('--out_dir', type=str, default=None,
                    help="output directory")
parser.add_argument('--image_type', type=str, default="SKYEXP",
                    choices='SKYEXP BIAS DARK FLAT'.split(),
                    help="image type")
parser.add_argument('--verbose', default=False, action='store_true',
                    help='print time to process the data for each sensor')
args = parser.parse_args()

repackager = desc.sim_utils.PhoSimRepackager()
repackager.process_visit(args.visit_dir, out_dir=args.out_dir,
                         image_type=args.image_type, verbose=args.verbose)
