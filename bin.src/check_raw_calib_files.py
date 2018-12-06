#!/usr/bin/env python
"""
Script to check the statistics of raw calibration input files.
"""
import argparse
import pickle
import lsst.daf.persistence as dp
import desc.sim_utils

parser = argparse.ArgumentParser(description="QA script for simulated calibration raw files")
parser.add_argument('repo', type=str, help='data repository containing the raw data')
parser.add_argument('--visits', nargs='*', help='visits to be considered')
parser.add_argument('--nsigma', type=float, default=5,
                    help='Number of sigmas to use in the threshold derived from overscan data to exclude bright pixels')
parser.add_argument('--outfile', type=str, default='raw_amp_stats.pkl',
                    help='Name of output file to contain the results')
parser.add_argument('--processes', type=int, default=1,
                    help='Number of cores to use')
args = parser.parse_args()

butler = dp.Butler(args.repo)

if args.visits[0].startswith('@'):
    with open(args.visits[0][1:], 'r') as fd:
        visits = [int(_.strip()) for _ in fd if not _.startswith('#')]
else:
    visits = [int(_) for _ in args.visits]

results = desc.sim_utils.process_visits(butler, visits, nsigma=args.nsigma,
                                        processes=args.processes)

with open(args.outfile, 'wb') as output:
    pickle.dump(results, output)
