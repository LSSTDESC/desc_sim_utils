#!/usr/bin/env python
import os
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

parser = argparse.ArgumentParser()
parser.add_argument('kernel_file', type=str,
                    help='Brighter-Fatter kernel file')
parser.add_argument('--outfile', type=str, default=None,
                    help='output png file for plot')
args = parser.parse_args()

outfile = args.outfile if args.outfile is not None \
          else os.path.basename(args.kernel_file.replace('.pkl', '.png'))

with open(args.kernel_file, 'rb') as fd:
    bfk = pickle.load(fd)
kernels = bfk.ampwiseKernels
fig = plt.figure()
for channel, data in kernels.items():
    plt.plot(range(len(data)), data[len(data)//2])
    print(channel, np.min(data))
plt.savefig(outfile)
