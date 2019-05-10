#!/usr/bin/env python
import os
import glob
import pickle
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('input_file', type=str, help='Input B/F kernel file')
parser.add_argument('--channel', type=str, default='average',
                    help='Channel to use for detector-wide kernel. '
                    'If "average" (the default), then average the kernels '
                    'over all channels.')
parser.add_argument('--outdir', type=str, default='.',
                    help='Output directory for kernel files. [.]')
parser.add_argument('--template_dir', type=str, default=None,
                    help='Directory of detector-wide template kernel files.')
args = parser.parse_args()

def get_kernel(amp_kernels, channel):
    if channel == 'average':
        kstack = np.array(list(amp_kernels.ampwiseKernels.values()))
        return np.mean(kstack, axis=0)
    else:
        return amp_kernels.ampwiseKernels[channel]


template_dir = args.template_dir
if template_dir is None:
    template_dir = os.path.join('/global/cscratch1/sd/jchiang8/desc/dev/bf_kernels_run2.1i', 'detector_kernels')

det_files = sorted(glob.glob(os.path.join(template_dir, 'bfKernel*')))

with open(args.input_file, 'rb') as fd:
    amp_kernels = pickle.load(fd)
os.makedirs(args.outdir, exist_ok=True)
for i in range(189):
    with open(det_files[i], 'rb') as fd:
        det_kernel = pickle.load(fd)
    det_kernel.kernel[i] = get_kernel(amp_kernels, args.channel)
    outfile = os.path.join(args.outdir, os.path.basename(det_files[i]))
    with open(outfile, 'wb') as output:
        pickle.dump(det_kernel, output)
