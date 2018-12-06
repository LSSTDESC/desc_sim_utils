"""
Module for QA checks of simulated calibration raw files and products.
"""
import os
import glob
import pickle
import multiprocessing
import numpy as np
import matplotlib.pyplot as plt
import lsst.afw.geom as afw_geom
import lsst.afw.math as afw_math
import lsst.daf.persistence as dp


__all__ = ['compute_amp_stats', 'plot_pixel_dist', 'process_visits',
           'print_qa_stats']


def get_overscan_stats(image, amp_info):
    """
    Compute the mean and standard deviation of the serial overscan
    pixels for an amp.  Corrections are made for non-standard overscan
    sizes.

    Parameters
    ----------
    image: lsst.afw.image.Image
        The full imaging section of the amplifier.
    amp_info: lsst.afw.cameraGeom.AmpInfoRecord
        Object containing the amplifier pixel geometry.

    Returns
    -------
    (float, float): Tuple of the (mean, stdev) of the overscan pixel values.
    """
    oscan_corners = amp_info.getRawHorizontalOverscanBBox().getCorners()
    image_corners = image.getBBox().getCorners()
    # Create a bounding box using the upper left corner of the full
    # segment to guard against non-standard overscan region sizes.
    bbox = afw_geom.Box2I(oscan_corners[0], image_corners[2])
    oscan = image.Factory(image, bbox)
    stats = afw_math.makeStatistics(oscan, afw_math.MEAN | afw_math.STDEV)
    return stats.getValue(afw_math.MEAN), stats.getValue(afw_math.STDEV)


def compute_amp_stats(butler, dataId, nsigma=5):
    """
    Loop over all amplifiers in all detectors in a visit and compute the
    clipped mean and clipped stdev for each.

    Parameters
    ----------
    butler: daf.persistence.Butler
        Butler object pointing at the data repo with the raw files.
    dataId: dict
        dataId of desired amplifier data to be considered.nnnnnnnnnnn
    nsigma: float [5]
        Number of overscan sigmas over overscan mean to use as the
        threshold for excluding bright pixels, i.e., cosmic rays or
        bright defects.

    Returns
    -------
    tuple of np.arrays: means, stdevs, dataIds
    """
    print("processing", dataId)
    datarefs = butler.subset('raw_amp', dataId=dataId)
    means, stdevs, dataIds = [], [], []
    for dataref in list(datarefs):
        amp_exp = butler.get('raw_amp', dataref.dataId)
        det = amp_exp.getDetector()
        channel = dataref.dataId['channel']
        amp_info = list(det)[channel-1]
        oscan_mean, oscan_stdev \
            = get_overscan_stats(amp_exp.getImage(), amp_info)
        threshold = oscan_mean + nsigma*oscan_stdev
        image = amp_exp.getImage().Factory(amp_exp.getImage(),
                                           amp_info.getRawDataBBox())
        index = np.where(image.array < threshold)
        mean = np.mean(image.array[index])
        stdev = np.std(image.array[index])
        means.append(mean)
        stdevs.append(stdev)
        dataIds.append(dataref.dataId)
    return np.array(means), np.array(stdevs), np.array(dataIds)


def plot_pixel_dist(butler, dataId, x_range=None, bins=40, label=None):
    """
    Plot the distribution of pixel values in the imaging section of
    the specified amp.

    Parameters
    ----------
    butler: daf.persistence.Butler
        Butler object pointing at the data repo with the raw files.
    dataId: dict
        Data ID pointing at the desired amplifier.
    x_range: (float, float) [None]
        Range of pixel values to plot.  If None, then min and max
        of the pixel values are used.
    bins: int [40]
        Number of bins to use for the histogram.
    label: str [None]
        Label to apply to the plotted histogram.

    Returns
    -------
    np.array: The flattened array of pixel values.
    """
    amp_exp = butler.get('raw_amp', dataId)
    det = amp_exp.getDetector()
    channel = dataId['channel']
    amp = list(det)[channel-1]
    image = amp_exp.getImage().Factory(amp_exp.getImage(), amp.getRawDataBBox())
    data = image.array.flatten()
    if x_range is None:
        x_range = min(data), max(data)
    plt.hist(data, range=x_range, bins=bins, histtype='step', label=label)
    return data


def print_qa_stats(infile):
    """
    Print the mean and stdev stats for each visit in a results file.
    """
    with open(infile, 'rb') as fd:
        results = pickle.load(fd)
    print("   visit           avg(mean)              avg(sigma)")
    for visit in results:
        means, stdevs, _ = results[visit]
        avg_mean = np.mean(means)
        avg_sigma = np.mean(stdevs)
        print("{:8d}   {:>10.4f} +/- {:>4.3f}%  {:>10.4f} +/- {:>4.3f}%"
              .format(visit, avg_mean, 100.*np.std(means)/avg_mean,
                      avg_sigma, 100.*np.std(stdevs)/avg_sigma))
    return results


def process_visits(butler, visits, nsigma=5, processes=1):
    """
    Use multiprocessing module to run compute_amp_stats in parallel
    over visits.

    Parameters
    ----------
    butler: lsst.daf.persistence.Butler
        Data butler pointing to the desired data repo.
    visits: list
        List of visits (as ints)
    nsigma: float [5]
        Number of overscan sigmas over overscan mean to use as the
        threshold for excluding bright pixels, i.e., cosmic rays or
        bright defects.
    processes: int [1]
        Number of processes to run in multiprocessing pool. If 1, then
        just run serially.

    Returns
    -------
    dict: A dictionary of means, stdevs, dataIds, keyed by visit.
    """
    results = dict()
    if processes == 1:
        for visit in visits:
            results[visit] \
                = compute_amp_stats(butler, dict(visit=visit), nsigma=nsigma)
    else:
        with multiprocessing.Pool(processes=processes) as pool:
            for visit in visits:
                results[visit] = pool.apply_async(compute_amp_stats,
                                                  (butler, dict(visit=visit)),
                                                  dict(nsigma=nsigma))
            pool.close()
            pool.join()
        results = {visit: results[visit].get() for visit in results}
    return results
