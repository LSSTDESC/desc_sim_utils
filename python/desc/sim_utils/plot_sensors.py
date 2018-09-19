"""
Functions to plot Run2.0 simulation region and sensors that intersect it.
"""
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib import patches
import lsst.afw.geom as afw_geom
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import lsst.sims.coordUtils
import desc.imsim


__all__ = ['make_patch', 'plot_sensors', 'plot_survey_region']


def make_patch(vertexList, wcs=None):
    """
    Return a Path in sky coords from vertex list in pixel coords.

    Parameters
    ----------
    vertexList: list of coordinates
        These are the corners of the region to be plotted either in
        pixel coordinates or sky coordinates.
    wcs: lsst.afw.geom.skyWcs.skyWcs.SkyWcs [None]
        The WCS object used to convert from pixel to sky coordinates.

    Returns
    -------
    matplotlib.path.Path: The encapsulation of the vertex info that
        matplotlib uses to plot a patch.
    """
    if wcs is not None:
        skyPatchList = [wcs.pixelToSky(pos).getPosition(afw_geom.degrees)
                        for pos in vertexList]
    else:
        skyPatchList = vertexList
    verts = [(coord[0], coord[1]) for coord in skyPatchList]
    verts.append((0, 0))
    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY,
             ]
    return Path(verts, codes)

def plot_sensors(sensors, obs_md, ax=None, color='red', figsize=(8, 6)):
    """
    Plot the CCDs in the LSST focal plane using CCD coordinates
    derived from the pointing info using the lsst.sims code.

    Parameters
    ----------
    sensors: list
        List of sensors to plot.
    obs_md: ObservationMetaData
        Instance of an lsst_sims observation metadata class.
    color: str ['red']
        Color to use for plotting the individual CCDs.
    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot: The subplot object used for plotting.
    """
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

    camera = desc.imsim.get_obs_lsstSim_camera()

    # Re-order the CCD vertex list returned by the lsst_sims code so
    # that a rectangle is plotted.
    corner_index = (np.array([0, 1, 3, 2]),)
    for detname in sensors:
        corners = np.array(lsst.sims.coordUtils.getCornerRaDec(detname, camera,
                                                               obs_md))
        path = make_patch(corners[corner_index])
        ccd = patches.PathPatch(path, alpha=0.2, lw=1, color=color)
        ax.add_patch(ccd)

    return ax

def plot_survey_region(ax=None):
    """Plot the Run2.0 WFD region."""
    ra = (73.79, 71.46, 52.25, 49.92, 73.79)
    dec = (-44.33, -27.25, -27.25, -44.33, -44.33)
    if ax is None:
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111)
    plt.errorbar(ra, dec, fmt='-', label='DC2 WFD region')
    return ax

def set_xylims(instcat, radius=2.1):
    """Set the plot limits so that it encompasses the focalplane."""
    with open(instcat, 'r') as fd:
        for line in fd:
            if line.startswith('rightascension'):
                ra = float(line.strip().split()[1])
            if line.startswith('declination'):
                dec = float(line.strip().split()[1])
    plt.ylim(dec - radius, dec + radius)
    dra = radius/np.cos(np.radians(dec))
    plt.xlim(ra - dra, ra + dra)

if __name__ == '__main__':
    from trim_sensors import Run20Region
    ax = plot_survey_region()
    run20_region = Run20Region()
    for instcat in ('phosim_cat_3040.txt', 'phosim_cat_2188.txt'):
        obs_md, _, _ \
            = desc.imsim.parsePhoSimInstanceFile(instcat, (), numRows=50)
        sensors = run20_region.trim_sensors(instcat)
        ax = plot_sensors(sensors, obs_md, ax=ax)
