"""
Functions to plot LSST sensors on the sky.
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

__all__ = ['make_patch', 'plot_sensors', 'plot_skymap_tract', 'set_xylims']


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


def plot_sensors(sensors, camera, obs_md, ax=None, color='red',
                 figsize=(8, 6), label_sensors=False):
    """
    Plot the CCDs in the LSST focal plane using CCD coordinates
    derived from the pointing info using the lsst.sims code.

    Parameters
    ----------
    sensors: list
        List of sensors to plot.
    camera: lsst.afw.cameraGeom.Camera
        Camera object that contains the camera info.
    obs_md: ObservationMetaData
        Object containing the pointing information for the
        desired visit.
    ax: matplotlib.axes.Axes
        matplotlib container class of figure elements.
    color: str ['red']
        Color to use for plotting the individual CCDs.
    figsize: tuple [(8, 6)]
        Size of the figure in inches
    label_sensors: bool [False]
        If True, then label the sensors on the map by their det_name.

    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot: The subplot object used for plotting.
    """
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

    # Re-order the CCD vertex list returned by the lsst_sims code so
    # that a rectangle is plotted.
    corner_index = (np.array([0, 1, 3, 2]),)
    for detname in sensors:
        corners = np.array(lsst.sims.coordUtils.getCornerRaDec(detname, camera,
                                                               obs_md))
        path = make_patch(corners[corner_index])
        ccd = patches.PathPatch(path, alpha=0.2, lw=1, color=color)
        ax.add_patch(ccd)
        if label_sensors:
            ras = [_[0] for _ in corners]
            decs = [_[1] for _ in corners]
            center = ((min(ras) + max(ras))/2., (min(decs) + max(decs))/2.)
            ax.text(center[0], center[1], detname, size=6,
                    ha="center", va="center")

    return ax


def plot_skymap_tract(skyMap, tract=0, title=None, ax=None,
                      patch_colors=None):
    """
    Plot a tract from a skyMap.

    Parameters
    ----------
    skyMap: lsst.skyMap.SkyMap
        The SkyMap object containing the tract and patch information.
    tract: int [0]
        The tract id of the desired tract to plot.
    title: str [None]
        Title of the tract plot.  If None, the use `tract <id>`.
    ax: matplotlib.axes._subplots.AxesSubplot [None]
        The subplot object to contain the tract plot.  If None, then
        make a new one.
    patch_colors: dict [None]
        Dictionary of colors keyed by patchId.

    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot: The subplot containing the
    tract plot.
    """
    if title is None:
        title = 'tract {}'.format(tract)
    tract_info = skyMap[tract]
    tractBox = afw_geom.Box2D(tract_info.getBBox())
    tractPosList = tractBox.getCorners()
    wcs = tract_info.getWcs()
    xNum, yNum = tract_info.getNumPatches()

    if ax is None:
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111)

    tract_center = wcs.pixelToSky(tractBox.getCenter())\
                      .getPosition(afw_geom.degrees)
    ax.text(tract_center[0], tract_center[1], '%d' % tract, size=16,
            ha="center", va="center", color='blue')
    for x in range(xNum):
        for y in range(yNum):
            patch_info = tract_info.getPatchInfo([x, y])
            patchBox = afw_geom.Box2D(patch_info.getOuterBBox())
            pixelPatchList = patchBox.getCorners()
            path = make_patch(pixelPatchList, wcs)
            try:
                color = patch_colors[(x, y)]
            except (TypeError, KeyError):
                color = 'blue'
            patch = patches.PathPatch(path, alpha=0.1, lw=1, color=color)
            ax.add_patch(patch)
            center = wcs.pixelToSky(patchBox.getCenter())\
                        .getPosition(afw_geom.degrees)
            ax.text(center[0], center[1], '%d,%d' % (x, y), size=6,
                    ha="center", va="center")

    skyPosList = [wcs.pixelToSky(pos).getPosition(afw_geom.degrees)
                  for pos in tractPosList]
    ax.set_xlim(max(coord[0] for coord in skyPosList) + 1,
                min(coord[0] for coord in skyPosList) - 1)
    ax.set_ylim(min(coord[1] for coord in skyPosList) - 1,
                max(coord[1] for coord in skyPosList) + 1)
    ax.grid(ls=':', color='gray')
    ax.set_xlabel("RA (deg.)")
    ax.set_ylabel("Dec (deg.)")
    ax.set_title(title)
    return ax


def set_xylims(instcat, radius=2.1):
    """Set the plot limits so that it encompasses the focalplane."""
    with open(instcat, mode='r') as fd:
        for _, line in zip(range(30), fd):
            if line.startswith('rightascension'):
                ra = float(line.strip().split()[1])
            if line.startswith('declination'):
                dec = float(line.strip().split()[1])
    plt.ylim(dec - radius, dec + radius)
    dra = radius/np.cos(np.radians(dec))
    plt.xlim(ra - dra, ra + dra)
