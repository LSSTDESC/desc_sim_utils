"""
Estimate the depth for skyMap patches for given raw files.
"""
import os
import pickle
import warnings
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import lsst.afw.geom as afw_geom
from lsst.afw import cameraGeom
import lsst.obs.lsst as obs_lsst
import lsst.sphgeom
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from lsst.sims.coordUtils import getCornerRaDec
    from lsst.sims.utils import ObservationMetaData

__all__ = ['SkyMapDepth', 'SkyMapPolygons', 'make_box_wcs_region']


class SkyMapDepth:
    """
    Class to estimate depth in skyMap patches for a set of sensor-visits.
    """
    def __init__(self, sky_map_file, tract_info_file='tract_info.pkl',
                 camera=None):
        with open(sky_map_file, 'rb') as fd:
            self.sky_map = pickle.load(fd)
        self.sky_map_polygons \
            = SkyMapPolygons(self.sky_map, tract_info_file=tract_info_file)
        self.camera \
            = obs_lsst.LsstCamMapper().camera if camera is None else camera
        self.obs_mds = dict()
        self.detectors = {det.getName(): det for det in self.camera
                          if det.getType() == cameraGeom.SCIENCE}
        columns = 'tract patch visit detname'.split()
        self.df = pd.DataFrame(columns=columns)

    def process_raw_files(self, raw_files):
        """
        Process a set of raw files.
        """
        rows = []
        print('processing:')
        for item in raw_files:
            print(item)
            # Assuming raw filenames follow the imSim naming convention.
            tokens = os.path.basename(item).split('_')
            visit = int(tokens[2])
            if visit not in self.obs_mds:
                self.obs_mds[visit] = get_obs_md_from_raw(item)
            detname = '_'.join(tokens[3:5])
            det = self.detectors[detname]
            polygon = get_det_polygon(det, self.camera, self.obs_mds[visit])
            tract_info = self.sky_map_polygons.find_overlaps(polygon)
            for tract, patches in tract_info:
                for patch in patches:
                    rows.append((tract, '%s,%s' % patch, visit, detname))
        self.df = self.df.append(pd.DataFrame(rows, columns=self.df.columns),
                                 ignore_index=True)

    def plot_visit_depths(self, bins=None, xy_range=None, ax=None,
                          figsize=None):
        """
        Plot a 2D histogram of visits per patch by histogramming patch
        centers on the sky, weighted by the number of visits per patch.
        """
        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)

        tracts = set(self.df['tract'])
        ras, decs, weights = [], [], []
        for tract in tracts:
            df_tract = self.df.query('tract==%d' % tract)
            patches = set(df_tract['patch'])
            for patch in patches:
                ra, dec = get_patch_center(self.sky_map[tract], eval(patch))
                ras.append(ra)
                decs.append(dec)
                weights.append(len(df_tract.query('patch=="%s"' % patch)))
        plt.hist2d(ras, decs, weights=weights, bins=bins, range=xy_range)
        return ax


def get_patch_center(tract_info, patch_id):
    patch_info = tract_info.getPatchInfo(patch_id)
    patch_box = afw_geom.Box2D(patch_info.getOuterBBox())
    return tract_info.getWcs().pixelToSky(patch_box.getCenter())\
                              .getPosition(afw_geom.degrees)


def get_obs_md_from_raw(raw_file):
    """
    Get the ObservationMetaData object from the raw file header info.
    """
    with fits.open(raw_file) as fd:
        hdr = fd[0].header
        obs_md = ObservationMetaData(pointingRA=hdr['RATEL'],
                                     pointingDec=hdr['DECTEL'],
                                     mjd=hdr['MJD-OBS'],
                                     rotSkyPos=hdr['ROTANGLE'])
        obs_md.OpsimMetaData = dict(obshistid=hdr['OBSID'])
    return obs_md

def get_det_polygon(detector, camera, obs_md):
    """
    Compute the convex polygon for a detector for a given observation.
    """
    corners = getCornerRaDec(detector.getName(), camera, obs_md)
    vertices = []
    for corner in corners:
        lonlat = lsst.sphgeom.LonLat.fromDegrees(*corner)
        vertices.append(lsst.sphgeom.UnitVector3d(lonlat))
    return lsst.sphgeom.ConvexPolygon(vertices)

class SkyMapPolygons:
    """
    Class to compute skymap tract-patch combinations that overlap with
    a specified wcs-projected bounding box.

    Code stolen from
    ImageProcessingPipelines/python/util/tract2visit_mapper.py to find
    overlaps of wcs-projected bounding boxes with skyMap tracts and
    patches.
    """
    def __init__(self, sky_map, tract_info_file='tract_info.pkl'):
        self.sky_map = sky_map
        self.tracts = {}
        self.patches = {}
        if os.path.isfile(tract_info_file):
            print('Retrieving tract info from %s' % tract_info_file)
            with open(tract_info_file, 'rb') as handle:
                self.tracts = pickle.load(handle)
        else:
            self._compute_tract_info(tract_info_file)

    def _compute_tract_info(self, tract_info_file):
        """
        Compute the convex polygons for each tract and persist in a
        dictionary to tract_info_file.
        """
        for _, tract_info in enumerate(self.sky_map):
            if _ % 100 == 0 and _ > 0:
                print("Prepping tract %d of %d" % (_, len(self.sky_map)))
            self.tracts[tract_info.getId()] \
                = make_box_wcs_region(tract_info.getBBox(), tract_info.getWcs())
        with open(tract_info_file, 'wb') as handle:
            pickle.dump(self.tracts, handle)

    def find_overlaps(self, polygon):
        """
        Find all patch-tract combinations that overlap with the specified
        bounding box.

        Parameters
        ----------
        box: lsst.afw.geom.Box2[ID]
            The bounding box of the region.
        wcs: lsst.afw.image.wcs
            The WCS to use to convert to a convex polygon.

        Returns
        -------
        list of tuples of (tract, list of patches) that overlap the box.
        """
#        polygon = make_box_wcs_region(box, wcs, margin=margin)
        results = []
        for tract, tract_poly in self.tracts.items():
            if polygon.relate(tract_poly) == lsst.sphgeom.DISJOINT:
                continue
            self._ensure_patches(tract)
            results.append(
                (tract,
                 [patch for patch, patch_poly in self.patches[tract].items()
                  if polygon.relate(patch_poly) != lsst.sphgeom.DISJOINT])
            )
        return results

    def _ensure_patches(self, tract):
        """
        Ensure the convex polygons for patches in this tract are
        available.
        """
        if tract not in self.patches:
            patches = {}
            tract_info = self.sky_map[tract]
            for patch_info in tract_info:
                patches[patch_info.getIndex()] \
                    = make_box_wcs_region(patch_info.getOuterBBox(),
                                          tract_info.getWcs())
            self.patches[tract] = patches


def make_box_wcs_region(box, wcs, margin=0.0):
    """
    Construct a spherical ConvexPolygon from a WCS and a bounding box.

    Parameters
    ----------
    box: lsst.afw.geom.Box2I or lsst.afw.geom.Box2D
        A box in the pixel coordinate system defined by the WCS.
    wcs: lsst.afw.image.Wcs
        A mapping from a pixel coordinate system to the sky.
    margin: float [0.0]
        A buffer in pixels to grow the box by (in all directions) before
        transforming it to sky coordinates.

    Returns
    -------
    lsst.sphgeom.ConvexPolygon
    """
    box = afw_geom.Box2D(box)
    box.grow(margin)
    vertices = []
    for point in box.getCorners():
        coord = wcs.pixelToSky(point)
        lonlat = lsst.sphgeom.LonLat.fromRadians(coord.getRa().asRadians(),
                                                 coord.getDec().asRadians())
        vertices.append(lsst.sphgeom.UnitVector3d(lonlat))
    return lsst.sphgeom.ConvexPolygon(vertices)
