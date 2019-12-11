"""
Estimate the depth for skyMap patches for given raw files.
"""
import os
import sys
import pickle
import warnings
from collections import defaultdict
import sqlite3
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import lsst.geom as lsst_geom
from lsst.afw import cameraGeom
import lsst.obs.lsst as obs_lsst
import lsst.sphgeom
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from lsst.sims.coordUtils import getCornerRaDec
    from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
    from lsst.sims.utils import ObservationMetaData, getRotSkyPos

__all__ = ['SkyMapDepth', 'SkyMapPolygons', 'make_box_wcs_region',
           'DescOpsimDbParams', 'get_det_polygon', 'process_registry_file']


class DescOpsimDbParams:
    def __init__(self, opsim_db_file):
        self.conn = sqlite3.connect(opsim_db_file)
    def update_obs_md(self, obs_md, visit):
        """
        Update the obs_md object with the DESC dithered versions
        of pointing parameters.
        """
        query = '''select descDitheredRA, descDitheredDec,
                   descDitheredRotTelPos from summary
                   where obshistid={}'''.format(visit)
        curs = self.conn.execute(query)
        ra, dec, rottelpos = [np.degrees(_) for _ in curs][0]
        obs_md.pointingRA = ra
        obs_md.pointingDec = dec
        obs_md.rotSkyPos = getRotSkyPos(ra, dec, obs_md, rottelpos)
        return obs_md

class SkyMapDepth:
    """
    Class to estimate depth in skyMap patches for a set of sensor-visits.
    """
    def __init__(self, sky_map_file, tract_info_file='tract_info.pkl',
                 camera=None):
        """
        Parameters
        ----------
        sky_map_file: str
            Pickle file containing the skyMap.
        tract_info_file: str ['tract_info.pkl']
            Pickle file containing the convex polygons for each tract.
            If it does not exist, then the convex polygons are computed
            from the skyMap data and then saved as this file.
        camera: lsst.afw.cameraGeom.camera.Camera [None]
            Camera object that provides sensor locations in focalplane.
            If None, then lsst.obs.lsst.lsstCam.LsstCam is used.
        """
        with open(sky_map_file, 'rb') as fd:
            self.sky_map = pickle.load(fd)
        self.sky_map_polygons \
            = SkyMapPolygons(self.sky_map, tract_info_file=tract_info_file)
        self.camera \
            = obs_lsst.LsstCamMapper().camera if camera is None else camera
        self.obs_mds = dict()
        self.detectors = {det.getName(): det for det in self.camera
                          if det.getType() == cameraGeom.SCIENCE}
        columns = 'band tract patch visit detname'.split()
        self.df = pd.DataFrame(columns=columns)
        self.ra, self.dec, self.visits = None, None, None

    def process_registry_file(self, registry_file, opsim_db, constraint=None,
                              row_bounds=None):
        """
        Process a data repo registry file, applying an optional
        constraint on the query for entries.

        Parameters
        ----------
        registry_file: str
            registry.sqlite3 file for a DM data repo.
        opsim_db: str
            Opsim database sqlite3 file.  Must be modified with DESC dithered
            pointing info.
        sky_map_polygons: desc.sim_utils.SkyMapPolygons
            Collection of convex polygons corresponding to the tracts and
            patches in the associated sky map.
        constraint: str [None]
            sqlite3 constraint on the selection of rows from the raw table
            in registry.sqlite3 file.
        row_bounds: (int, int) [None]
            Minimum and maximum row to process from the query of the raw
            table in the registry file. If None, then process all rows.
        """
        self.df = process_registry_file(registry_file, opsim_db,
                                        self.sky_map_polygons,
                                        constraint=constraint,
                                        row_bounds=row_bounds,
                                        camera=self.camera)

    def process_raw_files(self, raw_files):
        """
        Process a set of raw files, adding band, tract, patch, visit,
        detname info to the internal data frame containing the
        sensor-visit to patch mappings.

        Parameters
        ----------
        raw_files: sequence
            The raw files to process.  The filenames are assumed to follow
            the imSim naming convention, i.e.,
            `lsst_a_<visit>_<raftName>_<detectorName>_<band>.fits`.
        """
        rows = []
        for item in raw_files:
            # Assuming raw filenames follow the imSim naming convention.
            tokens = os.path.basename(item).split('_')
            visit = int(tokens[2])
            band = tokens[-1][0]
            if visit not in self.obs_mds:
                self.obs_mds[visit] = get_obs_md_from_raw(item)
            detname = '_'.join(tokens[3:5])
            det = self.detectors[detname]
            polygon = get_det_polygon(det, self.camera, self.obs_mds[visit])
            tract_info = self.sky_map_polygons.find_overlaps(polygon)
            for tract, patches in tract_info:
                for patch in patches:
                    rows.append((band, tract, '%s,%s' % patch, visit, detname))
        self.df = self.df.append(pd.DataFrame(rows, columns=self.df.columns),
                                 ignore_index=True)

    def compute_patch_visits(self):
        """
        Compute the number of visits per patch for each band.
        """
        self.ra = defaultdict(list)
        self.dec = defaultdict(list)
        self.visits = defaultdict(list)
        for band in set(self.df['band']):
            df_band = self.df.query('band=="%s"' % band)
            tracts = set(df_band['tract'])
            for tract in tracts:
                df_tract = df_band.query('tract==%d' % tract)
                patches = set(df_tract['patch'])
                for patch in patches:
                    ra, dec = get_patch_center(self.sky_map[tract], eval(patch))
                    self.ra[band].append(ra)
                    self.dec[band].append(dec)
                    visits = set(df_tract.query('patch=="%s"' % patch)['visit'])
                    self.visits[band].append(len(visits))

    def plot_visit_depths(self, band, title=None, ax=None, figsize=None,
                          markersize=None, cmap_name='plasma'):
        """
        Plot number of visits per patch.

        Parameters
        ----------
        band: str
            The ugrizy band to be plotted.
        title: str [None]
            The title of the plot. If None, then it will be set to
            `'# visits per patch, %s-band' % band`.
        ax: matplotlib.axes._subplots.AxesSubplot [None]
            The subplot object to contain the tract plot.  If None, then
            make a new one.
        figsize: (float, float) [None]
            Size of the figure in inches.
        markersize: float [None]
            If None, then the default is used: `rcParams['lines.markersize']**2`
        cmap_name: str ['plasma']
            The name of the colormap to use.

        Returns
        -------
        matplotlib.axes._subplots.AxesSubplot
        """
        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        if title is None:
            title = '# visits per patch, %s-band' % band

        cmap = matplotlib.cm.get_cmap(cmap_name)
        if self.ra is None:
            self.compute_patch_visits()
        plt.scatter(self.ra[band], self.dec[band], c=self.visits[band],
                    s=markersize, cmap=cmap)
        plt.colorbar()
        plt.title(title)
        plt.xlabel('RA (degrees)')
        plt.ylabel('Dec (degrees)')
        return ax


def get_patch_center(tract_info, patch_id):
    """.

    Compute the coordinates of the center of the patch.

    Parameters
    ----------
    tract_info: lsst.skymap.tractInfo.TractInfo
        The tract info object for the desired tract.
    patch_id: (int, int)
        The patch id, e.g., (1, 1).

    Returns
    -------
    (float, float): (RA, Dec) in degrees.
    """
    patch_info = tract_info.getPatchInfo(patch_id)
    patch_box = lsst_geom.Box2D(patch_info.getOuterBBox())
    return tract_info.getWcs().pixelToSky(patch_box.getCenter())\
                              .getPosition(lsst_geom.degrees)


def get_obs_md_from_raw(raw_file):
    """
    Get the ObservationMetaData object from the raw file header info.

    Parameters
    ----------
    raw_file: str
        The path to the raw file.

    Results
    -------
    lsst.sims.utils.ObservationMetaData
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

    Parameters
    ----------
    detector: lsst.afw.cameraGeom.detector.detector.Detector
        The detector at issue.
    camera: lsst.afw.cameraGeom.camera.Camera [None]
        Camera object that provides sensor locations in focalplane.
        If None, then lsst.obs.lsst.lsstCam.LsstCam is used.
    obs_md: lsst.sims.utils.ObservationMetaData
        The object containing the visit info.

    Returns
    -------
    lsst.sphgeom.ConvexPolygon
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
        polygon: lsst.sphgeom.ConvexPolygon
            The polygon representing the bbox at issue.

        Returns
        -------
        list of tuples of (tract, list of patches) that overlap the box.
        """
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
    box: lsst.geom.Box2I or lsst.geom.Box2D
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
    box = lsst_geom.Box2D(box)
    box.grow(margin)
    vertices = []
    for point in box.getCorners():
        coord = wcs.pixelToSky(point)
        lonlat = lsst.sphgeom.LonLat.fromRadians(coord.getRa().asRadians(),
                                                 coord.getDec().asRadians())
        vertices.append(lsst.sphgeom.UnitVector3d(lonlat))
    return lsst.sphgeom.ConvexPolygon(vertices)


def process_registry_file(registry_file, opsim_db, sky_map_polygons,
                          constraint=None, row_bounds=None, camera=None):
    """
    Process a data repo registry file to obtain the sensor visit overlaps
    sky map tracts and patches.

    Parameters
    ----------
    registry_file: str
        registry.sqlite3 file for a DM data repo.
    opsim_db: str
        Opsim database sqlite3 file.  Must be modified with DESC dithered
        pointing info.
    sky_map_polygons: desc.sim_utils.SkyMapPolygons
        Collection of convex polygons corresponding to the tracts and
        patches in the associated sky map.
    constraint: str [None]
        sqlite3 constraint on the selection of rows from the raw table
        in registry.sqlite3 file.
    row_bounds: (int, int) [None]
        Minimum and maximum row to process from the query of the raw table
        in the registry file. If None, then process all rows.
    camera: lsst.afw.cameraGeom.Camera [None]
        Camera object. If None, then use `lsst.obs.lsst.LsstMapper().camera`.

    Returns
    -------
    pandas DataFrame
    """
    if camera is None:
        camera = obs_lsst.LsstCamMapper().camera
    obs_gen = ObservationMetaDataGenerator(database=opsim_db, driver='sqlite')
    desc_opsim_db_params = DescOpsimDbParams(opsim_db)
    rows = []
    query = 'select visit, filter, raftName, detectorName from raw'
    if constraint is not None:
        query += ' where {}'.format(constraint)
    with sqlite3.connect(registry_file) as conn:
        registry = pd.read_sql(query, conn)
    nrows = len(registry)
    print("processing {} rows from {}".format(nrows, registry_file))
    sys.stdout.flush()
    obs_mds = dict()
    detectors = {det.getName(): det for det in camera
                 if det.getType() == cameraGeom.SCIENCE}

    df = pd.DataFrame(columns='band tract patch visit detname'.split())
    if row_bounds is not None:
        rmin = row_bounds[0]
        rmax = min(row_bounds[1], nrows)
    for iloc in range(rmin, rmax):
        print('{} / {}'.format(iloc, rmax - rmin))
        sys.stdout.flush()
        row = registry.iloc[iloc]
        visit = row['visit']
        band = row['filter']
        if visit not in obs_mds:
            obs_md = obs_gen.getObservationMetaData(obsHistID=visit,
                                                    boundType='circle',
                                                    boundLength=0)[0]
            obs_mds[visit] = desc_opsim_db_params.update_obs_md(obs_md, visit)
        detname = '_'.join((row['raftName'], row['detectorName']))
        det = detectors[detname]
        polygon = get_det_polygon(det,camera, obs_mds[visit])
        tract_info = sky_map_polygons.find_overlaps(polygon)
        for tract, patches in tract_info:
            for patch in patches:
                rows.append((band, tract, '%s,%s' % patch, visit, detname))
    df = df.append(pd.DataFrame(rows, columns=df.columns), ignore_index=True)
    return df
