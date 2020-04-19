import unittest
import lsst.sims.coordUtils
import lsst.sims.utils
import desc.sim_utils

class TrimSensorsTestCase(unittest.TestCase):
    """Test case class for trim_sensors function."""
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_trim_sensors(self):
        """Function to test trim_sensors."""
        camera = lsst.sims.coordUtils.lsst_camera()
        obs_md = lsst.sims.utils.ObservationMetaData(
            mjd=60624.111854,
            pointingRA=52.93798307449681,
            pointingDec=-27.466575126110218,
            rotSkyPos=149.45133356989606,
            bandpassName='z')
        ddf = desc.sim_utils.Run20Region(ra_mid=53.125,
                                         ne_corner=(53.764, -27.533),
                                         dec_range=(-28.667, -27.533))
        sensors = set(ddf.trim_sensors(obs_md))
        self.assertEqual(len(sensors), 34)
        sensors_to_exclude = set(['R:0,3 S:1,1', 'R:1,1 S:1,2',
                                  'R:2,4 S:0,1', 'R:3,2 S:0,2'])
        self.assertEqual(len(sensors.intersection(sensors_to_exclude)), 0)


if __name__ == '__main__':
    unittest.main()
