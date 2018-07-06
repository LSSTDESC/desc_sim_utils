"""
Example unit tests for desc_sim_utils package
"""
import unittest
import lsst.afw.geom as afw_geom
import desc.sim_utils

class PhosimRepackagerTestCase(unittest.TestCase):
    """TestCase class for phosim_repackager."""
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_phosim_repackager(self):
        """Unit tests for phosim_repackager"""
        # Test NOAO section keyword generator.
        xmin, xmax, ymin, ymax = 3, 512, 0, 2000
        bbox = afw_geom.Box2I(afw_geom.Point2I(xmin, ymin),
                              afw_geom.Point2I(xmax, ymax))
        value = desc.sim_utils.noao_section_keyword(bbox)
        self.assertEqual('[4:513,1:2001]', value)
        value = desc.sim_utils.noao_section_keyword(bbox, flipx=True)
        self.assertEqual('[513:4,1:2001]', value)
        value = desc.sim_utils.noao_section_keyword(bbox, flipy=True)
        self.assertEqual('[4:513,2001:1]', value)

        # Test multi-extension filename generator.
        phosim_amp_file = 'lsst_a_219976_f2_R22_S11_C00_E000.fits.gz'
        value = desc.sim_utils.PhoSimRepackager.mef_filename(phosim_amp_file)
        self.assertEqual('./lsst_a_219976_f2_R22_S11_E000.fits', value)

if __name__ == '__main__':
    unittest.main()
