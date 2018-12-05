"""
Code to repackage phosim amplifier files into single sensor
multi-extension FITS files, with an HDU per amplifier.
"""
import os
import sys
import glob
import time
import warnings
from collections import defaultdict
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning, AstropyUserWarning
import astropy.time
import lsst.obs.lsst as obs_lsst

__all__ = ['PhoSimRepackager', 'noao_section_keyword']


def noao_section_keyword(bbox, flipx=False, flipy=False):
    """
    Convert bounding boxes into NOAO section keywords.

    Parameters
    ----------
    bbox : lsst.afw.geom.Box2I
        Bounding box.
    flipx : bool
        Flag to indicate that data should be flipped in the x-direction.
    flipy : bool
        Flag to indicate that data should be flipped in the y-direction.
    """
    xmin, xmax = bbox.getMinX()+1, bbox.getMaxX()+1
    ymin, ymax = bbox.getMinY()+1, bbox.getMaxY()+1
    if flipx:
        xmin, xmax = xmax, xmin
    if flipy:
        ymin, ymax = ymax, ymin
    return '[%i:%i,%i:%i]' % (xmin, xmax, ymin, ymax)


class PhoSimRepackager:
    """
    Class to repackage phosim amplifier files into single sensor
    MEFs with one HDU per amp.
    """
    def __init__(self):
        self.amp_info_records = list(list(obs_lsst.LsstCamMapper().camera)[0])

    def process_visit(self, visit_dir, out_dir=None, image_type='SKYEXP',
                      verbose=False):
        """
        Parameters
        ----------
        visit_dir: str
            Directory containing the phosim amplifier for a given visit.
        out_dir: str [None]
            Output directory for MEF files. If None, then a directory
            with name v<visit #>-<band> will be created in the cwd.
        image_type: str ['SKYEXP']
            Image type, e.g., 'SKYEXP', 'BIAS', 'DARK', 'FLAT'.
        verbose: bool [False]
            Set to True to print out time for processing each sensor.
        """
        phosim_amp_files \
            = sorted(glob.glob(os.path.join(visit_dir, 'lsst_a_*')))
        amp_files = defaultdict(list)
        for item in phosim_amp_files:
            sensor_id = '_'.join(os.path.basename(item).split('_')[4:6])
            amp_files[sensor_id].append(item)
        if out_dir is None:
            tokens = os.path.basename(phosim_amp_files[0]).split('_')
            out_dir = 'v%07i-%s' % (int(tokens[2]), 'ugrizy'[int(tokens[3][1])])
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        for sensor_id in amp_files:
            if verbose:
                sys.stdout.write(sensor_id + '  ')
            t0 = time.time()
            # Skip already-processed and compressed files.
            gzip_file = self.mef_filename(amp_files[sensor_id][0],
                                          out_dir=out_dir) + '.gz'
            if os.path.isfile(gzip_file):
                continue
            print("repackaging", gzip_file)
            self.repackage(amp_files[sensor_id], out_dir=out_dir,
                           image_type=image_type)
            if verbose:
                print(time.time() - t0)
                sys.stdout.flush()

    def repackage(self, phosim_amp_files, out_dir='.', image_type='SKYEXP'):
        """
        Repackage a collection of phosim amplifier files for a
        single sensor into a multi-extension FITS file.

        Parameters
        ----------
        phosim_amp_files: list
            List of phosim amplifier filenames.
        """
        # PhoSim labels the amplifier channels incorrectly.  Here is a
        # mapping from the phosim channels to the correct ones.
        ch_map = {'00': '10',
                  '01': '11',
                  '02': '12',
                  '03': '13',
                  '04': '14',
                  '05': '15',
                  '06': '16',
                  '07': '17',
                  '17': '07',
                  '16': '06',
                  '15': '05',
                  '14': '04',
                  '13': '03',
                  '12': '02',
                  '11': '01',
                  '10': '00'}

        # Create the HDUList to contain the MEF data.
        phdu = fits.PrimaryHDU()
        sensor = fits.HDUList(phdu)

        # Extract the data for each segment from the FITS files
        # into a dictionary keyed by channel id.
        segments = dict()
        for fn in phosim_amp_files:
            # Get the channel id from the filename written by phosim.
            # This is the incorrect channel id, so correct it when
            # filling the segments dictionary.
            phosim_channel = os.path.basename(fn).split('_')[6][1:]
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=AstropyWarning,
                                        append=True)
                warnings.filterwarnings('ignore', category=AstropyUserWarning,
                                        append=True)
                with fits.open(fn) as phosim_amp:
                    # Use tile compression for the pixel data.
                    segments[ch_map[phosim_channel]] \
                        = fits.CompImageHDU(data=phosim_amp[0].data,
                                            compression_type='RICE_1')
                    # Copy the header keywords from the original file.
                    segments[ch_map[phosim_channel]]\
                        .header.update(phosim_amp[0].header)

        # Set the NOAO section keywords based on the pixel geometry
        # in the LsstCam object.
        for amp in self.amp_info_records:
            hdu = segments[amp.get('name')[1:]]
            hdu.header['EXTNAME'] = 'Segment%s' % amp.get('name')[1:]
            hdu.header['DATASEC'] = noao_section_keyword(amp.getRawDataBBox())
            hdu.header['DETSEC'] \
                = noao_section_keyword(amp.getBBox(),
                                       flipx=amp.get('raw_flip_x'),
                                       flipy=amp.get('raw_flip_y'))
            # Remove the incorrect BIASSEC keyword that phosim writes.
            try:
                hdu.header.remove('BIASSEC')
            except KeyError:
                pass

            sensor.append(hdu)

        # Set keywords in primary HDU, extracting most of the relevant
        # ones from the first phosim amplifier file.
        amp_hdr = sensor[1].header
        chip_id = amp_hdr['CHIPID']
        raft, ccd = chip_id.split('_')
        phdu.header['EXPTIME'] = amp_hdr['EXPTIME']
        phdu.header['DARKTIME'] = amp_hdr['DARKTIME']
        phdu.header['RUNNUM'] = amp_hdr['OBSID']
        phdu.header['MJD-OBS'] = amp_hdr['MJD-OBS']
        phdu.header['DATE-OBS'] \
            = astropy.time.Time(amp_hdr['MJD-OBS'], format='mjd').isot
        phdu.header['FILTER'] = amp_hdr['FILTER']
        phdu.header['LSST_NUM'] = chip_id
        phdu.header['CHIPID'] = chip_id
        phdu.header['OBSID'] = amp_hdr['OBSID']
        phdu.header['AIRMASS'] = amp_hdr['AIRMASS']
        phdu.header['TESTTYPE'] = 'PHOSIM'
        phdu.header['IMGTYPE'] = image_type
        phdu.header['MONOWL'] = -1
        phdu.header['RAFTNAME'] = raft
        phdu.header['SENSNAME'] = ccd
        # Add boresight pointing angles and rotskypos (angle of sky
        # relative to Camera coordinates) from which obs_lsst can
        # infer the CCD-wide WCS.
        phdu.header['RATEL'] = amp_hdr['RA_DEG']
        phdu.header['DECTEL'] = amp_hdr['DEC_DEG']
        phdu.header['ROTANGLE'] = amp_hdr['ROTANGZ']

        outfile = self.mef_filename(phosim_amp_files[0], out_dir=out_dir)
        sensor.writeto(outfile, overwrite=True)

    @staticmethod
    def mef_filename(phosim_amp_file, out_dir='.'):
        """
        Construct the filename of the output MEF file based
        on a phosim single amplifier filename.
        """
        tokens = os.path.basename(phosim_amp_file).split('_')
        # Remove the channel identifier.
        outfile = '_'.join(tokens[:6] + tokens[7:])
        outfile = os.path.join(out_dir, outfile)
        # Since these files use tile-compression, remove any .gz
        # extension from the computed output filename.
        if outfile.endswith('.gz'):
            outfile = outfile[:-len('.gz')]
        return outfile
