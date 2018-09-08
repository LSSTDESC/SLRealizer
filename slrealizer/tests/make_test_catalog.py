import pyfits
import os

catalog_f = os.path.join(os.environ['SLREALIZERDIR'], 'data', 'qso_mock.fits')
catalog = pyfits.open(catalog_f)[1].data

mask = catalog['LENSID'] == 6136045
test_catalog = catalog[mask]
hdu = pyfits.BinTableHDU(data=test_catalog)
hdu.writeto('test_catalog.fits')
