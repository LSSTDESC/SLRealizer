"""
script that runs to generate toy catalogs
"""

import lenspop
import om10
import desc.slrealizer
import format_sdss_to_om10

# initialize om10
db = om10.DB(vb=False, catalog='../../../data/qso_mock.fits')
# we select LSST-like lenses by quering with these values
db.select_random(maglim=23.3,area=20000.0,IQ=0.75)
# calculate the synthetic magnitude for the OM10 catalog -- need this
db.paint(synthetic=True)
# calculate the sizes for the OM10 catalog -- if you want to assume the point-source-like-galaxies, get rid of this
db.calculate_size()
# initialize realizer
realizer = desc.slrealizer.SLRealizer(catalog=db, observation="../../../data/twinkles_observation_history.csv")
realizer.make_source_catalog(dir='../../../data/source_catalog_galsim_noise_perfect.csv')
realizer.make_object_catalog(source_table_dir='../../../data/source_catalog_galsim_noise_perfect.csv', save_dir='../../../data/object_catalog_galsim_noise_perfect.csv')
