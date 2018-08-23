#from distutils.core import setup
from setuptools import setup

setup(
    name='slrealizer',
    version='0.2.1',
    author='Phil Marshall',
    author_email='dr.phil.marshall@gmail.com',
    packages=['slrealizer'],
    license='LICENSE.md',
    description='Catalog-level realization of simulated gravitational lenses.',
    long_description=open("README.md").read(),
    long_description_content_type='text/markdown',
    url='https://github.com/LSSTDESC/SLRealizer',
    package_dir={'slrealizer':'slrealizer'},
    include_package_data=True,
    package_data={'slrealizer': ['/data/*']},
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
    ],
    install_requires=[
        "matplotlib >= 2.2.2",
        "numpy >= 1.13.1",
        "astropy >= 0.18.1",
        "om10",
    ],
)
