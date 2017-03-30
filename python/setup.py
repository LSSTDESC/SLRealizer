from setuptools import setup

setup(# package information
      name="slrealizer",
      version="0.0",
      author="Phil Marshall",
      author_email="dr.phil.marshall@gmail.com",
      description="Catalog-level realization of simulated gravitational lenses",
      long_description=open("README.md").read(),
      url="https://github.com/LSSTDESC/SLRealizer",
      #packages=['desc.slrealizer'],
      packages=['slrealizer'],
      #package_dir={'desc.slrealizer':'python/desc'},
      package_dir={'slrealizer':'desc'},
      include_package_data=True,
      package_data={'slrealizer': ['../data/*']},
      classifiers=[
          "Development Status :: 4 - Beta",
          "License :: OSI Approved :: BSD License",
          "Intended Audience :: Developers",
          "Intended Audience :: Science/Research",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
      ],
      install_requires=["numpy", "scipy", "matplotlib", "astropy", "om10"],
)
