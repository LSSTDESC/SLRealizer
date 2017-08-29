# SLRealizer

![Logo](https://i.imgur.com/HAHxXNY.png)

<a href='https://travis-ci.org/jennykim1016/SLRealizer'>
<img src='https://secure.travis-ci.org/jennykim1016/SLRealizer.png?branch=master'></a>

Catalog-level simulation of LSST DM stack measurements of
gravitationally-lensed quasars.

In this package we are experimenting with predicting `Source`, `Object`,
`DIASource` and `DIAObject` properties directly, by realizing
model lenses as mixtures of Gaussians. Our goals are to

1. Be able to generate large-volume catalog level lens simulations, to enable catalog-level lens finder training;

2. Enable a model-based interpretation of the LSST catalogs, for use in lens finding and light curve extraction.


## Set-up and testing

To use (and not develop) `SLRealizer`, you can pip install it with
```
pip install git+git://github.com/LSSTDESC/SLRealizer.git#egg=slrealizer
```

To help develop it, fork and clone the repo, `cd` to its `python` folder, and then install the package with
```
python setup.py develop
```
To run the tests, do
```
nosetests
```

## Demo

* [Script to generate the source and the object catalog](https://github.com/jennykim1016/SLRealizer/blob/master/python/desc/slrealizer/script.py)
* [Null-deblending demo with SLRealizer](https://github.com/jennykim1016/SLRealizer/blob/master/notebooks/Null_Deblend_Demo.ipynb)
* [Feature Comparison between SDSS and OM10](https://github.com/jennykim1016/SLRealizer/blob/master/notebooks/SDSSvsOM10.ipynb)
* [Classification of lensed systems](https://github.com/jennykim1016/SLRealizer/blob/master/notebooks/Scikit_Classifying_Lensed_Systems.ipynb)


## Repository Structure

```
SLRealizer/
├── python/
│   ├── desc/
|   |    └── slrealizer/
|   |    |   ├── __init__.py
|   |    |   ├── SLRealizer.py // contains wrapper methods that ensembles the methods in different modules
|   |    |   ├── null_deblend.py // contains methods that null-deblends the input sources
|   |    |   ├── catalog.py // contains methods that format each entree of the source and the object catalog
|   |    |   ├── constant.py // defines constants that could be accessed globally in the package
|   |    |   ├── dropbox_manage.py // contains methods that sync dropbox repository with a local repository to save big data
|   |    |   ├── test.py // contains unit tests for the SLRealizer package
|   |    |   ├── plotting.py // _not called by SLRealizer.py_ given an OM10 system, make a simple plot of gaussian realization
|   |    |   ├── distance.py // _not yet used by SLRealizer.py_ calculate the image distance and evaluate the cost model
|   |    |   └── Makefile // _to be more updated_ currently, contains `make clean` that cleans up the junk files
│   └── setup.py
└── data/
    ├── qso_mock.fits // OM10 catalog
    ├── sdss_object.csv // SDSS object catalog, used to compare the features between OM10 and sdss
    └── twinkles_observation_history.csv // mock observation data
```

## People
* [Phil Marshall](https://github.com/LSSTDESC/SLRealizer/issues/new?body=@drphilmarshall) (SLAC)
* [Jenny Kim](https://github.com/LSSTDESC/SLRealizer/issues/new?body=@jennykim1016) (Stanford)
* [Mike Baumer](https://github.com/LSSTDESC/SLRealizer/issues/new?body=@mbaumer) (SLAC)
* [Steve Kahn](https://github.com/LSSTDESC/SLRealizer/issues/new?body=@stevkahn) (SLAC)
* [Rahul Biswas](https://github.com/LSSTDESC/SLRealizer/issues/new?body=@rbiswas4) (UW)


## License, etc.

This is research in progress, using open source software which is available for
re-use under the BSD license. If you make use of any of the ideas or code in
this repository in your own research, please do cite us as "(Marshall et al, in
prep.)" and provide a link to this repository - but before then, we hope you
will get in touch! Please do drop us a line via the hyperlinked contact names
above, or by [writing us an
issue](https://github.com/LSSTDESC/SLRealizer/issues/new).
