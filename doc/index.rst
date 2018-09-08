SLRealizer
==========

The SLRealizer package enables the realization, or emulation, of LSST "Object" and "Source" catalogs, given an input "truth" catalog. The output catalogs have (some of) the same columns as the actual data release production (DRP) catalogs will have. Galaxies and PSFs are modeled as mixtures of Gaussians, and convolution is carried out arithmetically; mock measured quantities are estimated from the moments of the resulting GMMs. A very simple Gaussian error model is used to add noise to the emulated measurements. Please find below the the API documentation for the code.

**Contents:**

.. toctree::
    :maxdepth: 2

    index


Basic Strong Lens Realization
-----------------------------

.. automodule:: slrealizer.realize_sl
    :members:
    :undoc-members:


Realizing OM10 Lensed Quasars
-----------------------------

.. automodule:: slrealizer.realize_om10
    :members:
    :undoc-members:


Realizing Non-Lenses: SDSS Galaxies
-----------------------------------

.. automodule:: slrealizer.realize_sdss
    :members:
    :undoc-members:


Utilities
---------

Accessing Standard Data Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: slrealizer.data.dataloader
    :members:
    :undoc-members:

Useful Functions
^^^^^^^^^^^^^^^^

.. automodule:: slrealizer.utils.utils
    :members:
    :undoc-members:

Useful Constants
^^^^^^^^^^^^^^^^

.. automodule:: slrealizer.utils.constants
    :members:
    :undoc-members:
