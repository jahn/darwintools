.. _darwintools:

Darwintools
===========

This Python package includes a number of helpful functions and scripts for
dealing with MITgcm and Darwin output.  You can download it from `github
<https://github.com/jahn/darwintools>`_.

The following functions are exposed at the package level:

- from module mds: :meth:`~darwintools.mds.rdmds` and
  :meth:`~darwintools.mds.wrmds`
- from module mnc: :meth:`~darwintools.mnc.rdmnc` and
  :meth:`~darwintools.mnc.mnc_files`
- from module ptracers: :meth:`~darwintools.ptracers.iolabel` and
  :meth:`~darwintools.ptracers.iolabel2num`
- from module diagnostics: :meth:`~darwintools.diagnostics.readstats`
- from module conversion: :meth:`~darwintools.conversion.pfromz`
- from module utils:
  :meth:`~darwintools.utils.gen_blanklist`,
  :meth:`~darwintools.utils.hfac`,
  :meth:`~darwintools.utils.readbin` and
  :meth:`~darwintools.utils.writebin`

The following functions require `Matplotlib <https://matplotlib.org/>`_
and are exposed in a separate module plotting_:

- from module plotutils: :meth:`~darwintools.plotutils.tilecmap`

The following modules are automatically imported: llc_, density_ (as dens) and
mds_.  The plotting module makes the cs_ and plotutils_ submodules available.

The package also includes a standalone script for joining tiled mnc files:
gluemncbig_.

For more functions, see the individual modules:

.. _mds:

mds
---

.. automodule:: darwintools.mds
    :members:

mnc
---

.. automodule:: darwintools.mnc
    :members:

diagnostics
-----------

.. automodule:: darwintools.diagnostics
    :members:

ptracers
--------

.. automodule:: darwintools.ptracers
    :members:

.. _density:

density
-------

.. automodule:: darwintools.density
    :members:

miscellaneous utilities
-----------------------

.. automodule:: darwintools.utils
    :members:

.. _plotutils:

plotting utilities
-----------------------

.. automodule:: darwintools.plotutils
    :members:

conversion
----------

.. automodule:: darwintools.conversion
    :members:

.. _cs:

cs
--

.. automodule:: darwintools.cs
    :members:

.. _llc:

llc
---

.. automodule:: darwintools.llc
    :members:

.. _plotting:

plotting tools
--------------

.. automodule:: darwintools.plotting
    :members:

.. _gluemncbig:

gluemncbig
----------

This command line script is part of darwintools and provides a convenient
method for stitching together NetCDF files into a single file covering the
model domain. Be careful though - the resulting files can get very large.

.. program-output:: ../scripts/gluemncbig --help
