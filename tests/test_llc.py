from __future__ import print_function
import sys
from os.path import join as pjoin, dirname
import numpy as np
try:
    import matplotlib as mpl
except ImportError:
    pass
else:
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.testing.compare import compare_images
    try:
        from mpl_toolkits.basemap import Basemap
    except ImportError:
        havebasemap = False
    else:
        havebasemap = True
    import darwintools as mit
    from darwintools.llc.plotting import contourf, contour, pcol

    try:
        import matplotlib.style
    except ImportError:
        pass
    else:
        mpl.style.use('classic')

    TEST_DATA_PATH = pjoin(dirname(__file__), 'data')
    BASELINE_PATH = pjoin(dirname(__file__), 'data', 'baseline_images')
    figsize = (12., 6.)

    def test_contourf(tmp_path):
        x = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'XC'))
        y = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'YC'))
        e = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'Eta'), 8)
        lm = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'hFacC'), lev=0) == 0
        e[lm] = 0.

        fig = plt.figure(figsize=figsize)
        plt.clf()
        h = contourf(x, y, e, np.linspace(-2., 1., 31), cmap = 'jet')
        plt.xlim(-180, 180)
        plt.ylim(-90, 90)
        pngname = tmp_path / 'llc_contourf.png'
        plt.savefig(pngname)
        err = compare_images(pjoin(BASELINE_PATH, pngname), pngname, 13)
        if err:
            raise AssertionError(err)

    def test_contour(tmp_path):
        x = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'XC'))
        y = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'YC'))
        e = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'Eta'), 8)
        lm = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'hFacC'), lev=0) == 0
        e[lm] = 0.

        fig = plt.figure(figsize=figsize)
        plt.clf()
        h = contour(x, y, e, np.linspace(-2., 1., 31), cmap = 'jet')
        plt.xlim(-180, 180)
        plt.ylim(-90, 90)
        pngname = tmp_path / 'llc_contour.png'
        plt.savefig(pngname)
        err = compare_images(pjoin(BASELINE_PATH, pngname), pngname, 13)
        if err:
            raise AssertionError(err)

    def test_pcol(tmp_path):
        x = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'XG'))
        y = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'YG'))
        e = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'Eta'), 8)
        lm = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'hFacC'), lev=0) == 0
        e[lm] = np.nan

        fig = plt.figure(figsize=figsize)
        plt.clf()
        h = pcol(x, y, e, vmin=-2., vmax=1., cmap = 'jet')
        pngname = tmp_path / 'llc_pcol.png'
        plt.savefig(pngname)
        err = compare_images(pjoin(BASELINE_PATH, pngname), pngname, 13)
        if err:
            raise AssertionError(err)

    def test_pcol_xg(tmp_path):
        x = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'XG'))
        y = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'YG'))

        fig = plt.figure(figsize=figsize)
        plt.clf()
        h = pcol(x, y, x, vmin=-180., vmax=180., cmap = 'jet')
        pngname = tmp_path / 'llc_pcol_xg.png'
        plt.savefig(pngname)
        err = compare_images(pjoin(BASELINE_PATH, pngname), pngname, 13)
        if err:
            raise AssertionError(err)

    if havebasemap:
        def test_pcol_stere(tmp_path):
            x = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'XG'))
            y = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'YG'))
            e = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'Eta'), 8)
            lm = mit.rdmds(pjoin(TEST_DATA_PATH, 'llc90', 'hFacC'), lev=0) == 0
            e[lm] = np.nan

            fig = plt.figure(figsize=(8., 8.))
            plt.clf()
            mp = Basemap(-45, 30, 135, 30, projection='stere', lon_0=0., lat_0=90.,
                         resolution = 'l', area_thresh = 1000.)
            h = pcol(x, y, e, mp, vmin=-2., vmax=1., cmap = 'jet')
            mp.fillcontinents(color = 'grey')
            mp.drawmapboundary()
            mp.drawmeridians(np.arange(0, 360, 30))
            mp.drawparallels(np.arange(-90, 90, 30))
            pngname = tmp_path / 'llc_pcol_stere.png'
            plt.savefig(pngname)
            err = compare_images(pjoin(BASELINE_PATH, pngname), pngname, 13)
            if err:
                raise AssertionError(err)

