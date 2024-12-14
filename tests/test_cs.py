import sys
from packaging import version
import pytest
from os.path import join as pjoin, dirname
import numpy as np
try:
    import matplotlib as mpl
except ImportError:
    havematplotlib = False
else:
    havematplotlib = True
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
    from darwintools.cs import pcol

    try:
        import matplotlib.style
    except ImportError:
        pass
    else:
        mpl.style.use('classic')

TEST_DATA_PATH = pjoin(dirname(__file__), 'data')
BASELINE_PATH = pjoin(dirname(__file__), 'data', 'baseline_images')

if havematplotlib:
    def test_pcol(tmp_path):
        ds = mit.rdmnc(pjoin(TEST_DATA_PATH, 'global_ocean.cs32x15',
                             'state.0000072000.t*.nc'),
            ['XG', 'YG', 'Eta'])
        x = ds['XG']
        y = ds['YG']
        e = ds['Eta'][-1]
        e = np.squeeze(e)
        e = np.ma.masked_where(e==0., e)

        fig = plt.figure(figsize=(6.4, 4.8))
        plt.clf()
        h = pcol(x, y, e, cmap = 'jet')
        pngname = tmp_path / 'cs_pcol.png'
        plt.savefig(pngname)
        err = compare_images(pjoin(BASELINE_PATH, 'cs_pcol.png'), pngname, 13)
        if err:
            raise AssertionError(err)


    def test_pcol_sphere(tmp_path):
        ds = mit.rdmnc(pjoin(TEST_DATA_PATH, 'global_ocean.cs32x15',
                             'state.0000072000.*.nc'),
                       ['XG', 'YG', 'Eta'])
        x = ds['XG']
        y = ds['YG']
        e = ds['Eta'][-1]
        e = np.squeeze(e)
        e = np.ma.masked_where(e==0., e)

        fig = plt.figure(figsize=(6.4, 4.8))
        plt.clf()
        h = pcol(x, y, e, projection = 'sphere', cmap = 'jet')
        pngname = tmp_path / 'cs_pcol_sphere.png'
        if version.parse(mpl.__version__) < version.parse('3.3.0'):
            pngname = tmp_path / 'cs_pcol_sphere_pre330.png'
        plt.savefig(pngname)
        err = compare_images(pjoin(BASELINE_PATH, pngname), pngname, 13)
        if err:
            raise AssertionError(err)

    if havebasemap:
        def test_pcol_basemap(tmp_path):
            ds = mit.rdmnc(pjoin(TEST_DATA_PATH, 'global_ocean.cs32x15',
                                 'state.0000072000.*.nc'),
                           ['XG', 'YG', 'Eta'])
            x = ds['XG']
            y = ds['YG']
            e = ds['Eta'][-1]
            e = np.squeeze(e)
            e = np.ma.masked_where(e==0., e)

            fig = plt.figure(figsize=(6.4, 4.8))
            mp = Basemap(projection='moll', lon_0 = 0.,
                         resolution = 'l', area_thresh = 1000.)
            plt.clf()
            h = pcol(x, y, e, projection = mp, cmap = 'jet')
            mp.fillcontinents(color = 'grey')
            mp.drawmapboundary()
            mp.drawmeridians(np.arange(0, 360, 30))
            mp.drawparallels(np.arange(-90, 90, 30))
            pngname = tmp_path / 'cs_pcol_basemap.png'
            plt.savefig(pngname)
            err = compare_images(pjoin(BASELINE_PATH, pngname), pngname, 13)
            if err:
                raise AssertionError(err)

