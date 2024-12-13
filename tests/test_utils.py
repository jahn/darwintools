# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 08:34:30 2023

@author: EGavilan Pascual-Ahuir
"""
from os.path import dirname, join as pjoin
import numpy as np
import darwintools as mit

TEST_DATA_PATH = pjoin(dirname(__file__), 'data')
BASELINE_PATH = pjoin(dirname(__file__), 'data', 'baseline_images')

bathy = np.loadtxt(pjoin(TEST_DATA_PATH, 'utils', 'bathy.txt'), dtype='f4')
rF = np.array([0, -10, -20, -30, -40, -50, -60, -70, -80.01, -90.04, -100.15, -110.47, -121.27,
               -133.03, -146.45, -162.49, -182.31, -207.16, -238.26, -276.68, -323.18,
               -378.18, -441.68, -513.26, -592.16, -677.31, -767.49, -861.45, -958.03, 
               -1056.28, -1155.53, -1255.54, -1356.87, -1461.43, -1572.76, -1695.59, -1834.68,
               -1993.62, -2174.45, -2378.00, -2604.5, -2854, -3126.5], 'f4')

blank_ref = np.r_[10, 11, 12, 24, 32, 43, 44, 45, 51, 56, 57, 62, 63, 64, 68, 79, 90, 91, 102, 103]

def test_blanklist(tmp_path):
    """Example blanklist generator
    """

    # Example 1: Output blanklist without tilemap
    blank=mit.gen_blanklist(bathy, 5, 5, tilemap=False)
    assert np.all(blank == blank_ref)


def test_hfac():
    """Example grid mask generator
    """

    # Example: Output vertical grid mask for the C grid
    [hFacC]=mit.hfac(bathy,rF,0.3,50,'C')
    hFacC_ref = np.fromfile(pjoin(TEST_DATA_PATH, 'utils', 'hFacC.bin'), 'f8')
    hFacC_ref.shape = (42,60,60)
    assert np.all(hFacC == hFacC_ref)


try:
    import matplotlib as mpl
except ImportError:
    have_matplotlib = False
else:
    have_matplotlib = True
    mpl.use('Agg')
    from matplotlib.testing.compare import compare_images
    import darwintools.plotting as mitplot

    def test_blanklist_map(tmp_path):
        """Example blanklist generator
        """

        # Example: Output blanklist with tilemap
        [blank,fig1]=mit.gen_blanklist(bathy, 5, 5, tilemap=True)
        pngname = tmp_path / 'utils_tilemap_blanklist.png'
        fig1.savefig(pngname)
        err = compare_images(pjoin(BASELINE_PATH, pngname), pngname, 13)
        if err:
            raise AssertionError(err)
        assert np.all(blank == blank_ref)


    def test_tilemap(tmp_path):
        """Example tilemap plot distribution
        """

        # Example 1: Output tilemap without specific tile
        fig = mitplot.tilecmap(bathy, 5, 5)
        pngname = tmp_path / 'utils_tilemap.png'
        fig.savefig(pngname)
        err = compare_images(pjoin(BASELINE_PATH, pngname), pngname, 13)
        if err:
            raise AssertionError(err)


    def test_tilemap_zoom(tmp_path):
        """Example tilemap plot distribution
        """

        # Example 2: Output tilemap without specific tile
        fig = mitplot.tilecmap(bathy, 5, 5, 66, sel_zoom=4)
        pngname = tmp_path / 'utils_tilemap_zoom.png'
        fig.savefig(pngname)
        err = compare_images(pjoin(BASELINE_PATH, pngname), pngname, 13)
        if err:
            raise AssertionError(err)
