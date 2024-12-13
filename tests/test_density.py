from os.path import join as pjoin, dirname
import numpy as np
from darwintools import density

TEST_DATA_PATH = pjoin(dirname(__file__), 'data')

def test_jmd95():
    rho = density.jmd95(35.5, 3., 3000.)
    assert abs(rho - 1041.83267) < .00001
    # 1041.8326696373254

def test_mdjwf():
    rho = density.mdjwf(35., 25., 2000.)
    assert abs(rho - 1031.654229) < .000001
