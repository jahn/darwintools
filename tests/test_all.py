try:
    import matplotlib
except ImportError:
    pass
else:
    from darwintools import *

    def test_iolabel():
        for i, s in zip([99, 105, 130, 999, 3800, 3843],
                        ['99', '0f', '0E', 'g7', 'Zi', 'ZZ']):
            assert iolabel(i) == s

        for i in range(1, 3844):
            assert iolabel2num(iolabel(i)) == i
