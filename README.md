# Darwintools

Python tools for the Darwin ecosystem model, see the
[darwinproject](https://darwinproject.mit.edu/), and the
MIT General Circulation Model, [MITgcm](http://mitgcm.org), in general.

The Darwin model is documented
[here](http://darwin3.rtfd.io/en/latest/phys_pkgs/darwin.html)
as a part of the MITgcm.

This package is developed on [github](https://github.com/jahn/darwintools).

To test any changes to darwintools, do the following
1. Check out the submodule with reference data:
```
    git submodule init
    git submodule update
```
2. Install the needed python versions (python3.7 through python3.12).
3. Install tox: https://tox.wiki/en/4.15.0/installation.html.
   On a mac, it may be necessary to install pipx via macports/brew
   or use the "pip" method.
4. Run "tox".