from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from .llc import *

def contourf(*arguments, **kwargs):
    """
    Create a contourf plot of a 2-D llc array (with tricontour).

    Call signatures::

        contourf(X, Y, C, N, **kwargs)

        contourf(X, Y, C, V, **kwargs)

    Parameters
    ----------
    X : array-like
        x coordinates of the grid points

    Y : array-like
        y coordinates of the grid points

    C : array-like
        array of color values.

    N : int
        number of levels

    V : list of float
        list of levels

    kwargs
        passed to tricontour.

    """

    arglen = len(arguments)
    h = []
    if arglen >= 3:
        data = np.copy(arguments[2].flatten())
        x = arguments[0].flatten()
        y = arguments[1].flatten()

        # Create the Triangulation;
        # no triangles so Delaunay triangulation created.
        triang = tri.Triangulation(x, y)
        ntri = triang.triangles.shape[0]

        # Mask off unwanted triangles.
        mask = np.where(data[triang.triangles].prod(axis=1)==0., 1, 0)
        triang.set_mask(mask)

        if arglen == 3:
            h = plt.tricontourf(triang, data, **kwargs)
        elif arglen == 4:
            h = plt.tricontourf(triang, data, arguments[3], **kwargs)
        else:
            print("wrong number of arguments")
            print("need at least 3 or 4 arguments")
            sys.exit(__doc__)

        # show the triangles for debugging
        #plt.triplot(triang, color='0.7')

    else:
        print("wrong number of arguments")
        print("need at least x,y,fld")
        sys.exit(__doc__)

    return h


def contour(*arguments, **kwargs):
    """
    Create a contour plot of a 2-D llc array (with tricontour).

    Call signatures::

        contour(X, Y, C, N, **kwargs)

        contour(X, Y, C, V, **kwargs)

    Parameters
    ----------
    X : array-like
        x coordinates of the grid points

    Y : array-like
        y coordinates of the grid points

    C : array-like
        array of color values.

    N : int
        number of levels

    V : list of float
        list of levels

    kwargs
        passed to tricontour.

    """

    arglen = len(arguments)
    h = []
    if arglen >= 3:
        data = arguments[2].flatten()
        x = arguments[0].flatten()
        y = arguments[1].flatten()

        # Create the Triangulation;
        # no triangles so Delaunay triangulation created.
        triang = tri.Triangulation(x, y)
        ntri = triang.triangles.shape[0]

        # Mask off unwanted triangles.
        mask = np.where(data[triang.triangles].prod(axis=1)==0., 1, 0)
        triang.set_mask(mask)

        if arglen == 3:
            h = plt.tricontour(triang, data, **kwargs)
        elif arglen == 4:
            h = plt.tricontour(triang, data, arguments[3], **kwargs)
        else:
            print("wrong number of arguments")
            print("need at least 3 or 4 arguments")
            sys.exit(__doc__)

        # show the triangles for debugging
        #plt.triplot(triang, color='0.7')

    else:
        print("wrong number of arguments")
        print("need at least x,y,fld")
        sys.exit(__doc__)

    return h


def _sqCoord(a):
    b = np.squeeze(a)
    return b


def _sqData(a):
    b = np.copy(np.squeeze(a))
    b = np.ma.masked_where(b==0., b)
    b = np.ma.masked_where(np.isnan(b), b)
    return b


def _pcolormesh(x, y, d, **kwargs):
    '''
    Replace non-finite coordinates and mask data before calling pcolormesh
    to deal with some Basemap projection issues.
    '''
    j, i = np.where(~(np.isfinite(x)&(np.isfinite(y))))
    if len(i):
        ny, nx = d.shape
        iclip = np.clip(i,0,nx-1)
        jclip = np.clip(j,0,ny-1)
        im1clip = np.clip(i-1,0,nx-1)
        jm1clip = np.clip(j-1,0,ny-1)
        d[jclip,iclip] = np.NaN
        d[jclip,im1clip] = np.NaN
        d[jm1clip,iclip] = np.NaN
        d[jm1clip,im1clip] = np.NaN
        x[j,i] = 0
        y[j,i] = 0

    return plt.pcolormesh(x, y, d, **kwargs)


def pcol(*arguments, **kwargs):
    """
    Create a pseudo-color plot of a 2-D llc array (with plt.pcolormesh).

    Call signatures::

        pcol(X, Y, C, **kwargs)

        pcol(X, Y, C, m, **kwargs)

    Parameters
    ----------
    X : array-like
        x coordinates of the grid point corners (G-points)

    Y : array-like
        y coordinates of the grid point corners (G-points)

    C : array-like
        array of color values.

    m : Basemap instance, optional
        map projection to use.
        NOTE: currently not all projections work

    kwargs
        passed to plt.pcolormesh.

    """

    arglen = len(arguments)
    h = []
    mapit = False
    if arglen < 3:
        print("wrong number of arguments")
        print("need at least x,y,fld")
        sys.exit(__doc__)
    elif arglen > 3:
        mapit = True
        m = arguments[3]

    if mapit:
        # not all projections work, catch few of these here
        if ( (m.projection == 'hammer') |
             (m.projection == 'robin')  |
             (m.projection == 'moll')   |
             (m.projection == 'cea') ):
            sys.exit("selected projection '"+m.projection
                     +"' is not supported")

        # these projections use simple code for the Arctic face;
        # all others require more complicted methods
        stereographicProjection = (m.projection == 'npaeqd')  | \
                                  (m.projection == 'spaeqd')  | \
                                  (m.projection == 'nplaea')  | \
                                  (m.projection == 'splaea')  | \
                                  (m.projection == 'npstere') | \
                                  (m.projection == 'spstere') | \
                                  (m.projection == 'stere')
    else:
        stereographicProjection = False


    xg = arguments[0]
    yg = arguments[1]
    data = arguments[2]

    nx = data.shape[-1]
    ny = data.shape[-2]
    n = ny//nx//4

    # color range
    cax = [data.min(),data.max()]
    # overwrite if necessary
    if 'vmin' in kwargs: cax[0] = kwargs.pop('vmin','')
    if 'vmax' in kwargs: cax[1] = kwargs.pop('vmax','')
    # divide into faces
    f0 = []
    f0.append(faces(xg))
    f0.append(faces(yg))
    f0.append(faces(data))
    # fill holes in coordinate arrays
#    for t in [0,1,3,4]:
#        inan = f0[2][t]==0 # _sqCoord(f0[2][t])==np.nan]
#        f0[0][t][inan]=np.nan
#        f0[1][t][inan]=np.nan

#    for t in [0,1]:
#        for i in range(nx):
#            for j in range(n*nx):
#                if f0[0][t][j,i]==0:f0[0][t][100,i]
#                if f0[1][t][j,i]==0:f0[1][t][100,i]
#
#    for t in [3,4]:
#        for i in range(n*nx):
#            for j in range(nx):
#                if f0[0][t][j,i]==0:f0[0][t][j,239]
#                if f0[1][t][j,i]==0:f0[1][t][j,239]

    # find the missing corners by interpolation
    fo = []
    fo.append( (f0[0][0][-1,0]+f0[0][2][-1,0]+f0[0][4][-1,0])/3. )
    fo.append( (f0[1][2][-1,0]+f0[1][2][-1,0]+f0[1][4][-1,0])/3. )
    fo.append( np.nan )
    fe = []
    fe.append( (f0[0][1][0,-1]+f0[0][3][0,-1])/2. )
    fe.append( (f0[1][1][0,-1]+f0[1][3][0,-1])/2. )
    fe.append( np.nan )
    f = np.array(f0, dtype=object)
    # fill some gaps at the face boundaries, but only for the coordinate arrays (k=0,1)
    for t in [0,2,4]:
        tp = 2*(t//2)
        tpp = tp
        if tp==4: tpp = tp-6
        for k in [0,1]:
            tp = min(tp,3)
            f[k][t] = np.concatenate((f0[k][t],f0[k][1+tp][:,:1]),axis=1)
            if k==2: tmp = np.atleast_2d(np.append(f0[k][2+tpp][::-1,:1],fo[k]))
            else:    tmp = np.atleast_2d(np.append(fo[k],f0[k][2+tpp][::-1,:1]))
            f[k][t] = np.concatenate((f[k][t],tmp),axis=0)

    for t in [1,3]:
        tp = 2*(t//2)
        for k in [0,1]:
            f[k][t] = np.concatenate((f0[k][t],f0[k][2+tp][:1,:]),axis=0)
            if k==2: tmp = np.atleast_2d(np.append(f0[k][3+tp][:1,::-1],fe[k]))
            else:    tmp = np.atleast_2d(np.append(fe[k],f0[k][3+tp][:1,::-1]))
            f[k][t] = np.concatenate((f[k][t],tmp.transpose()),axis=1)

    # we do not really have a sixth face so we overwrite the southernmost row
    # of face 4 and 5 by a hack:
    for t in [3,4]:
        f[0][t][:,-1] = f[0][t][:,-2]
        f[1][t][:,-1] = -90. # degree = south pole

    # make sure that only longitudes of one sign are on individual lateral faces
    i0 = f[0][3]<0.
    f[0][3][i0] = f[0][3][i0]+360.
    # plot the lateral faces
    ph = []
    for t in [0,1,3,4]:
        if mapit: x, y = m(_sqCoord(f[0][t]), _sqCoord(f[1][t]))
        else:     x, y =   _sqCoord(f[0][t]), _sqCoord(f[1][t])
        ph.append(_pcolormesh(x,y,_sqData(f[2][t]), **kwargs))
    # plot more lateral faces to be able to select the longitude range later
    for t in [1,3,4]:
        f[0][t] = f[0][t]+ (-1)**t*360.
        if mapit: x, y = m(_sqCoord(f[0][t]), _sqCoord(f[1][t]))
        else:     x, y =   _sqCoord(f[0][t]), _sqCoord(f[1][t])
        ph.append(_pcolormesh(x,y,_sqData(f[2][t]), **kwargs))

    # Arctic face is special, because of the rotation of the grid by
    # rangle = 7deg (seems to be the default)
    t = 2

    if mapit & stereographicProjection:
        x, y = m(_sqCoord(f[0][t]),_sqCoord(f[1][t]))
        ph.append(_pcolormesh(x,y,_sqData(f[2][t]), **kwargs))
    else:
        rangle = 7.
        # first half of Arctic tile
        nn = nx//2+1
        xx = np.copy(f[0][t][:nn,:])
        yy = np.copy(f[1][t][:nn,:])
        # make sure that the data have one columns/rows fewer that the coordinates
        zz = np.copy(f[2][t][:nn-1,:])
        xx = np.where(xx<rangle,xx+360,xx)
        if mapit: x, y = m(_sqCoord(xx),_sqCoord(yy))
        else:     x, y =   _sqCoord(xx),_sqCoord(yy)
        ph.append(_pcolormesh(x,y,_sqData(zz), **kwargs))
        # repeat for xx-360
        xx = xx-360.
        if mapit: x, y = m(_sqCoord(xx),_sqCoord(yy))
        else:     x, y =   _sqCoord(xx),_sqCoord(yy)
        ph.append(_pcolormesh(x,y,_sqData(zz), **kwargs))
        # second half of Arctic tile
        nn = nx//2
        xx = np.copy(f[0][t][nn:,:])
        yy = np.copy(f[1][t][nn:,:])
        zz = np.copy(f[2][t][nn:,:])
        #
        if mapit: x, y = m(_sqCoord(xx),_sqCoord(yy))
        else:     x, y =   _sqCoord(xx),_sqCoord(yy)
        ph.append(_pcolormesh(x,y,_sqData(zz), **kwargs))
        # repeat for xx+360
        xx = xx + 360.
        if mapit: x, y = m(_sqCoord(xx),_sqCoord(yy))
        else:     x, y =   _sqCoord(xx),_sqCoord(yy)
        ph.append(_pcolormesh(x,y,_sqData(zz), **kwargs))

    if not mapit:
        plt.xlim([-170,190])
        plt.ylim([-90,90])

    for im in ph:
        im.set_clim(cax[0],cax[1])

    return ph
