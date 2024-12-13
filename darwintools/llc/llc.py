from __future__ import print_function
import sys
import numpy as np

def flat(fld, **kwargs):
    """convert mds data into global 2D field
    only fields with 2 to 5 dimensions are allowed"""

    ndims = len(fld.shape)
    if ndims == 2:
        gfld = _flat2D(fld, **kwargs)
    elif ndims == 3:
        gfld = [ _flat2D(fld[a,:,:], **kwargs)
                 for a in range(fld.shape[0]) ]
    elif ndims == 4:
        gfld = [ [ _flat2D(fld[a,b,:,:], **kwargs)
                   for b in range(fld.shape[1]) ]
                 for a in range(fld.shape[0]) ]
    elif ndims == 5:
        gfld = [ [ [ _flat2D(fld[a,b,c,:,:], **kwargs)
                     for c in range(fld.shape[2]) ]
                   for b in range(fld.shape[1]) ]
                 for a in range(fld.shape[0]) ]
    else:
        print("wrong number of dimensions")
        print("only 2 to 5 dimensions are allowed")
        sys.exit(__doc__)

    gfld = np.array(gfld)

    return gfld


def _flat2D(fld, center='Atlantic'):
    """convert mds 2D data into global 2D field"""

    nx = fld.shape[1]
    ny = fld.shape[0]
    n = ny//nx//4

    # eastern and western hemispheres
    eastern=np.concatenate((fld[:n*nx,:],fld[n*nx:2*(n*nx)]),axis=1)
    tmp    = fld[2*(n*nx)+nx:,        ::-1]
    western=np.concatenate((tmp[2::n,:].transpose(),
                            tmp[1::n,:].transpose(),
                            tmp[0::n,:].transpose()))
    # Arctic face is special
    arctic  = fld[2*(n*nx):2*(n*nx)+nx,:]
    arctice = np.concatenate((np.triu(arctic[::-1,:nx//2].transpose()),
                              np.zeros((nx//2,nx))),axis=1)
    # arcticw = np.concatenate((arctic[:,nx:nx//2-1:-1].transpose(),
    #                           np.zeros((nx//2,nx//2)),
    #                           arctic[nx:nx//2-1:-1,nx//2-1::-1]),axis=1)
    mskr = np.tri(nx//2)[::-1,:]
    arcticw = np.concatenate((arctic[0:nx//2,nx:nx//2-1:-1].transpose(),
                              arctic[nx//2:nx,nx:nx//2-1:-1].transpose()*mskr,
                              np.triu(arctic[nx:nx//2-1:-1,nx:nx//2-1:-1]),
                              arctic[nx:nx//2-1:-1,nx//2-1::-1]*mskr),axis=1)
    #
    if center == 'Pacific':
        gfld = np.concatenate( ( np.concatenate((eastern,arctice)),
                                 np.concatenate((western,arcticw)) ), axis=1)
    else:
        gfld = np.concatenate( ( np.concatenate((western,arcticw)),
                                 np.concatenate((eastern,arctice)) ), axis=1)

    return gfld


def _mds2D(fld,center='Atlantic'):
    """convert global 2D 'flat field' to mds 2D data"""

    ni = fld.shape[-1]
    nj = fld.shape[-2]
    nx = ni//4
    ny = nx*(3*4+1)
    n = ny//nx//4

    # arctic face
    arcticw = fld[n*nx:,:nx]
    arctice = fld[n*nx:,2*nx:3*nx]
    arctic = np.concatenate((arctice,arcticw[::-1,::-1]),axis=0)

    # eastern and western hemispheres
    eastern=fld[:n*nx,2*nx:]
    # this is tricky
    western=fld[:n*nx,:2*nx]

    mdsfld = np.concatenate((eastern[:,:nx],
                             eastern[:,nx:],
                             arctic[:,::-1].transpose(),
                             western[::-1,:].transpose().reshape((2*n*nx,nx))),
                             axis=0)
    return mdsfld


def mds(fld,center='Atlantic'):
    """convert global 'flat' field into mds data;
    only fields with 2 to 5 dimensions are allowed"""

    ndims = len(fld.shape)
    if ndims == 2:
        mdsfld = _mds2D(fld, **kwargs)
    elif ndims == 3:
        mdsfld = [ _mds2D(fld[a,:,:], **kwargs)
                 for a in range(fld.shape[0]) ]
    elif ndims == 4:
        mdsfld = [ [ _mds2D(fld[a,b,:,:], **kwargs)
                   for b in range(fld.shape[1]) ]
                 for a in range(fld.shape[0]) ]
    elif ndims == 5:
        mdsfld = [ [ [ _mds2D(fld[a,b,c,:,:], **kwargs)
                     for c in range(fld.shape[2]) ]
                   for b in range(fld.shape[1]) ]
                 for a in range(fld.shape[0]) ]
    else:
        print("wrong number of dimensions")
        print("only 2 to 5 dimensions are allowed")
        sys.exit(__doc__)

    mdsfld = np.array(mdsfld)

    return mdsfld


def faces(fld):
    """convert mds multidimensional data into a list with 6 faces"""

    ndim = len(fld.shape)
    if ndim == 2:
        f = _faces2D(fld)
    else:
        # use list for dynamical memory allocation, because it is fast
        if ndim == 3:
            ff = []
            nk = fld.shape[0]
            for k in range(nk):
                fld2D = fld[k,:,:]
                f2D = _faces2D(fld2D)
                ff.append(f2D)
        elif ndim == 4:
            ff = []
            nk = fld.shape[1]
            nl = fld.shape[0]
            for l in range(nl):
                for k in range(nk):
                    fld2D = fld[l,k,:,:]
                    f2D = _faces2D(fld2D)
                    ff.append(f2D)
        elif ndim == 5:
            ff = []
            nk = fld.shape[2]
            nl = fld.shape[1]
            nm = fld.shape[0]
            for m in range(nm):
                for l in range(nl):
                    for k in range(nk):
                        fld2D = fld[m,l,k,:,:]
                        f2D = _faces2D(fld2D)
                        ff.append(f2D)
        # permute list indices so that face index is the first
        ff = np.transpose(ff)
        f  = []
        for listIndex in range(len(ff)):
            # for each face turn list into array, glue together and reshape
            nx = ff[listIndex][0].shape[-1]
            ny = ff[listIndex][0].shape[-2]
            if   ndim == 3: rshp =       (nk,ny,nx)
            elif ndim == 4: rshp =    (nl,nk,ny,nx)
            elif ndim == 5: rshp = (nm,nl,nk,ny,nx)
            f.append(np.concatenate(np.array(ff[listIndex])).reshape(rshp))

    return f


def faces2mds(ff):
    """convert 6 faces to mds 2D data,
    inverse opertation of llc.faces"""

    ndims = len(ff[0].shape)
    shp = list(ff[0].shape)
    shp[-2]=2*ff[0].shape[-2]
    wd = np.concatenate( (ff[3],ff[4]),axis=-2 ).reshape(shp)
    f  = np.concatenate( (ff[0],ff[1],ff[2],wd),axis=-2)

    return f


def _faces2D(fld):
    """convert mds 2D data into a list with 6 faces"""

    nx = fld.shape[-1]
    ny = fld.shape[-2]
    n = ny//nx//4

    # divide into faces
    f = []
    f.append(fld[:n*nx,:])
    f.append(fld[n*nx:2*(n*nx),:])
    # arctic face
    f.append(fld[2*(n*nx):2*(n*nx)+nx,:])
    # western hemisphere
    wd = fld[2*(n*nx)+nx:,:].reshape(2*nx,n*nx)
    f.append(wd[:nx,:])
    f.append(wd[nx:,:])
    # pseudo-sixth face
    f.append(np.zeros((nx,nx)))

    return f


def _getDims(u,v):
    """
    Retrieve dimensions of input fields and do some sanity checks
    """

    lenu = len(np.shape(u))
    nt = 1
    nk = 1
    ntt = 1
    nkk = 1
    if lenu == 2:
        nju, niu = np.shape(u)
        njv, niv = np.shape(v)
    elif lenu == 3:
        nk, nju, niu = np.shape(u)
        nkk,njv, niv = np.shape(v)
    elif lenu == 4:
        nt, nk, nju, niu = np.shape(u)
        ntt,nkk,njv, niv = np.shape(v)
    else:
        raise ValueError('Can only handle 2 to 4 dimensions')

    if nju!=njv or niu!=niv or nk!=nkk or nt!=ntt:
        raise ValueError('input arrays must have the same dimensions.\n'
                         + 'u.shape = (%i,%i,%i,%i)\n'%(nt,nk,nju,niu)
                         + 'v.shape = (%i,%i,%i,%i)'%(ntt,nkk,njv,niv)
                         )

    if nju!=13*niu:
        raise ValueError('nju=%i not equal 13*niu, niu=%i\n'%(nju,niu) +
                         'This is not and NOT an llc grid.')

    return nt, nk, nju, niu

def div(u, v, dxg=None, dyg=None, rac=None, hfw=None, hfs=None):
    """
    Compute divergence of vector field (U,V) on llc grid

    Call signatures::

       divergence = div(U, V, DXG, DYG, RAC, HFW, HFS)
       divergence = div(U, V)
       divergence = div(U, V, DXG, DYG)
       divergence = div(U, V, DXG, DYG, RAC)
       divergence = div(U, V, DXG, DYG, hfw=HFW, hfs=HFS)

    Parameters
    ----------
    u   : array-like (timelevel,depthlevel,jpoint,ipoint)
          x-component of vector field at u-point

    v   : array-like (timelevel,depthlevel,jpoint,ipoint)
          y-component of vector field at v-point

    dxg : array-like (jpoint,ipoint), optional
          grid spacing in x across v-point, defaults to one

    dyg : array-like (jpoint,ipoint), optional
          grid spacing in y across u-point, defaults to one

    rac : array-like (jpoint,ipoint), optional
          grid cell area, defaults to dxg*dyg

    hfw : array-like (depthlevel,jpoint,ipoint), optional
          hFac at u-point, defaults to one

    hfs : array-like (depthlevel,jpoint,ipoint), optional
          hFac at v-point, defaults to one
    """

    nt, nk, nj, ni =  _getDims(u,v)

    ushape = u.shape

    if dxg is None:
        dxg = np.ones((nj,ni))

    if dyg is None:
        dyg = np.ones((nj,ni))

    if rac is None:
        rac = dxg*dyg

    if hfw is None:
        hfw = np.ones((nk,nj,ni))

    if hfs is None:
        hfs = np.ones((nk,nj,ni))

    u   = u.reshape(nt,nk,nj,ni)
    v   = v.reshape(nt,nk,nj,ni)
    hfw = hfw.reshape(nk,nj,ni)
    hfs = hfs.reshape(nk,nj,ni)

    recip_racf = faces(1./np.where(rac==0.,np.inf,rac))
    divergence = np.zeros(u.shape)
    for t in range(nt):
        for k in range(nk):
            uflx = faces(u[t,k,:,:]*hfw[k,:,:]*dyg)
            vflx = faces(v[t,k,:,:]*hfs[k,:,:]*dxg)
            divf = faces(np.zeros((nj,ni)))
            for iface in range(len(uflx)-1):
                du = np.roll(uflx[iface],-1,axis=-1)-uflx[iface]
                dv = np.roll(vflx[iface],-1,axis=-2)-vflx[iface]
                # now take care of the connectivity
                if iface==0:
                    du[:,-1] = uflx[1][:,   0] - uflx[iface][:,-1]
                    dv[-1,:] = uflx[2][::-1,0] - vflx[iface][-1,:]
                if iface==1:
                    du[:,-1] = vflx[3][0,::-1] - uflx[iface][:,-1]
                    dv[-1,:] = vflx[2][0,:]    - vflx[iface][-1,:]
                if iface==2:
                    du[:,-1] = uflx[3][:,   0] - uflx[iface][:,-1]
                    dv[-1,:] = uflx[4][::-1,0] - vflx[iface][-1,:]
                if iface==3:
                    du[:,-1] = 0. # hack
                    dv[-1,:] = vflx[4][0,:]    - vflx[iface][-1,:]
                if iface==4:
                    du[:,-1] = 0. # hack
                    dv[-1,:] = uflx[0][::-1,0] - vflx[iface][-1,:]
                # putting it all together
                divf[iface] = (du + dv)*recip_racf[iface]

            divergence[t,k,:,:] = faces2mds(divf)

    return divergence.reshape(ushape)


def uv2c(u,v):
    """
    Average vector component (u,v) to center points on llc grid

    Call signatures::

       uc,vc = uv2c(U,V)

    Parameters
    ----------
    U   : array-like (timelevel,depthlevel,jpoint,ipoint)
          x-component of vector field at u-point

    V   : array-like (timelevel,depthlevel,jpoint,ipoint)
          y-component of vector field at v-point

    """

    nt, nk, nj, ni =  _getDims(u,v)
    ushape = u.shape

    u   = u.reshape(nt,nk,nj,ni)
    v   = v.reshape(nt,nk,nj,ni)

    uc = np.zeros(u.shape)
    vc = np.zeros(v.shape)
    for t in range(nt):
        for k in range(nk):
            uf  = faces(u[t,k,:,:])
            vf  = faces(v[t,k,:,:])
            ucf = faces(np.zeros((nj,ni)))
            vcf = faces(np.zeros((nj,ni)))
            for iface in range(len(uf)-1):
                uk = np.roll(uf[iface],-1,axis=-1)+uf[iface]
                vk = np.roll(vf[iface],-1,axis=-2)+vf[iface]
                # now take care of the connectivity
                if iface==0:
                    uk[:,-1] = uf[1][:,   0] + uf[iface][:,-1]
                    vk[-1,:] = uf[2][::-1,0] + vf[iface][-1,:]
                if iface==1:
                    uk[:,-1] = vf[3][0,::-1] + uf[iface][:,-1]
                    vk[-1,:] = vf[2][0,:]    + vf[iface][-1,:]
                if iface==2:
                    uk[:,-1] = uf[3][:,   0] + uf[iface][:,-1]
                    vk[-1,:] = uf[4][::-1,0] + vf[iface][-1,:]
                if iface==3:
                    uk[:,-1] = 0. # hack
                    vk[-1,:] = vf[4][0,:]    + vf[iface][-1,:]
                if iface==4:
                    uk[:,-1] = 0. # hack
                    vk[-1,:] = uf[0][::-1,0] + vf[iface][-1,:]
                #
                ucf[iface] = 0.5*uk
                vcf[iface] = 0.5*vk

            uc[t,k,:,:] = faces2mds(ucf)
            vc[t,k,:,:] = faces2mds(vcf)

    return uc.reshape(ushape), vc.reshape(ushape)


def grad(X, dxc=None, dyc=None, hfw=None, hfs=None):
    """
    Compute horizontal gradient of scalar field X on llc grid

    Call signatures::

       dXdx, dXdy = div(X, DXC, DYC, HFW, HFS)
       dXdx, dXdy = div(X)
       dXdx, dXdy = div(X, DXC, DYC)
       dXdx, dXdy = div(X, hfw=HFW, hfs=HFS)

    Parameters
    ----------
    X   : array-like (timelevel,depthlevel,jpoint,ipoint)
          scalar field at c-point

    dxc : array-like (jpoint,ipoint), optional
          grid spacing in x across u-point, defaults to one

    dyc : array-like (jpoint,ipoint), optional
          grid spacing in y across v-point, defaults to one

    hfw : array-like (depthlevel,jpoint,ipoint), optional
          hFac at u-point, defaults to one

    hfs : array-like (depthlevel,jpoint,ipoint), optional
          hFac at v-point, defaults to one
    """

    nt, nk, nj, ni =  _getDims(X,X)
    if dxc is None:
        dxc = np.ones((nj,ni))

    if dyc is None:
        dyc = np.ones((nj,ni))

    if hfw is None:
        hfw = np.ones((nk,nj,ni))

    if hfs is None:
        hfs = np.ones((nk,nj,ni))

    mskw=np.where(hfw>0.,1.,0.).reshape(nk,nj,ni)
    msks=np.where(hfs>0.,1.,0.).reshape(nk,nj,ni)

    xshape = X.shape

    X   = X.reshape(nt,nk,nj,ni)

    dXdx = np.zeros(X.shape)
    dXdy = np.zeros(X.shape)

    rdxc = faces(1./np.where(dxc==0.,np.inf,dxc))
    rdyc = faces(1./np.where(dyc==0.,np.inf,dyc))
    for t in range(nt):
        for k in range(nk):
            xf  = faces(X[t,k,:,:])
            duf = faces(np.zeros((nj,ni)))
            dvf = faces(np.zeros((nj,ni)))
            for iface in range(len(xf)-1):
                du = (xf[iface] - np.roll(xf[iface],1,axis=-1))*rdxc[iface]
                dv = (xf[iface] - np.roll(xf[iface],1,axis=-2))*rdyc[iface]
                # now take care of the connectivity
                if iface==0:
                    du[:,0] = (xf[0][:,0] - xf[4][-1,::-1])*rdxc[0][:,0]
                    dv[0,:] = 0. # hack
                if iface==1:
                    du[:,0] = (xf[1][:,0] - xf[0][:,   -1])*rdxc[1][:,0]
                    dv[0,:] = 0. # hack
                if iface==2:
                    du[:,0] = (xf[2][:,0] - xf[0][-1,::-1])*rdxc[2][:,0]
                    dv[0,:] = (xf[2][0,:] - xf[1][-1,   :])*rdxc[2][0,:]
                if iface==3:
                    du[:,0] = (xf[3][:,0] - xf[2][:   ,-1])*rdxc[3][:,0]
                    dv[0,:] = (xf[3][0,:] - xf[1][::-1,-1])*rdyc[3][0,:]
                if iface==4:
                    du[:,0] = (xf[4][:,0] - xf[2][-1,::-1])*rdxc[4][:,0]
                    dv[0,:] = (xf[4][0,:] - xf[3][-1,   :])*rdyc[4][0,:]

                duf[iface]=du
                dvf[iface]=dv

            dXdx[t,k,:,:] = faces2mds(duf)*mskw[k,:,:]
            dXdy[t,k,:,:] = faces2mds(dvf)*msks[k,:,:]

    return dXdx.reshape(xshape), dXdy.reshape(xshape)
