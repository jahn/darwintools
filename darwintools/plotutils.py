# Created by EGavilan Pascual-Ahuir on 2023-04-07
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1.inset_locator import (mark_inset,inset_axes,
                                                   zoomed_inset_axes)
__doc__ = """
Plotting utilities for MITgcm.
"""

cmap_lm =  ListedColormap(["lightsteelblue"])

def tilecmap(arr, sNx, sNy, tilen=None, sel_zoom=5, fill_value=0):
    """
    Pseudocolor plot of land mask with tiles superimposed, optionally
    showing the values of arr for a single tile.

    Parameters
    ----------
    arr        : 2D array_like
                 values to plot, land mask is taken as arr==fill_value
    sNx        : int
                 number of x points in each tile
    sNy        : int
                 number of y points in each tile
    tilen      : int or None
                 plot a specific tile, default None
    sel_zoom   : int
                 zooming range, default 5
    fill_value : float
                 default 0

    Returns
    -------
    figure     : matplotlib figure
                 Plot of land mask, tiles and arr values

    Usage
    -----
    >>> [fig]=tilecmap(bathy, 5, 5)
    >>> [fig]=tilecmap(bathy, 5, 5, 66, sel_zoom=4)

    """

    #Check dimensions of arr
    assert arr.ndim==2,'check_stp: array must be 2D'

    [Ny,Nx]=arr.shape

    nTx = Nx//sNx
    nTy = Ny//sNy

    mland = np.copy(arr)
    mland[mland!=fill_value] = 1
    mland[mland==fill_value] = np.nan
    arr = arr*mland

    tile_order = np.zeros([nTy, nTx], dtype=int)

    for n in range(0, nTy):
        for m in range(0, nTx):
            tile_order[n, m] = int(n*nTx+m+1)

    [cn_x, cn_y] = np.meshgrid(np.arange(sNx//2, Nx, sNx),
                               np.arange(sNy//2, Ny, sNy))

    major_xticks = np.arange(0, Nx+sNx, sNx)
    major_yticks = np.arange(0, Ny+sNy, sNy)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.pcolormesh(mland,cmap=cmap_lm)

    for a, b, c in zip(cn_x.flat, cn_y.flat, tile_order.flat):
        ax.annotate(str(c), (a, b), color='black',
                    ha='center', va='center')

    if tilen is not None:

        [Tind]=np.argwhere(tile_order==tilen)

        #Select position for the zoom inseting
        if (Tind[0]>nTy/2):
            if(Tind[1]>nTx//2):locz=2; locm1=1; locm2=4; cbar_pos=-0.1
            else:locz=1; locm1=2; locm2=3; cbar_pos=1.05
        else:
            if(Tind[1]>nTx//2):locz=1; locm1=3; locm2=4; cbar_pos=1.05
            else:locz=2; locm1=3; locm2=4; cbar_pos=-0.1
   
        #Select colorbar range for zoom
        Tix = Tind[1]*sNx;  Tiy = Tind[0]*sNy
        arrmin = np.nanmin(arr[Tiy:Tiy+sNy,Tix:Tix+sNx])
        arrmax = np.nanmax(arr[Tiy:Tiy+sNy,Tix:Tix+sNx])
   
        ax2 = zoomed_inset_axes(ax, zoom=sel_zoom, loc=locz, borderpad=-1)
        pc=ax2.pcolormesh(arr,vmin=arrmin,vmax=arrmax,cmap=plt.cm.jet)
   
        ax2.set_xlim([major_xticks[Tind[1]],major_xticks[Tind[1]]+sNx])
        ax2.set_ylim([major_yticks[Tind[0]],major_xticks[Tind[0]]+sNy])
        mark_inset(ax, ax2, loc1=locm1, loc2=locm2, fc="none", lw=1.5, ec='0')
        plt.xticks(visible=False)
        plt.yticks(visible=False)
        cax = inset_axes(ax2,
                         width="5%",
                         height="100%",
                         loc="lower left",
                         bbox_to_anchor=(cbar_pos,0,1,1),
                         bbox_transform=ax2.transAxes,
                         borderpad=0,
                         )
        cbar=plt.colorbar(pc,cax=cax, orientation='vertical')
        if cbar_pos<0:
            cax.yaxis.tick_left()

    ax.set_xticks(major_xticks)
    ax.set_yticks(major_yticks)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([0,Nx])
    ax.set_ylim([0,Ny])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.grid()

    return fig
