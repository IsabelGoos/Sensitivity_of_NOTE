import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from constants import *
from math_utils import *

def save_fig(figure, folder_out, filename_out):
    """Save a matplotlib figure to the specified output folder.

    Parameters
    ----------
    figure : matplotlib.figure.Figure
        Figure object to save.
    folder_out : str
        Output folder path.
    filename_out : str
        Name of the output file.

    Returns
    -------
    None
    """

    figure.savefig(f"{folder_out}{filename_out}", dpi=DPI, bbox_inches='tight', transparent=True)

    return None

def truncate_cmap(cmap_name, minval=0.0, maxval=1.0, n=256):
    """Create a truncated version of a matplotlib colormap.

    Parameters
    ----------
    cmap_name : str
        Name of the matplotlib colormap to truncate.
    minval : float, optional (default=0.0)
        Lower bound of the colormap range, between 0.0 and 1.0.
    maxval : float, optional (default=1.0)
        Upper bound of the colormap range, between 0.0 and 1.0.
    n : int, optional (default=256)
        Number of discrete colors used to construct the truncated colormap.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
        Truncated colormap spanning the specified range.
    """

    cmap = plt.get_cmap(cmap_name)
    colors_map = cmap(np.linspace(minval, maxval, n)) 

    return colors.LinearSegmentedColormap.from_list(f"trunc_{cmap_name}", colors_map)

def plot_histo(x, y, z, 
               cmap='cividis', 
               cbartxt=False, 
               cbrticks=False, 
               geolabel=False,
               xlim=False, 
               xlog=False, 
               xticks=False, 
               xticklabels=False, 
               Delta=False):
    """Create a 2D histogram plot with an optional colorbar and twin axis.

    Parameters
    ----------
    x, y, z : array-like
        Values defining the x-axis, y-axis, z-values of the histogram.
    cmap : str, optional (default='cividis')
        Name of the matplotlib colormap used for the histogram.
    cbartxt : str or bool, optional (default=False)
        Label for the colorbar. If False, no label is added.
    cbrticks : array-like or bool, optional (default=False)
        Tick locations for the colorbar. If False, default ticks are used.
    geolabel : bool, optional (default=False)
        If True, the x-axis is time, otherwise it is neutrino energy.
    xlim : tuple or list or bool, optional (default=False)
        Lower and upper limits of the x-axis. If False, default limits are used.
    xlog : bool, optional (default=False)
        If True, use a logarithmic scale for the x-axis.
    xticks : array-like or bool, optional (default=False)
        Locations of the x-axis ticks. If False, default ticks are used.
    xticklabels : array-like or bool, optional (default=False)
        Labels for the x-axis ticks. If False, default labels are used.
    Delta : bool, optional (default=False)
        If True, add a secondary y-axis showing epicentral distance
        in degrees.

    Returns
        -------
    matplotlib.figure.Figure
        Figure containing the 2D histogram.
    """

    # figure
    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
    im      = ax.pcolormesh(x, y, z, cmap=cmap)

    # colorbar
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', location='top')
    if cbartxt:  cbar.set_label(cbartxt, labelpad=10)
    if cbrticks: cbar.set_ticks(cbrticks)
    ##### Probabilities are forced to go from 0 to 1
    ######if (len(x.shape) == 2): im.set_clim(0, 1) 

    # labels
    if geolabel:
        ax.set_xlabel(r"Time [min]")
    else:
        ax.set_xlabel(r"$E_\nu$ [GeV]")
    ax.set_ylabel(r"$\cos \theta_z$")

    # ticks
    if xlim: ax.set_xlim(xlim[0], xlim[1])
    if xlog: ax.set_xscale('log')
    if xticks: ax.set_xticks(xticks)
    if xticklabels: ax.set_xticklabels(xticklabels)

    # twin axis
    if Delta: # Delta = epicentral distance
        ax2 = ax.twinx()
        ax2.set_ylabel(r"$\Delta$ [$^\circ$]")
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks([epi2costh(0),  epi2costh(30),  epi2costh(60),
                        epi2costh(90), epi2costh(120), epi2costh(180)])
        ax2.set_yticklabels(["0", "30", "60", "90", "120", "180"])
    return fig

def plot_binstudy(z, n_pois_norm=25):
    """Plot the results of the binning study.

    Parameters
    ----------
    z : array-like
        2D array containing the logarithm (base 10) of the lowest bin counts.
        The first dimension corresponds to the number of cos(theta) bins
        and the second dimension to the number of energy bins.
    n_pois_norm : float, optional (default=25)
        Reference Poisson count used to define the threshold separating
        the two color scales. The threshold is given by
        ``log10(n_pois_norm)``.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the binning study plot.

    Notes
    -----
    Bins with values above ``log10(n_pois_norm)`` are displayed using a
    truncated ``Blues`` colormap, while bins at or below the threshold
    are displayed using a truncated ``Reds`` colormap.
    """

    # figure
    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(6.5, 7.0), constrained_layout=True)
    (n_ct_bins, n_enu_bins) = np.shape(z)
    X, Y = np.meshgrid(np.arange(n_enu_bins+1), np.arange(n_ct_bins+1), indexing='ij')
    threshold = np.log10(n_pois_norm)
    z_smaller = np.ma.masked_greater(z, threshold)
    z_greater = np.ma.masked_less_equal(z, threshold)
    blue_cmap = truncate_cmap('Blues', 0.4, 1.0)
    red_cmap  = truncate_cmap('Reds',  0.0, 0.4)
    im1 = ax.pcolormesh(X, Y, z_greater.T, cmap=blue_cmap)
    im2 = ax.pcolormesh(X, Y, z_smaller.T, cmap=red_cmap)

    # colorbar
    cbar1 = fig.colorbar(im1, ax=ax, orientation='horizontal', location='top')
    cbar2 = fig.colorbar(im2, ax=ax, orientation='horizontal', location='top')
    cbar1.set_label(r"$\log_{10}$ (Lowest bin count)", labelpad=10)
    cbar2.set_label("", labelpad=20)

    # labels
    ax.set_xlabel(r"Number of energy bins")
    ax.set_ylabel(r"Number of $\cos(\theta)$ bins")
    return fig

