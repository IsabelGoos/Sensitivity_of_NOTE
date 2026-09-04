import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from constants import DPI, wt_rocks, elements
from math_utils import epi2costh

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
    minval, maxval : float, optional (default=0.0, 1.0)
        Lower, upper bound of the colormap range, between 0.0 and 1.0.
    n : int, optional (default=256)
        Number of discrete colors used to construct the truncated colormap.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
        Truncated colormap spanning the specified range.
    """

    cmap = plt.get_cmap(cmap_name)
    trunc_map = cmap(np.linspace(minval, maxval, n)) 

    return colors.LinearSegmentedColormap.from_list(f"trunc_{cmap_name}", trunc_map)

def plot_histo(enu_edges, ct_edges, counts, 
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
    enu_edges, ct_edges, counts : array-like
        Values defining the x-axis, y-axis, z-values of the histogram.
    cmap : str, optional (default='cividis')
        Matplotlib colormap used for the histogram.
    cbartxt : str or bool, optional (default=False)
        Label for the colorbar. If False, no label is added.
    cbrticks : array-like or bool, optional (default=False)
        Tick locations for the colorbar. If False, default ticks are used.
    geolabel : bool, optional (default=False)
        If True, the x-axis shows time, otherwise it shows neutrino energy.
    xlim : tuple or list or bool, optional (default=False)
        Lower and upper limits of the x-axis. If False, default limits are used.
    xlog : bool, optional (default=False)
        If True, use a logarithmic scale for the x-axis.
    xticks : array-like or bool, optional (default=False)
        Tick locations for the x-axis. If False, default ticks are used.
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

    # Create figure
    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
    im      = ax.pcolormesh(enu_edges, ct_edges, counts, cmap=cmap)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', location='top')
    if cbartxt:  cbar.set_label(cbartxt, labelpad=10)
    if cbrticks: cbar.set_ticks(cbrticks)
    ##### Probabilities are forced to go from 0 to 1
    ######if (len(x.shape) == 2): im.set_clim(0, 1) 

    # Labels
    if geolabel:
        ax.set_xlabel(r"Time [min]")
    else:
        ax.set_xlabel(r"$E_\nu$ [GeV]")
    ax.set_ylabel(r"$\cos \theta_z$")

    # x-axis
    if xlim: ax.set_xlim(xlim[0], xlim[1])
    if xlog: ax.set_xscale('log')
    if xticks:      ax.set_xticks(xticks)
    if xticklabels: ax.set_xticklabels(xticklabels)

    # Twin axis
    if Delta: # Delta = epicentral distance
        ax2 = ax.twinx()
        ax2.set_ylabel(r"$\Delta$ [$^\circ$]")
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks([epi2costh(0),  epi2costh(30),  epi2costh(60),
                        epi2costh(90), epi2costh(120), epi2costh(180)])
        ax2.set_yticklabels(["0", "30", "60", "90", "120", "180"])

    return fig

def plot_binstudy(counts, n_pois_norm=25):
    """Plot the results of the binning study.

    Parameters
    ----------
    counts : array-like
        2D array containing the logarithm (base 10) of the lowest bin counts.
        The first, second dimension corresponds to the number of cos(theta), 
        energy bins.
    n_pois_norm : float, optional (default=25)
        Reference Poisson count used to define the threshold separating
        the two color scales. The threshold is given by ``log10(n_pois_norm)``.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the binning study plot.

    Notes
    -----
    Bins with values above ``log10(n_pois_norm)`` are displayed using a
    truncated ``Blues`` colormap, otherwise using a truncated ``Reds`` colormap.
    """

    # Bin coordinates
    (n_ct_bins, n_enu_bins) = np.shape(counts)
    X, Y = np.meshgrid(np.arange(n_enu_bins+1), np.arange(n_ct_bins+1), indexing='ij')

    # Split histogram at the Poisson threshold
    threshold = np.log10(n_pois_norm)
    z_smaller = np.ma.masked_greater(counts, threshold)
    z_greater = np.ma.masked_less_equal(counts, threshold)

    # Truncated colormaps
    blue_cmap = truncate_cmap('Blues', 0.4, 1.0)
    red_cmap  = truncate_cmap('Reds',  0.0, 0.4)

    # Create figure
    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(6.5, 7.0), constrained_layout=True) 
    im1 = ax.pcolormesh(X, Y, z_greater.T, cmap=blue_cmap)
    im2 = ax.pcolormesh(X, Y, z_smaller.T, cmap=red_cmap)

    # Colorbar
    cbar1 = fig.colorbar(im1, ax=ax, orientation='horizontal', location='top')
    cbar2 = fig.colorbar(im2, ax=ax, orientation='horizontal', location='top')
    cbar1.set_label(r"$\log_{10}$ (Lowest bin count)", labelpad=10)
    cbar2.set_label("", labelpad=20)

    # Labels
    ax.set_xlabel(r"Number of energy bins")
    ax.set_ylabel(r"Number of $\cos(\theta)$ bins")
    
    return fig

def plot_Ye(Z, A, labels, xlim=(0, 36), ylim_l=(0.43, 0.505), ylim_h=(0.98, 1.0)):
    """Plot the electron yield Y_e for a set of elements and composites."""

    # Extract rock weights and nuclear properties
    wt_BE   = np.array([v["BE"]   for v in wt_rocks.values()])
    wt_core = np.array([v["core"] for v in wt_rocks.values()])
    wt_BSE  = np.array([v["BSE"]  for v in wt_rocks.values()])
    Z_rocks = np.array([elements[element][0] for element in wt_rocks])
    A_rocks = np.array([elements[element][1] for element in wt_rocks])

    # Compute weighted atmoic numbers and mean electron yields for each composite
    Z_BE   = np.average(Z_rocks, weights=wt_BE)
    Z_core = np.average(Z_rocks, weights=wt_core)
    Z_BSE  = np.average(Z_rocks, weights=wt_BSE)
    Ye_BE   = np.average(Z_rocks/A_rocks, weights=wt_BE)
    Ye_core = np.average(Z_rocks/A_rocks, weights=wt_core)
    Ye_BSE  = np.average(Z_rocks/A_rocks, weights=wt_BSE)

    # Create figure
    plt.rcParams.update({'font.size': 16})
    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        sharex=True,
        figsize=(8, 6),
        constrained_layout=True,
        gridspec_kw={'height_ratios': [1, 5]}
    )

    # Elements, the first has to be hydrogen
    for i, (z, a, label) in enumerate(zip(Z, A, labels)):
        ax = ax1 if i == 0 else ax2
        ax.text(z, z/a, label, fontsize=14,
                ha='center', va='center')

    # Plot composite points
    composites = {
        "BE":   (Z_BE,   Ye_BE),
        "core": (Z_core, Ye_core),
        "BSE":  (Z_BSE,  Ye_BSE),
    }

    for label, (z, ye) in composites.items():
        ax2.scatter(
            z, ye,
            color='tab:red',
            s=40,
            zorder=3,
            marker='x'
        )

    # Axis limits
    ax1.set_xlim(*xlim)
    ax1.set_ylim(*ylim_h)
    ax2.set_ylim(*ylim_l)

    # Broken-axis appearance
    ax1.spines.bottom.set_visible(False)
    ax2.spines.top.set_visible(False)
    ax1.tick_params(axis='x', bottom=False)
    ax2.tick_params(axis='x', top=False)

    # Diagonal marks indicating the broken y-axis
    d = 0.008
    range1 = np.diff(ax1.get_ylim())[0]
    range2 = np.diff(ax2.get_ylim())[0]

    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot(
        (-d, d),
        (-d*range2/range1, d*range2/range1),
        **kwargs
    )
    ax1.plot(
        (1-d, 1+d),
        (-d*range2/range1, d*range2/range1),
        **kwargs
    )

    kwargs.update(transform=ax2.transAxes)
    ax2.plot((-d, d), 
             (1-d, 1+d), 
             **kwargs)
    ax2.plot((1-d, 1+d), 
             (1-d, 1+d), 
             **kwargs
    )

    # Labels
    ax2.set_xlabel(r"$Z$")
    ax2.set_ylabel(r"$Y_e$")

    return fig