import numpy as np
import uproot

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from math_utils import *

class SyntheticData:
    def __init__(self, folder, filename, channel):
        """Initialize a synthetic dataset.

        Parameters
        ----------
        folder : str
            Path to the folder containing the ROOT file.
        filename : str
            Name of the ROOT file.
        channel : str
            Name of the ROOT object corresponding to the desired
            interaction channel (see README for the possible options).
        """
        self.folder   = folder
        self.filename = filename
        self.channel  = channel

    def get_rootdata(self):
        """Read the synthetic data from a ROOT file.

        Returns
        -------
        tuple
            Data arrays returned by ``uproot`` through ``to_numpy()``. 
        """
        file = uproot.open(self.folder + self.filename)
        data = file[self.channel]
        return data.to_numpy()
    
    def get_histo(self):
        """Extract the 2D histogram and its bin edges.

        Returns
        -------
        x, y : numpy.ndarray
            Energy, cos(theta) bin edges.
        var_z : numpy.ndarray
            2D histogram values, with shape
            ``(n_cos_theta_bins, n_energy_bins)``.
        """
        data = self.get_rootdata()
        # x, y, z(, b) = number, energy, costheta(, Bjorken-Y)
        x, y, z = data[:3] 
        var_z   = np.squeeze(x).T 
        return y, z, var_z

    def rebin_histo(self, ct_rebin=20, enu_rebin=5, firsts=True):
        """Rebin the 2D histogram by summing neighboring bins.

        Parameters
        ----------
        ct_rebin : int, optional (default=20)
            New number of cos(theta) bins.
        enu_rebin : int, optional (default=5)
            New number of energy bins.
        firsts : bool, optional (default=True)
            If True, retain the first set of bins when trimming down to
            divisible sizes. If False, retain the last set of bins.

        Returns
        -------
        x_rebinned, y_rebinned : numpy.ndarray
            Energy, cos(theta) bin edges after rebinning.
        histo_rebinned : numpy.ndarray
            Rebinned 2D histogram with shape
            ``(ct_rebin, enu_rebin)``.

        Note
        -----
        Only complete groups of ``ct_rebin`` and ``enu_rebin`` bins
        are retained. Any remaining bins that do not form a complete
        group are discarded.
        """
        x, y, z = self.get_histo()
        ct_bins, enu_bins = z.shape

        # trim to divisible sizes
        ct_bins_new   = (ct_bins  // ct_rebin)  * ct_rebin
        enu_bins_new  = (enu_bins // enu_rebin) * enu_rebin
        if firsts:
            # keep the first ct_rebin and enu_rebin bins
            x_trimmed = x[:(enu_bins_new+1)] 
            y_trimmed = y[:(ct_bins_new+1)]
            histo_trimmed = z[:ct_bins_new, :enu_bins_new] 
        else:
            # keep the last ct_rebin and enu_rebin bins
            x_trimmed = x[-(enu_bins_new+1):] 
            y_trimmed = y[-(ct_bins_new+1):]
            histo_trimmed = z[-ct_bins_new:, -enu_bins_new:]

        # reshape
        x_rebinned = x_trimmed[::(enu_bins // enu_rebin)]
        y_rebinned = y_trimmed[::(ct_bins  // ct_rebin)]
        histo_rebinned = histo_trimmed.reshape(ct_bins_new  // ct_rebin,  ct_rebin,
                                               enu_bins_new // enu_rebin, enu_rebin)
        
        # sum
        histo_rebinned = histo_rebinned.sum(axis=(0, 2))

        return x_rebinned, y_rebinned, histo_rebinned

    def regroup_bins_DG(self, n_pois_norm=25, ct_rebin=80):
        """Regroup energy bins according to a Poisson-count threshold.

        Parameters
        ----------
        n_pois_norm : float, optional (default=25)
            Minimum number of events required in every cos(theta) bin
            before an energy bin is finalized.
        ct_rebin : int, optional (default=80)
            New number of cos(theta) bins.

        Returns
        -------
        x_regrouped : numpy.ndarray
            Energy bin edges after regrouping.
        y_regrouped : numpy.ndarray
            cos(theta) bin edges after rebinning.
        histo_regrouped : numpy.ndarray
            Regrouped 2D histogram.

        Notes
        -----
        The regrouping procedure follows the method described in
        arXiv:2408.07015.

        Energy bins are accumulated sequentially. A new energy bin
        is created once the accumulated number of events is at least
        ``n_pois_norm`` in every cos(theta) bin.

        If a final group does not reach the threshold but contains
        non-zero counts in every cos(theta) bin, it is retained as
        the final bin.
        """

        x, y, z = self.rebin_histo(ct_rebin=ct_rebin, enu_rebin=100)
        ct_bins, enu_bins = z.shape

        # Stores the columns of the regrouped histogram (each column corresponds to a merged energy bin).
        histo_regrouped = []

        # Contains the current event counts of all cos(theta) bins, up to the i-th energy bin.
        counts_allct_uptoithe = np.zeros(ct_bins)

        # Upper edges of the regrouped energy bins.
        enu_upper_edges = [x[0]]
        for i in range(enu_bins):
            counts_allct_uptoithe += z[:, i]
            # Regrouping condition
            if (np.amin(counts_allct_uptoithe) >= n_pois_norm):
                histo_regrouped.append(counts_allct_uptoithe.copy())
                enu_upper_edges.append(x[i+1]) 
                counts_allct_uptoithe = np.zeros(ct_bins)

        # Handle leftover last bin 
        if (np.amin(counts_allct_uptoithe) > 0) & (np.amin(counts_allct_uptoithe) < n_pois_norm):
            histo_regrouped.append(counts_allct_uptoithe.copy())
            enu_upper_edges.append(x[-1])

        # Create the regrouped histogram 
        x_regrouped = np.array(enu_upper_edges)
        y_regrouped = y
        histo_regrouped = np.column_stack(histo_regrouped)

        return x_regrouped, y_regrouped, histo_regrouped

    def get_chi2(self, z_expected, n_pois_norm=25, ct_rebin=80):
        """Compute chi-squared values relative to an expected histogram.

        Parameters
        ----------
        z_expected : numpy.ndarray
            Expected 2D histogram used as the reference distribution.
        n_pois_norm : float, optional (default=25)
            Minimum number of events required in every cos(theta) bin
            when regrouping the observed histogram.
        ct_rebin : int, optional (default=80)
            New number of cos(theta) bins.

        Returns
        -------
        numpy.ndarray
            Chi-squared values calculated by comparing the expected
            histogram with the regrouped observed histogram.

        """
        _, _, z_observed = self.regroup_bins_DG(n_pois_norm=n_pois_norm, ct_rebin=ct_rebin)
        chi2_vals = total_chi2(z_expected, z_observed)
        return chi2_vals




















