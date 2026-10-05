import numpy as np
import uproot

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from math_utils import total_chi2

class SyntheticNeutrinoData:
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
        Tuple of NumPy arrays with counts, energy, cos(theta), Bjorken-Y
        """
        file = uproot.open(self.folder + self.filename)
        data = file[self.channel]
        return data.to_numpy()
    
    def get_histo(self):
        """Extract the 2D histogram and its bin edges.

        Returns
        -------
        enu_edges, ct_edges : numpy.ndarray
            Energy, cos(theta) bin edges.
        counts : numpy.ndarray
            2D histogram values.  
        """
        data = self.get_rootdata()
        # counts, enu_edges, ct_edges(, Bjorken-Y)
        counts, enu_edges, ct_edges = data[:3] 
        counts = np.squeeze(counts).T 
        return enu_edges, ct_edges, counts

    def rebin_histo(self, ct_rebin=20, enu_rebin=5, firsts=True):
        """Rebin the 2D histogram by summing neighboring bins.

        Parameters
        ----------
        ct_rebin : int, optional (default=20)
            New number of cos(theta) bins.
        enu_rebin : int, optional (default=5)
            New number of energy bins.
        firsts : bool, optional (default=True)
            If True, retain the first bins when trimming down to
            divisible sizes. If False, retain the last bins.

        Returns
        -------
        x_rebinned, y_rebinned : numpy.ndarray
            Energy, cos(theta) bin edges after rebinning.
        histo_rebinned : numpy.ndarray
            Rebinned 2D histogram with shape
            ``(ct_rebin, enu_rebin)``.

        Note
        -----
        Any remaining bins that do not form a complete
        group are discarded.
        """
        x, y, z = self.get_histo()
        ct_bins, enu_bins = z.shape

        # trim to divisible sizes
        ct_bins_new   = (ct_bins  // ct_rebin)  * ct_rebin
        enu_bins_new  = (enu_bins // enu_rebin) * enu_rebin
        if firsts:
            # keep the first ct_bins_new and enu_bins_new bins
            x_trimmed = x[:(enu_bins_new+1)] 
            y_trimmed = y[:(ct_bins_new+1)]
            histo_trimmed = z[:ct_bins_new, :enu_bins_new] 
        else:
            # keep the last ct_bins_new and enu_bins_new bins
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
        """Regroup energy bins to meet a Poisson-count threshold.

        Parameters
        ----------
        n_pois_norm : float, optional (default=25)
            Minimum number of events required in every bin.
        ct_rebin : int, optional (default=80)
            New number of cos(theta) bins.

        Returns
        -------
        x_regrouped, y_regrouped : numpy.ndarray
            Energy, cos(theta) bin edges after regrouping.
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

        # Store the columns of the regrouped histogram 
        # (each column corresponds to a merged energy bin).
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

class SyntheticSeismicData:
    def __init__(self, max_deg=180, step=1, model="ak135f_5s", cache_name="seismic_cache.mseed"):
        """Initialize synthetic seismic datasets from Syngine.

        Parameters
        ----------
        max_deg : int, optional (default=180)
            Maximum angular distance in degrees.
        step : int, optional (default=1)
            Angular step size in degrees.
        model : str, optional (default="ak135f_5s")
            Syngine Earth model.
        cache_name : str, optional (default="seismic_cache.mseed")
            File path to save/load downloaded waveforms to prevent network timeouts.
        """
        self.max_deg    = max_deg
        self.step       = step
        self.model      = model
        self.cache_name = cache_name
        self.st_synth   = None

    def fetch_waveforms(self, force_download=False, max_retries=3):
        """Fetch synthetic seismic waveforms via Syngine."""
        
        if os.path.exists(self.cache_filename) and not force_download:
            print(f"Loading cached seismic waveforms from '{self.cache_filename}'...")
            self.st_synth = obspy.read(self.cache_filename)
            return self.st_synth

        print("Downloading seismic waveforms from Syngine...")
        c_s = SyngineClient()
        self.st_synth = obspy.Stream()
        deg_range = np.arange(0, self.max_deg + 1, self.step)

        for k in deg_range:
            if k % 10 == 0:
                print(f"Fetching distance: {k}°")

            # Retry loop to handle transient network timeouts
            for attempt in range(max_retries):
                try:
                    self.st_synth += c_s.get_waveforms(
                        model=self.model,
                        receiverlatitude=0,
                        receiverlongitude=k,
                        sourcelatitude=0,
                        sourcelongitude=0,
                        sourcedepthinmeters=30000,
                        sourcedoublecouple=[145, 43, 61, 3.51e+24],
                        dt="0.1",
                        units="displacement",
                        components="Z"
                    )
                    break
                except RequestException as e:
                    if attempt < max_retries - 1:
                        time.sleep(2)
                    else:
                        raise e

        # Cache waveforms locally
        self.st_synth.write(self.cache_filename, format="MSEED")
        print(f"Waveforms cached to '{self.cache_filename}'.")
        return self.st_synth

    def process_envelopes(self):
        """Compute the Hilbert transform envelope matrix from the loaded waveforms.

        Returns
        -------
        time_min : numpy.ndarray
            Array of time points in minutes.
        costh : numpy.ndarray
            Array of cos(theta_z) coordinates.
        envelope_mat : numpy.ndarray
            Normalized envelope matrix corresponding to time vs cos(theta_z).
        """
        if self.st_synth is None:
            self.fetch_waveforms()

        n_traces = len(self.st_synth)
        n_pts = len(self.st_synth[0].times())
        envelope_mat = np.zeros([n_traces, n_pts])

        time_min = self.st_synth[0].times() / 60.0
        costh = np.zeros(n_traces)

        for k in range(n_traces):
            data = self.st_synth[k].data
            st_demean = data - np.mean(data)
            env = np.abs(hilbert(st_demean))
            envelope_mat[k, :] = env
            costh[k] = np.cos(np.pi / 2.0 + k * np.pi / 180.0 / 2.0)

        # Scale and clip envelope
        if envelope_mat.max() > 0:
            envelope_mat = (envelope_mat * 5000) / envelope_mat.max()
            envelope_mat = np.clip(envelope_mat, 0, 1)

        return time_min, costh, envelope_mat



















