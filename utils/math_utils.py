import numpy as np

def epi2costh(epi):
    """ 
    Return cos(theta) for a neutrino with epicentral distance epi (degrees). 
    """
    epi_radians  = epi * np.pi / 180.0
    theta_intern = np.pi - (epi_radians + np.pi) / 2.0
    theta        = np.pi - theta_intern
    return np.cos(theta)

def total_chi2(z_expected, z_observed):
    """
    Return total chi2, summed over all energies and incidence angles.
    """
    diff = z_observed - z_expected
    div  = np.divide(np.power(diff, 2), z_expected)
    # The summation over all energies and incidence angles happens here:
    chi2 = np.sum(div[div>0])
    return chi2

