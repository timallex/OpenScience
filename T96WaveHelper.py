
from pyspedas import  get_data, tplot
from pyspedas.projects.themis import fft
import numpy as np




def getBWVals(probe, selectedTime, trange, noUpdate, plotBool,timeW1,dataW1,freqs,timeW2,dataW2,eW2):
    """Calculate wave magnetic field amplitudes from THEMIS spectral data.

    Parameters
    ----------
    probe : str
        THEMIS probe identifier: "a", "b", "c", "d", or "e".
    selectedTime : float
        Selected time as a Unix timestamp in seconds.
    trange : list of str
        Start and end dates for data retrieval.
    noUpdate : bool
        Whether to avoid downloading updated data.
    plotBool : bool
        Whether to display the wave spectral data plot.
    timeW1 : array-like
        Timestamps for the first magnetic field spectral dataset.
    dataW1 : array-like
        Power spectral density values for the first dataset in nT^2/Hz.
    freqs : array-like
        Center frequencies of the spectral bins in Hz.
    timeW2 : array-like
        Timestamps for the second magnetic field spectral dataset.
    dataW2 : array-like
        Power spectral density values for the second dataset in nT^2/Hz.
    eW2 : array-like
        Frequency values for the second spectral dataset in Hz.

    Returns
    -------
    freqs : numpy.ndarray
        Center frequencies of the spectral bins in Hz.
    Bw : numpy.ndarray
        Wave magnetic field amplitudes in nT.
    omegas : numpy.ndarray
        Wave angular frequencies in rad/s.
    """


    print(len(timeW1),len(timeW2),len(dataW2),len(eW2), len(dataW2[0]))

    print(timeW1)

    idx_timeW1 = np.argmin(np.abs(np.array(timeW1) - selectedTime))

    print(idx_timeW1)
    print("Closest time:", timeW1[idx_timeW1], selectedTime,(timeW1[idx_timeW1] - selectedTime))
    print(dataW2[idx_timeW1])
    print(len(dataW2[idx_timeW1]))

    dataW = np.array(dataW1[idx_timeW1]) + np.array(dataW2[idx_timeW1])


    print(dataW)

    print(freqs,eW2)


    # freqs = the 32 frequency centers, like [4, 20, 36, ..., 3836]
    freqs = np.array(freqs, dtype=float)

    # make 33 edges for 32 bins
    edges = np.zeros(len(freqs) + 1)

    # first bin starts at 0 Hz
    edges[0] = 0

    # middle edges are halfway between neighboring frequency centers
    edges[1:-1] = 0.5 * (freqs[:-1] + freqs[1:])

    # last bin extends past the final center by half the last spacing
    edges[-1] = freqs[-1] + 0.5 * (freqs[-1] - freqs[-2])

    # bin sizes
    bin_widths = np.diff(edges)

    print(bin_widths)

    # convert PSD to Bw
    Bw = np.sqrt(dataW * bin_widths)


    omegas = freqs * np.pi * 2

    if plotBool: tplot(f'th{probe}_fff_32_scm2')

    return freqs, Bw, omegas


def select_wave_values(
        measured_omega_values,
        measured_bw_values,
        omega_ce_eq,
        omega_ratios
):
    """Select requested frequencies and interpolate their wave amplitudes.

    Parameters
    ----------
    measured_omega_values : array-like
        Measured wave angular frequencies in rad/s.
    measured_bw_values : array-like
        Measured wave magnetic field amplitudes in nT.
    omega_ce_eq : float
        Equatorial electron cyclotron angular frequency in rad/s.
    omega_ratios : array-like
        Ratios of the requested wave angular frequencies to the
        equatorial electron cyclotron angular frequency.

    Returns
    -------
    requested_omega_values : numpy.ndarray
        Requested wave angular frequencies in rad/s.
    selected_bw_values : numpy.ndarray
        Interpolated wave magnetic field amplitudes in nT
        corresponding to the requested frequencies.

    Raises
    ------
    ValueError
        If any requested wave angular frequency falls outside
        the measured frequency range.
    """

    measured_omega_values = np.asarray(measured_omega_values)
    measured_bw_values = np.asarray(measured_bw_values)

    sorted_indices = np.argsort(measured_omega_values)
    measured_omega_values = measured_omega_values[sorted_indices]
    measured_bw_values = measured_bw_values[sorted_indices]

    requested_omega_values = np.asarray([
        ratio * omega_ce_eq for ratio in omega_ratios
    ])

    outside_range = (
            (requested_omega_values < measured_omega_values.min())
            | (requested_omega_values > measured_omega_values.max())
    )

    if np.any(outside_range):
        invalid_values = requested_omega_values[outside_range]
        raise ValueError(f"Requested ω values are outside measured FFT range: {invalid_values}")

    selected_bw_values = np.interp(
        requested_omega_values,
        measured_omega_values,
        measured_bw_values
    )

    return requested_omega_values, selected_bw_values
