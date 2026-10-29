
from pyspedas import  get_data, tplot
from pyspedas.projects.themis import fft
import numpy as np




def get_bw_values(probe, seed_time_epoch, time_range, no_update, plot_enabled, scm2_times, scm2_psd, frequencies_hz, scm3_times, scm3_psd, scm3_frequencies):
    """Calculate wave magnetic field amplitudes from THEMIS spectral data.

    Parameters
    ----------
    probe : str
        THEMIS probe identifier: "a", "b", "c", "d", or "e".
    seed_time_epoch : float
        Selected time as a Unix timestamp in seconds.
    time_range : list of str
        Start and end dates for data retrieval.
    no_update : bool
        Whether to avoid downloading updated data.
    plot_enabled : bool
        Whether to display the wave spectral data plot.
    scm2_times : array-like
        Timestamps for the first magnetic field spectral dataset.
    scm2_psd : array-like
        Power spectral density values for the first dataset in nT^2/Hz.
    frequencies_hz : array-like
        Center frequencies of the spectral bins in Hz.
    scm3_times : array-like
        Timestamps for the second magnetic field spectral dataset.
    scm3_psd : array-like
        Power spectral density values for the second dataset in nT^2/Hz.
    scm3_frequencies : array-like
        Frequency values for the second spectral dataset in Hz.

    Returns
    -------
       frequencies_hz : numpy.ndarray
           Center frequencies of the spectral bins in Hz.
       measured_bw_values : numpy.ndarray
           Wave magnetic field amplitudes in nT.
       measured_omega_values : numpy.ndarray
           Wave angular frequencies in rad/s.
    """


    print(len(scm2_times), len(scm3_times), len(scm3_psd), len(scm3_frequencies), len(scm3_psd[0]))

    print(scm2_times)

    idx_timeW1 = np.argmin(np.abs(np.array(scm2_times) - seed_time_epoch))

    print(idx_timeW1)
    print("Closest time:", scm2_times[idx_timeW1], seed_time_epoch, (scm2_times[idx_timeW1] - seed_time_epoch))
    print(scm3_psd[idx_timeW1])
    print(len(scm3_psd[idx_timeW1]))

    dataW = np.array(scm2_psd[idx_timeW1]) + np.array(scm3_psd[idx_timeW1])


    print(dataW)

    print(frequencies_hz, scm3_frequencies)


    # freqs = the 32 frequency centers, like [4, 20, 36, ..., 3836]
    frequencies_hz = np.array(frequencies_hz, dtype=float)

    # make 33 edges for 32 bins
    edges = np.zeros(len(frequencies_hz) + 1)

    # first bin starts at 0 Hz
    edges[0] = 0

    # middle edges are halfway between neighboring frequency centers
    edges[1:-1] = 0.5 * (frequencies_hz[:-1] + frequencies_hz[1:])

    # last bin extends past the final center by half the last spacing
    edges[-1] = frequencies_hz[-1] + 0.5 * (frequencies_hz[-1] - frequencies_hz[-2])

    # bin sizes
    bin_widths = np.diff(edges)

    print(bin_widths)

    # convert PSD to Bw
    measured_bw_values = np.sqrt(dataW * bin_widths)


    measured_omega_values = frequencies_hz * np.pi * 2

    if plot_enabled: tplot(f'th{probe}_fff_32_scm2')

    return frequencies_hz, measured_bw_values, measured_omega_values


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
