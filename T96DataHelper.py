import pyspedas
import numpy as np
from pyspedas.projects.themis import fft, state
from pyspedas.geopack import get_tsy_params
from pyspedas import (
    get_data,
    join_vec,
    store_data,
    tinterpol,
    tplot
)

def load_required_data(probe, time_range, no_update):
    """Load THEMIS and OMNI data required for the program.

        Parameters
        ----------
        probe : str
            THEMIS probe identifier: "a", "b", "c", "d", or "e".
        time_range : list of str
            Start and end dates for data retrieval in YYYY-MM-DD format.
        no_update : bool
            Whether to avoid downloading updated data.

        Returns
        -------
        scm2_times : numpy.ndarray
            Timestamps for the SCM2 spectral data in Unix seconds.
        scm2_psd : numpy.ndarray
            Power spectral density values for SCM2 in nT^2/Hz.
        frequencies_hz : numpy.ndarray
            Frequency values for SCM2 in Hz.
        scm3_times : numpy.ndarray
            Timestamps for the SCM3 spectral data in Unix seconds.
        scm3_psd : numpy.ndarray
            Power spectral density values for SCM3 in nT^2/Hz.
        scm3_frequencies : numpy.ndarray
            Frequency values for SCM3 in Hz.

        Raises
        ------
        RuntimeError
            If any required THEMIS or OMNI data variables are missing
            for the selected time range.
        """
    state_variables = [f'th{probe}_pos_gsm']
    omni_variables = [
        'BX_GSE',
        'BY_GSM',
        'BZ_GSM',
        'SYM_H',
        'proton_density',
        'flow_speed',
        'Pressure'
    ]
    fft_variables = [f'th{probe}_fff_32_scm2',
                     f'th{probe}_fff_32_scm3'
                     ]

    state_data = state(
        probe=probe,
        trange=time_range,
        no_update=no_update,
        varnames=state_variables
    )

    omni_data = pyspedas.projects.omni.data(
        trange=time_range,
        datatype='1min',
        level='hro',
        time_clip=True,
        no_update=no_update,
        varnames=omni_variables
    )

    fft_data = fft(
        probe=probe,
        trange=time_range,
        no_update=no_update,
        varnames=fft_variables
    )

    if (
            len(fft_data) != len(fft_variables)
            or len(omni_data) != len(omni_variables)
            or len(state_data) != len(state_variables)
    ):
        raise RuntimeError("Required data are missing for the selected time")

    scm2_times, scm2_psd, frequencies_hz = get_data(f'th{probe}_fff_32_scm2')
    scm3_times, scm3_psd, scm3_frequencies = get_data(f'th{probe}_fff_32_scm3')

    return (
        scm2_times,
        scm2_psd,
        frequencies_hz,
        scm3_times,
        scm3_psd,
        scm3_frequencies
    )



def prepare_trace_seed(probe, selected_time_epoch, plot_enabled=False):
    """Prepare T96 parameters and get selected GSM seed position.

    Parameters
    ----------
    probe : str
        THEMIS probe identifier: "a", "b", "c", "d", or "e".
    selected_time_epoch : float
        Requested time as a Unix timestamp in seconds.
    plot_enabled : bool, optional
        Whether to display diagnostic plots comparing T96 model
        magnetic field data with satellite data. (Default: False)

    Returns
    -------
    seed_time_epoch : float
        Unix timestamp in seconds of the spacecraft position
        closest to the requested time.
    seed_parameters : numpy.ndarray
        T96 model parameters corresponding to the selected seed time.
    """

    join_vec(['BX_GSE', 'BY_GSM', 'BZ_GSM'])

    # Get necessary params for T96 model
    parameter_variable = get_tsy_params(
        'SYM_H',  # dst_tvar
        'BX_GSE-BY_GSM-BZ_GSM_joined',  # imf_tvar
        'proton_density',  # Np_tvar
        'flow_speed',  # Vp_tvar
        't96',
        pressure_tvar="Pressure",
        speed=True,
    )

    # Interpolate parameters onto the position time stream
    tinterpol(f'th{probe}_pos_gsm', parameter_variable)

    # Plot T96 plot with real satellite B field data for comparison
    if plot_enabled:
        tplot([
            f'th{probe}_pos_gsm-itrp_bt96',
            f'th{probe}_pos_gsm-itrp'
        ])

    # get arrays of time/location, params
    position_data = get_data(f'th{probe}_pos_gsm-itrp')
    parameter_data = get_data(parameter_variable)

    # Choose time
    selected_index = np.argmin(
        np.abs(position_data.times - selected_time_epoch)
    )

    # get specific time, location, params at specific point

    seed_time_epoch = position_data.times[selected_index]

    seed_xyz_gsm = np.atleast_2d(
        position_data.y[selected_index]
    )

    seed_parameters = np.atleast_2d(
        parameter_data.y[selected_index]
    )

    # store seed data as a size: 1 matrix
    store_data(
        'seed_pos_gsm',
        data={
            'x': np.array([seed_time_epoch]),
            'y': seed_xyz_gsm
        }
    )

    # Set coordinates to GSM
    pyspedas.set_coords('seed_pos_gsm', 'GSM')

    return seed_time_epoch, seed_parameters

