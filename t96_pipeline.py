"""T96 magnetic-field tracing and wave-particle calculation pipeline.

This module downloads THEMIS and OMNI data, traces a magnetic field line
using the T96 model, and calculates numerous global quantities and quantities along the field line.
"""

import time
import matplotlib.pyplot as plt
import numpy as np
from pyspedas import (
    get_data
)
from T96EquationHelper import (
    calculate_global_values,
    add_trace_equation_values,
    calculate_wave_results
)
import T96plotsHelper
from T96TraceLineHelper import trace_field_line,process_trace_results
from T96WaveHelper import get_bw_values, select_wave_values
from T96inputHelper import get_inputs
from T96DataHelper import prepare_trace_seed, load_required_data


# CONSTANTS
EARTH_RADIUS_KM = 6371.2
ELECTRON_CYCLOTRON_FREQUENCY_HZ_PER_NT = 28.0
ELECTRON_REST_ENERGY_KEV = 511.0
OMEGA_RATIOS = (1 / 2, 1 / 5, 1 / 10, 1 / 20)


def run_pipeline(
        probe,
        time_range,
        selected_time,
        energy_kev,
        theta,
        max_step=0.05,
        plot_enabled=False,
        no_update=False,
        output_csv=None
):


    """Run the complete T96 field-line calculation pipeline.

    Parameters
    ----------
    probe : str
        THEMIS probe identifier: "a", "b", "c", "d", or "e".
    time_range : sequence of str
        Start and end times used when retrieving THEMIS and OMNI data.
    selected_time : datetime.datetime
        UTC time used to select the spacecraft position and model inputs.
    energy_kev : float
        Electron kinetic energy in keV.
    theta : float
        Wave-normal angle in radians.
    max_step : float, optional
        Maximum field-line tracing step in Earth radii.
    plot_enabled : bool, optional
        Whether to display diagnostic and result plots.
    no_update : bool, optional
        Whether PySPEDAS should avoid downloading updated data.
    output_csv : str or pathlib.Path, optional
        File path where the completed results DataFrame should be saved.
        If omitted, no CSV file is created.
    Returns
    -------
    full_df : pandas.DataFrame
        Values calculated at every trace point for every selected wave
        frequency.
    global_values : dict
        Values that remain constant across the field line, including the
        L-shell, equatorial magnetic-field strength, equatorial electron
        cyclotron frequency, equatorial plasma density, and Lorentz factor.
    """

    selected_time_epoch = selected_time.timestamp()

    start_time = time.perf_counter()

    (
        scm2_times,
        scm2_psd,
        frequencies_hz,
        scm3_times,
        scm3_psd,
        scm3_frequencies
    ) = load_required_data(probe, time_range, no_update)

    data_loaded_time = time.perf_counter()

    seed_time_epoch, seed_parameters = prepare_trace_seed(
        probe=probe,
        selected_time_epoch=selected_time_epoch,
        plot_enabled=plot_enabled
    )

    _, measured_bw_values, measured_omega_values = get_bw_values(
        probe,
        seed_time_epoch,
        time_range,
        no_update,
        plot_enabled,
        scm2_times,
        scm2_psd,
        frequencies_hz,
        scm3_times,
        scm3_psd,
        scm3_frequencies
    )

    if plot_enabled:
        T96plotsHelper.omegaVsBws(
            measured_omega_values,
            measured_bw_values,
            0,
            25000
        )

    # Compute Trace
    (
        s_full,
        r_full_gsm,
        b_full,
        bmag_full
    ) = trace_field_line(
        seed_parameters,
        max_step,
        plot_enabled
    )

    df, seed_time = process_trace_results(
        s_full=s_full,
        r_full_gsm=r_full_gsm,
        b_full=b_full,
        bmag_full=bmag_full,
        plot_enabled=plot_enabled
    )
    #print("time: ", seed_time)
    # GLOBAL VALS

    l_shell = get_data('tha_pos_lshell').y[0]
    equatorial_bmag = df.loc[np.isclose(df['s'], 0.0), 'Bmag'].iloc[0]

    f_ce_eq, omega_ce_eq, n_pe_eq, gamma = calculate_global_values(
        l_shell=l_shell,
        equatorial_bmag=equatorial_bmag,
        energy_kev=energy_kev
    )

    global_values = {
        'l_shell': l_shell,
        'B_0': equatorial_bmag,
        'f_ce_eq': f_ce_eq,
        'omega_ce_eq': omega_ce_eq,
        'n_pe_eq': n_pe_eq,
        'gamma': gamma,
    }

    omega_ratios = OMEGA_RATIOS

    omega_values, bw_values = select_wave_values(
        measured_omega_values=measured_omega_values,
        measured_bw_values=measured_bw_values,
        omega_ce_eq=omega_ce_eq,
        omega_ratios=omega_ratios
    )

    df = add_trace_equation_values(
        df=df,
        n_pe_eq=n_pe_eq
    )

    full_df = calculate_wave_results(
        df=df,
        omega_values=omega_values,
        bw_values=bw_values,
        theta=theta,
        gamma=gamma,
        omega_ce_eq=omega_ce_eq
    )

    if plot_enabled:
        plt.figure()
        plt.scatter(full_df['omega'], full_df['S'])
        plt.xlabel('ω [rad/s]')
        plt.ylabel('S')
        plt.title(fr'S vs $\omega$')
        plt.grid(False)
        plt.tight_layout()

    if plot_enabled:
        plt.show()

    if output_csv is not None:
        full_df.to_csv(output_csv, index=False)


    end_time = time.perf_counter()

    print(f'DATA RETRIEVAL/DOWNLOAD TIME: {(data_loaded_time - start_time):.2f}')

    print(f'PROGRAM TIME: {(end_time - data_loaded_time):.2f}')

    if plot_enabled:
        T96plotsHelper.plotNewOmegaAlphaEQ(full_df, omega_values, omega_ratios, energy_kev, seed_time)

    print(full_df["S"][500])



    return full_df, global_values


def main():
    """Collect predefined inputs and run the pipeline.

    Returns
    -------
    full_df : pandas.DataFrame
        Values calculated at every field line position for each
        selected wave frequency.
    global_values : dict
        Values that remain constant across the field line, including
        the L-shell, equatorial magnetic field magnitude, equatorial
        electron cyclotron frequency, equatorial electron number
        density, and Lorentz factor.
    """


    (
        probe,
        time_range,
        selected_time,
        energy_kev,
        theta,
        max_step,
        output_csv
    ) = get_inputs(False)


    return run_pipeline(
        probe=probe,
        time_range=time_range,
        selected_time=selected_time,
        energy_kev=energy_kev,
        theta=theta,
        max_step=max_step,
        output_csv = output_csv
    )

if __name__ == "__main__":
    main()
