import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import pandas as pd
import spacepy.time as spt
import spacepy.coordinates as spc

from pyspedas.geopack import ttrace2endpoint
from pyspedas import (
    get_data,
    tplotxy3,
    calculate_lshell,
    set_units,
    cotrans
)
from T96plotsHelper import plotVectorDistances

EARTH_RADIUS_KM = 6371.2

def get_trace_distance_and_B(pos_var, b_var):
    """Calculate cumulative distance and magnetic field magnitude along a field line.

    Parameters
    ----------
    pos_var : str
        Name of the PySPEDAS variable containing field line positions.
    b_var : str
        Name of the PySPEDAS variable containing magnetic field vectors.

    Returns
    -------
    s_km : numpy.ndarray
        Cumulative distance along the field line from the seed point in km.
    R : numpy.ndarray
        Cartesian positions along the field line in km.
    B : numpy.ndarray
        Magnetic field vectors along the field line in nT.
    Bmag : numpy.ndarray
        Magnetic field magnitudes along the field line in nT.

    Raises
    ------
    ValueError
        If the number of position points does not match the number
        of magnetic field vectors.
    """

    # Load position and B data
    R = np.array(get_data(pos_var).y)

    B = np.array(get_data(b_var).y)

    # Remove leading time dimension: (1, n, 3) -> (n, 3)
    if R.ndim == 3:
        R = R[0]
    if B.ndim == 3:
        B = B[0]

    if R.shape[0] != B.shape[0]:
        raise ValueError(f"Mismatch: {pos_var} has {R.shape[0]} points, {b_var} has {B.shape[0]} points")

    # Cumulative distance from seed
    dR = np.diff(R, axis=0)
    ds = np.sqrt(np.sum(dR**2, axis=1))
    s_km = np.concatenate(([0.0], np.cumsum(ds)))

    # B magnitude
    Bmag = np.sqrt(B[:, 0]**2 + B[:, 1]**2 + B[:, 2]**2)

    return s_km, R, B, Bmag


def plot_full_field_line_equator_zero(pos_south, b_south, pos_north, b_north, pos_eq, pl=False, title='Full field line'):
    """Combine northern and southern field line traces with the equator at s = 0.

    Parameters
    ----------
    pos_south : str
        Name of the PySPEDAS variable containing southern trace positions.
    b_south : str
        Name of the PySPEDAS variable containing southern trace magnetic field vectors.
    pos_north : str
        Name of the PySPEDAS variable containing northern trace positions.
    b_north : str
        Name of the PySPEDAS variable containing northern trace magnetic field vectors.
    pos_eq : str
        Name of the PySPEDAS variable containing the equatorial trace positions.
    pl : bool, optional
        Whether to plot magnetic field components along the field line
        (default is False).
    title : str, optional
        Title of the magnetic field plot (default is 'Full field line').

    Returns
    -------
    s_full : numpy.ndarray
        Distances along the combined field line in km, with s = 0
        at the sampled point closest to the magnetic equator.
    R_full : numpy.ndarray
        Cartesian positions along the combined field line in km.
    B_full : numpy.ndarray
        Magnetic field vectors along the combined field line in nT.
    Bmag_full : numpy.ndarray
        Magnetic field magnitudes along the combined field line in nT.
    """

    # Load traces
    _, R_south, B_south_arr, Bmag_south = get_trace_distance_and_B(pos_south, b_south)
    _, R_north, B_north_arr, Bmag_north = get_trace_distance_and_B(pos_north, b_north)

    # Combine full trace: south endpoint -> seed -> north endpoint
    R_full = np.concatenate((R_south[::-1], R_north[1:]), axis=0)

    # Combine B values in same order
    B_full = np.concatenate((B_south_arr[::-1], B_north_arr[1:]), axis=0)
    Bmag_full = np.concatenate((Bmag_south[::-1], Bmag_north[1:]))

    # Calculate cumulative distance along full trace from south endpoint
    dR = np.diff(R_full, axis=0)
    ds = np.sqrt(np.sum(dR**2, axis=1))
    s_from_south = np.concatenate(([0.0], np.cumsum(ds)))

    # Load actual equator trace point
    R_eq = np.array(get_data(pos_eq).y)
    if R_eq.ndim == 3:
        R_eq = R_eq[0]

    R_eq_point = R_eq[-1]

    # Find closest point on full trace to actual equator point
    idx_eq = np.argmin(np.linalg.norm(R_full - R_eq_point, axis=1))

    # Shift so equator is s = 0
    s_full = s_from_south - s_from_south[idx_eq]

    print("Equator index:", idx_eq)
    print("s at equator:", s_full[idx_eq])

    print(f'IS TRUE: {R_south[0] == R_north[0]}')
    if pl:
        plt.figure()
        plt.plot(s_full, B_full[:, 0], label='Bx', color='blue')
        plt.plot(s_full, B_full[:, 1], label='By', color='green')
        plt.plot(s_full, B_full[:, 2], label='Bz', color='red')
        plt.xlabel('Distance along field line from equator [km]')
        plt.ylabel('Magnetic field [nT]')
        plt.title(title)
        plt.legend()
        plt.grid(False)

        plt.show()

    return s_full, R_full, B_full, Bmag_full


def trace_field_line(seed_parameters, max_step, plot_enabled=False):
    """Trace a magnetic field line using the T96 model.

       Traces from the seed position to the magnetic equator and both
       ionospheres, then combines the traces with the equator at s = 0.

       Parameters
       ----------
       seed_parameters : numpy.ndarray
           T96 model parameters for the selected time and seed position.
       max_step : float
           Maximum field-line tracing step size in Earth radii.
       plot_enabled : bool, option
           Whether to display the field line traces and magnetic field plots. (Default = False)

       Returns
       -------
       s_full : numpy.ndarray
           Distances along the combined field line in km, with s = 0
           at the sampled point closest to the magnetic equator.
       R_full : numpy.ndarray
           Cartesian positions along the combined field line in km.
       B_full : numpy.ndarray
           Magnetic field vectors along the combined field line in nT.
       Bmag_full : numpy.ndarray
           Magnetic field magnitudes along the combined field line in nT.
       """

    # Trace seed position to the equator and both south and north ionosphere

    ttrace2endpoint(
        'seed_pos_gsm',
        't96',
        'equator',
        trace_name='seed_trace_t96',
        bvec_name='seed_trace_b_t96',
        foot_name='seed_equator_foot_t96',
        diag_reached_name='seed_trace_reached_t96',
        diag_s_max_name='seed_trace_smax_t96',
        diag_npts_name='seed_trace_npts_t96',
        parmod=seed_parameters,
        km=True,
        max_step=max_step,
    )

    ttrace2endpoint(
        'seed_pos_gsm',
        't96',
        'ionosphere-north',
        trace_name='seed_trace_north_t96',
        bvec_name='seed_trace_b_north_t96',
        diag_s_max_name='seed_trace_smax_north_t96',
        parmod=seed_parameters,
        km=True,
        max_step=max_step,

    )

    ttrace2endpoint(
        'seed_pos_gsm',
        't96',
        'ionosphere-south',
        trace_name='seed_trace_south_t96',
        bvec_name='seed_trace_b_south_t96',
        diag_s_max_name='seed_trace_smax_south_t96',
        parmod=seed_parameters,
        km=True,
        max_step=max_step,
    )

    # Plot all three traces on same plots
    if plot_enabled:
        tplotxy3(
            ['seed_trace_north_t96', 'seed_trace_south_t96', 'seed_trace_t96'],
            legend_names=['North trace', 'South trace', 'Equator trace'],
            show_centerbody=True
        )

    return plot_full_field_line_equator_zero(
        'seed_trace_south_t96',
        'seed_trace_b_south_t96',
        'seed_trace_north_t96',
        'seed_trace_b_north_t96',
        'seed_trace_t96',
        plot_enabled,
        title='T96 field line from seed point'
    )

def convert_gsm_to_mag(gsm_positions, timestamps):
    """Convert GSM Cartesian positions to MAG Cartesian coordinates and magnetic latitude.

       Parameters
       ----------
       gsm_positions : array-like
           Cartesian positions in GSM coordinates in km.
       timestamps : float or array-like
           Unix timestamps corresponding to the GSM positions in seconds.

       Returns
       -------
       mag_positions : numpy.ndarray
           Cartesian positions in MAG coordinates in km.
       magnetic_latitudes : numpy.ndarray
           Magnetic latitudes of the positions in degrees.
       """

    gsm_positions = np.atleast_2d(gsm_positions)
    timestamps = np.atleast_1d(timestamps)

    mag_positions = cotrans(
        time_in=timestamps,
        data_in=gsm_positions,
        coord_in='gsm',
        coord_out='mag'
    )

    datetimes = [
        dt.datetime.fromtimestamp(timestamp, dt.UTC) for timestamp in timestamps
    ]

    ticks = spt.Ticktock(datetimes, 'UTC')

    mag_coordinates = spc.Coords(
        mag_positions,
        'MAG',
        'car'
    )

    mag_coordinates.ticks = ticks

    # MAG Cartesian -> MAG spherical
    mag_spherical = mag_coordinates.convert(
        'MAG',
        'sph'
    )

    # Get radius, MLAT, MLONG

    magnetic_latitudes = mag_spherical.data[:, 1]

    return mag_positions, magnetic_latitudes


def build_trace_dataframe(
        s_full,
        r_full_gsm,
        b_full,
        bmag_full,
        r_full_mag,
        trace_magnetic_latitudes,
        equator_xyz_gsm,
        equator_b,
        equator_xyz_mag,
        equator_mlat
):
    """Build field line DataFrame and convert distances.

    Parameters
    ----------
    s_full : numpy.ndarray
        Distances along the field line in km, with the equator at s = 0.
    r_full_gsm : numpy.ndarray
        Cartesian positions along the field line in GSM coordinates in km.
    b_full : numpy.ndarray
        Magnetic field vectors along the field line in nT.
    bmag_full : numpy.ndarray
        Magnetic field magnitudes along the field line in nT.
    r_full_mag : numpy.ndarray
        Cartesian positions along the field line in MAG coordinates in km.
    trace_magnetic_latitudes : numpy.ndarray
        Magnetic latitudes along the field line in degrees.
    equator_xyz_gsm : array-like
        Equatorial position in GSM coordinates in km.
    equator_b : array-like
        Magnetic field vector at the equator in nT.
    equator_xyz_mag : array-like
        Equatorial position in MAG coordinates in km.
    equator_mlat : float
        Magnetic latitude of the equatorial position in degrees.

    Returns
    -------
    df : pandas.DataFrame
        Field line data sorted by distance from the equator, containing
        distances (s) in meters, GSM and MAG positions in Earth radii,
        magnetic latitude (MLAT) in degrees, and magnetic field
        components and magnitudes in nT.
    """

    df = pd.DataFrame({
        's': s_full,
        'MLAT': trace_magnetic_latitudes,
        'x_gsm': r_full_gsm[:, 0],
        'y_gsm': r_full_gsm[:, 1],
        'z_gsm': r_full_gsm[:, 2],
        'Bx': b_full[:, 0],
        'By': b_full[:, 1],
        'Bz': b_full[:, 2],
        'Bmag': bmag_full,
        'x_mag': r_full_mag[:, 0],
        'y_mag': r_full_mag[:, 1],
        'z_mag': r_full_mag[:, 2],
    })

    equator_xyz_gsm = np.asarray(equator_xyz_gsm).squeeze()

    equator_b = np.asarray(equator_b).squeeze()

    equator_xyz_mag = np.asarray(equator_xyz_mag).squeeze()

    row_s0 = {
        's': 0.0,
        'MLAT': equator_mlat,
        'x_gsm': equator_xyz_gsm[0],
        'y_gsm': equator_xyz_gsm[1],
        'z_gsm': equator_xyz_gsm[2],
        'Bx': equator_b[0],
        'By': equator_b[1],
        'Bz': equator_b[2],
        'Bmag': np.linalg.norm(equator_b),
        'x_mag': equator_xyz_mag[0],
        'y_mag': equator_xyz_mag[1],
        'z_mag': equator_xyz_mag[2],
    }

    df = df[~np.isclose(df['s'], 0.0)].copy()
    df = pd.concat(
        [df, pd.DataFrame([row_s0])],
        ignore_index=True
    )

    # Sort by s
    df = df.sort_values('s').reset_index(drop=True)

    # KM TO M
    df['s'] = df['s'] * 1000

    # KM TO RE
    coordinate_columns = [
        'x_gsm',
        'y_gsm',
        'z_gsm',
        'x_mag',
        'y_mag',
        'z_mag'
    ]

    df[coordinate_columns] /= EARTH_RADIUS_KM
    return df



def process_trace_results(
        s_full,
        r_full_gsm,
        b_full,
        bmag_full,
        plot_enabled=False
):
    """Convert trace coordinates and build field line DataFrame.

    Parameters
    ----------
    s_full : numpy.ndarray
        Distances along the field line in km, with the equator at s = 0.
    r_full_gsm : numpy.ndarray
        Cartesian positions along the field line in GSM coordinates in km.
    b_full : numpy.ndarray
        Magnetic field vectors along the field line in nT.
    bmag_full : numpy.ndarray
        Magnetic field magnitudes along the field line in nT.
    plot_enabled : bool, optional
        Whether to display diagnostic plots (default is False).

    Returns
    -------
    trace_df : pandas.DataFrame
        Field line data containing distances in meters, GSM and MAG
        positions in Earth radii, magnetic latitudes in degrees,
        and magnetic field values in nT.
    equator_time : datetime.datetime
        UTC time associated with the equatorial trace position.
    """

    # Set seed position to km
    set_units('seed_pos_gsm', 'km')

    # Calculate l-shell using seed position/time
    calculate_lshell('seed_pos_gsm', 'tha_pos_lshell')

    # Access the equatorial data for the selected satellite point

    equator_data = get_data('seed_equator_foot_t96')

    # Access coordinates in shape = [x,y,z]
    equator_xyz_gsm = np.array(equator_data.y[0])

    # Access time of data
    equator_time_epoch = equator_data.times[0]

    equator_time = dt.datetime.fromtimestamp(
        equator_time_epoch,
        dt.UTC
    )

    equator_xyz_mag, equator_magnetic_latitudes = convert_gsm_to_mag(
        equator_xyz_gsm,
        equator_time_epoch
    )

    equator_mlat = equator_magnetic_latitudes[0]

    # Get MLAT for each trace position

    if plot_enabled:
        plotVectorDistances(
            s_full,
            r_full_gsm
        )

    trace_times = np.full(
        r_full_gsm.shape[0],
        equator_time_epoch
    )

    r_full_mag, trace_magnetic_latitudes = convert_gsm_to_mag(
        r_full_gsm,
        trace_times
    )

    if plot_enabled:
        plt.figure()
        plt.plot(s_full, trace_magnetic_latitudes, label='MLAT [°]', color='red')
        plt.grid(False)
        plt.xlabel('Distance along field line from equator [km]')
        plt.ylabel('MLAT [°]')
        plt.title("Distance along field line from equator Vs. MLAT")
        plt.show()

    equator_b = np.asarray(
        get_data('seed_trace_b_t96').y
    ).squeeze()[-1]

    trace_df = build_trace_dataframe(
        s_full=s_full,
        r_full_gsm=r_full_gsm,
        b_full=b_full,
        bmag_full=bmag_full,
        r_full_mag=r_full_mag,
        trace_magnetic_latitudes=trace_magnetic_latitudes,
        equator_xyz_gsm=equator_xyz_gsm,
        equator_b=equator_b,
        equator_xyz_mag=equator_xyz_mag,
        equator_mlat=equator_mlat,
    )

    return trace_df, equator_time




