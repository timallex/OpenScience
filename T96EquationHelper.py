import numpy as np
import pandas as pd

ELECTRON_CYCLOTRON_FREQUENCY_HZ_PER_NT = 28.0
ELECTRON_REST_ENERGY_KEV = 511.0

def calculate_global_values(l_shell, equatorial_bmag, energy_kev):
    """Calculate values that remain constant along every point of a field line.

    Parameters
    ----------
    l_shell : float
        L-shell value of the magnetic field line.
    equatorial_bmag : float
        Magnetic field magnitude at the equator in nT.
    energy_kev : float
        Electron kinetic energy in keV.

    Returns
    -------
    f_ce_eq : float
        Equatorial electron cyclotron frequency in Hz.
    omega_ce_eq : float
        Equatorial electron cyclotron angular frequency in rad/s.
    n_pe_eq : float
        Equatorial electron number density in cm^-3.
    gamma : float
        Electron Lorentz factor.
    """


    f_ce_eq = ELECTRON_CYCLOTRON_FREQUENCY_HZ_PER_NT * equatorial_bmag
    omega_ce_eq = 2 * np.pi * f_ce_eq

    n_pe_eq = 100 / ((l_shell / 3) ** 4)

    gamma = (energy_kev / ELECTRON_REST_ENERGY_KEV) + 1

    return f_ce_eq, omega_ce_eq, n_pe_eq, gamma


def add_trace_equation_values(df, n_pe_eq):
    """Calculate values that depend on field line position but not omega or Bw.

    Parameters
    ----------
    df : pandas.DataFrame
        Field line data containing magnetic field magnitude (Bmag),
        distance along the field line (s), and magnetic latitude (MLAT).
    n_pe_eq : float
        Equatorial electron number density in cm^-3.

    Returns
    -------
    df : pandas.DataFrame
        Copy of the input DataFrame with additional columns for electron
        cyclotron frequency (f_ce), electron cyclotron angular frequency
        (omega_ce), its gradient along the field line (dOce_ds),
        electron plasma frequency (f_pe), and electron plasma angular
        frequency (omega_pe).
    """

    df = df.copy()

    df['f_ce'] = getFce(df['Bmag'])
    df['omega_ce'] = getOmegaCE(df['f_ce'])
    df['dOce_ds'] = getd0ce_ds(df['omega_ce'], df['s'])
    df['f_pe'] = getFpe(n_pe_eq, df['MLAT'])
    df['omega_pe'] = getOmegape(df['f_pe'])

    return df


def calculate_wave_results(
        df,
        omega_values,
        bw_values,
        theta,
        gamma,
        omega_ce_eq
):
    """Calculate values that depend on Bw and/or omega.

    Parameters
    ----------
    df : pandas.DataFrame
        Field line data containing previously calculated plasma properties,
        including omega_pe, omega_ce, and dOce_ds.
    omega_values : array-like of float
        Wave angular frequencies in rad/s.
    bw_values : array-like of float
        Wave magnetic field amplitudes in nT corresponding to omega_values.
    theta : float
        Wave-normal angle in radians.
    gamma : float
        Electron Lorentz factor.
    omega_ce_eq : float
        Equatorial electron cyclotron angular frequency in rad/s.

    Returns
    -------
    full_df : pandas.DataFrame
        DataFrame containing calculated wave-particle interaction values
        for every field line position and wave frequency, including
        wave number (k), refractive index (N), pitch angles (alpha,
        alpha_eq), velocities (v_perp, v_R), and S factor.

    Raises
    ------
    ValueError
        If no omega and Bw pairs are provided, or if the number of
        omega values does not match the number of Bw values.
    """
    omega_values = np.asarray(omega_values)
    bw_values = np.asarray(bw_values)

    if omega_values.size == 0:
        raise ValueError("At least one ω and Bw pair is required")

    if omega_values.size != bw_values.size:
        raise ValueError(
            "The number of ω values must match the number of Bw values"
        )

    result_dataframes = []

    for omega, bw in zip(omega_values, bw_values):
        wave_df = df.copy()

        wave_df['omega'] = omega
        wave_df['B_w'] = bw
        wave_df['omega_w'] = 2 * np.pi * ELECTRON_CYCLOTRON_FREQUENCY_HZ_PER_NT * wave_df['B_w']

        wave_df['k'] = getK(
            wave_df['omega_pe'],
            wave_df['omega_ce'],
            theta,
            wave_df['omega']
        )

        wave_df['N'] = getN(
            wave_df['k'],
            wave_df['omega']
        )

        wave_df['cos_alpha'] = getAlphaWithResCond(
            wave_df['omega'],
            gamma,
            wave_df['omega_ce'],
            wave_df['k'],
            theta
        )

        valid_alpha = (wave_df['cos_alpha'] > 0) & (wave_df['cos_alpha'] < 1)

        wave_df['alpha'] = np.nan

        wave_df.loc[valid_alpha, 'alpha'] = np.degrees(
            np.arccos(wave_df.loc[valid_alpha, 'cos_alpha'])
        )

        wave_df['alpha_eq'] = getAlphaEq(
            wave_df['alpha'],
            wave_df['omega_ce'],
            omega_ce_eq
        )

        wave_df['v_perp'] = getVperp(
            gamma,
            wave_df['alpha']
        )

        wave_df['v_R'] = getVr(
            wave_df['k'],
            wave_df['omega'],
            wave_df['omega_ce'],
            gamma
        )

        wave_df['S'] = getS(
            wave_df['N'],
            wave_df['k'],
            wave_df['v_perp'],
            wave_df['omega_w'],
            gamma,
            wave_df['omega_ce'],
            wave_df['v_R'],
            wave_df['omega'],
            wave_df['dOce_ds']
        )

        result_dataframes.append(wave_df)

    full_df = pd.concat(result_dataframes, ignore_index=True)

    return full_df



def getFce(Bmag):
    """Calculate the electron cyclotron frequency.

    Parameters
    ----------
    Bmag : float or array-like
        Magnetic field magnitude in nT.

    Returns
    -------
    float or array-like
        Electron cyclotron frequency in Hz.
    """

    return 28 * Bmag

def getOmegaCE(Fce):
    """Calculate the electron cyclotron angular frequency.

    Parameters
    ----------
    Fce : float or array-like
        Electron cyclotron frequency in Hz.

    Returns
    -------
    float or array-like
        Electron cyclotron angular frequency in rad/s.
    """

    return 2 * np.pi * Fce

def getFpe(n_pe_eq, MLAT):
    """Calculate the electron plasma frequency along a field line.

    Parameters
    ----------
    n_pe_eq : float
        Equatorial electron number density in cm^-3.
    MLAT : float or array-like
        Magnetic latitude in degrees.

    Returns
    -------
    float or array-like
        Electron plasma frequency in Hz.
    """

    return 8980 * np.sqrt(n_pe_eq) * ((np.cos(np.radians(MLAT)))**(-5/2))


def getOmegape(f_pe):
    """Calculate the electron plasma angular frequency.

    Parameters
    ----------
    f_pe : float or array-like
        Electron plasma frequency in Hz.

    Returns
    -------
    float or array-like
        Electron plasma angular frequency in rad/s.
    """

    return 2 * np.pi * f_pe

def getK(omega_pe, omega_ce, theta, omega, c=2.99792458e8):
    """Calculate the wave number using the plasma dispersion relation.

    Parameters
    ----------
    omega_pe : float or array-like
        Electron plasma angular frequency in rad/s.
    omega_ce : float or array-like
        Electron cyclotron angular frequency in rad/s.
    theta : float
        Wave-normal angle in radians.
    omega : float or array-like
        Wave angular frequency in rad/s.
    c : float, optional
        Speed of light in m/s (default is 2.99792458e8).

    Returns
    -------
    float or array-like
        Wave number in rad/m.
    """

    return ( omega_pe / c) * ((((omega_ce * np.cos(theta))/omega) - 1)**(-1/2))

def getN(k, omega, c=2.99792458e8):
    """Calculate the wave refractive index.

    Parameters
    ----------
    k : float or array-like
        Wave number in rad/m.
    omega : float or array-like
        Wave angular frequency in rad/s.
    c : float, optional
        Speed of light in m/s (default is 2.99792458e8).

    Returns
    -------
    float or array-like
        Wave refractive index.
    """

    return (k * c) / omega

def getAlpha(alpha_eq, omega_ce, omega_ce_eq):
    """Calculate the local electron pitch angle along a field line.

    Parameters
    ----------
    alpha_eq : float or array-like
        Equatorial electron pitch angle in degrees.
    omega_ce : float or array-like
        Local electron cyclotron angular frequency in rad/s.
    omega_ce_eq : float
        Equatorial electron cyclotron angular frequency in rad/s.

    Returns
    -------
    float or array-like
        Local electron pitch angle in degrees.
    """

    return np.degrees(np.arcsin( np.sin(np.radians(alpha_eq)) * np.sqrt(omega_ce/omega_ce_eq) ))

def getAlphaWithResCond(omega, gamma, omega_ce, k, theta, c=2.99792458e8):
    """Calculate the cosine of the electron pitch angle using the resonance condition.

        Parameters
        ----------
        omega : float or array-like
            Wave angular frequency in rad/s.
        gamma : float
            Electron Lorentz factor.
        omega_ce : float or array-like
            Electron cyclotron angular frequency in rad/s.
        k : float or array-like
            Wave number in rad/m.
        theta : float
            Wave-normal angle in radians.
        c : float, optional
            Speed of light in m/s (default is 2.99792458e8).

        Returns
        -------
        float or array-like
            Cosine of the electron pitch angle satisfying the resonance condition.
        """

    return -1 * ((omega * gamma) - omega_ce) / (c * k * np.cos(theta) * np.sqrt((gamma**2)-1))



def getResCheck(omega, gamma, k, theta, alpha, omega_ce, c=2.99792458e8):
    """Calculate the residual of the electron resonance condition.

        Parameters
        ----------
        omega : float or array-like
            Wave angular frequency in rad/s.
        gamma : float
            Electron Lorentz factor.
        k : float or array-like
            Wave number in rad/m.
        theta : float
            Wave-normal angle in radians.
        alpha : float or array-like
            Local electron pitch angle in degrees.
        omega_ce : float or array-like
            Electron cyclotron angular frequency in rad/s.
        c : float, optional
            Speed of light in m/s (default is 2.99792458e8).

        Returns
        -------
        float or array-like
            Resonance condition residual in rad/s, where a value of zero
            indicates that the resonance condition is satisfied.
        """

    return (( (omega * gamma)
                     + (c * k * np.cos(theta) * np.sqrt((gamma ** 2)-1) * np.cos(np.radians(alpha))) )
                   - omega_ce)

def getAlphaEq(alpha_deg, omega_ce,  omega_ce_eq):
    """Calculate the equatorial electron pitch angle.

        Parameters
        ----------
        alpha_deg : float or array-like
            Local electron pitch angle in degrees.
        omega_ce : float or array-like
            Local electron cyclotron angular frequency in rad/s.
        omega_ce_eq : float
            Equatorial electron cyclotron angular frequency in rad/s.

        Returns
        -------
        float or array-like
            Equatorial electron pitch angle in degrees.
        """
    return np.degrees(np.arcsin((np.sin(np.radians(alpha_deg)))/(np.sqrt(omega_ce/omega_ce_eq))))

def getVperp(gamma,alpha,c=2.99792458e8):
    """Calculate the electron velocity perpendicular to the magnetic field.

    Parameters
    ----------
    gamma : float
        Electron Lorentz factor.
    alpha : float or array-like
        Local electron pitch angle in degrees.
    c : float, optional
        Speed of light in m/s (default is 2.99792458e8).

    Returns
    -------
    float or array-like
        Perpendicular electron velocity in m/s.
    """

    return c * np.sqrt(1 - (gamma ** (-2))) * np.sin(np.radians(alpha))

def getVr(k, omega, omega_ce, gamma):
    """Calculate the electron resonance velocity.

        Parameters
        ----------
        k : float or array-like
            Wave number in rad/m.
        omega : float or array-like
            Wave angular frequency in rad/s.
        omega_ce : float or array-like
            Electron cyclotron angular frequency in rad/s.
        gamma : float
            Electron Lorentz factor.

        Returns
        -------
        float or array-like
            Electron resonance velocity in m/s.
        """
    return (1 / k) * (omega - (omega_ce / gamma))

def getS(N, k, v_perp, omega_w, gamma, omega_ce, v_R, omega, d0ce_ds):
    """Calculate the S factor for electron-wave interactions.

        Parameters
        ----------
        N : float or array-like
            Wave refractive index.
        k : float or array-like
            Wave number in rad/m.
        v_perp : float or array-like
            Electron velocity perpendicular to the magnetic field in m/s.
        omega_w : float or array-like
            Angular frequency associated with the wave magnetic field in rad/s.
        gamma : float
            Electron Lorentz factor.
        omega_ce : float or array-like
            Local electron cyclotron angular frequency in rad/s.
        v_R : float or array-like
            Electron resonance velocity in m/s.
        omega : float or array-like
            Wave angular frequency in rad/s.
        d0ce_ds : float or array-like
            Gradient of the electron cyclotron angular frequency along
            the field line in rad/(s*m).

        Returns
        -------
        float or array-like
            S factor describing electron-wave interactions.
        """

    return (((1 / (1- (N)**(-2))) * (1 / (2 * k * v_perp * omega_w)))
           * (((gamma * k * (v_perp)**2) / (omega_ce))
              - (v_R * (2 + (1 - (N)**(-2)) * ((omega_ce - (gamma * omega)) / (omega_ce - omega)))))
           * d0ce_ds
           )

def getd0ce_ds(omega_ce,s):
    """Calculate the gradient of electron cyclotron angular frequency along a field line.

    Parameters
    ----------
    omega_ce : array-like
        Electron cyclotron angular frequencies in rad/s.
    s : array-like
        Distances along the field line in meters.

    Returns
    -------
    numpy.ndarray
        Gradient of electron cyclotron angular frequency along
        the field line in rad/(s*m).
    """
    return np.gradient(omega_ce, s)