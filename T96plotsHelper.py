import matplotlib.pyplot as plt
from fractions import Fraction



def omegaVsBws(omegas, Bws, min = 10000, max = 25000):
    """Plot wave magnetic field amplitude against wave angular frequency.

    Parameters
    ----------
    omegas : array-like
        Wave angular frequencies in rad/s.
    Bws : array-like
        Wave magnetic field amplitudes in nT.
    min : float, optional
        Minimum wave angular frequency displayed in rad/s
        (default is 10000).
    max : float, optional
        Maximum wave angular frequency displayed in rad/s
        (default is 25000).
    """

    plt.figure()
    plt.plot(omegas, Bws)
    plt.xlabel('ω [rad/s]')
    plt.ylabel(r'$B_w$ [nT]')
    plt.title(r'$B_w$ vs $\omega$')
    plt.grid(False)
    plt.xlim(min, max)
    plt.show()


def plotVectorDistances(s_full,R_full_gsm):
    """Plot GSM Cartesian position components against distance along a field line.

    Parameters
    ----------
    s_full : numpy.ndarray
        Distances along the field line in km, with the equator at s = 0.
    R_full_gsm : numpy.ndarray
        Cartesian positions along the field line in GSM coordinates in km.
    """

    plt.figure()
    plt.plot(s_full, R_full_gsm[:, 0], label='x', color='blue')
    plt.plot(s_full, R_full_gsm[:, 1], label='y', color='green')
    plt.plot(s_full, R_full_gsm[:, 2], label='z', color='red')
    plt.legend()
    plt.grid(False)
    plt.show()



def plotNewOmegaAlphaEQ(df,omegas,ratios,E,seed):
    """Plot equatorial electron pitch angle against magnetic latitude.

    Parameters
    ----------
    df : pandas.DataFrame
        Field line data containing wave angular frequencies (omega),
        magnetic latitudes (MLAT), and equatorial pitch angles (alpha_eq).
    omegas : array-like
        Wave angular frequencies to plot in rad/s.
    ratios : array-like
        Ratios of the wave angular frequencies to the equatorial
        electron cyclotron angular frequency.
    E : float
        Electron kinetic energy in keV.
    seed : datetime.datetime
        UTC time associated with the selected field line.
    """

    plt.figure()
    ratios_str = [str(Fraction(x).limit_denominator()) for x in ratios]
    for i,j in zip(omegas,ratios_str):
        temp_df = df[df['omega'] == i]
        temp_df2 = temp_df[temp_df['MLAT'] <= 40]
        temp_df3 = temp_df2[temp_df2['MLAT'] >= -40]


        plt.plot(temp_df3['MLAT'], temp_df3['alpha_eq'],label=fr'${j} \Omega_{{ce,eq}}$')


    plt.xlabel(r'$\lambda$ (degrees)')
    plt.ylabel(r'$\alpha_{eq}$ (degrees)')
    plt.title(fr'$\alpha_{{eq}}$ Vs. $\lambda$ at {seed} (E = {E} keV)')
    plt.legend()
    plt.show()