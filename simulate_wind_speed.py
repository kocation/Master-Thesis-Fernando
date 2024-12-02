import numpy as np


def simulate_wind_speed(I=0.14, l=340.2, V_10=8, f_cutout=0.5, T_max=600, dt=0.05):
    """
    Simulates wind speed time series using a Kaimal spectrum method.
    
    Parameters:
    I (float): Turbulence intensity.
    l (float): Length scale (related to turbulence).
    V_10 (float): Mean wind speed at 10 meters height.
    f_cutout (float): Cutoff frequency (Hz).
    T_max (float): Total time for simulation (in seconds).
    dt (float): Time step for simulation.
    
    Returns:
    np.array: Simulated wind speed time series.
    np.array: Corresponding time array.
    """
    
    df = 1 / T_max  # Frequency resolution
    f_p = np.arange(df, f_cutout, df)  # Frequency range
    time = np.arange(0, T_max, dt)  # Time array

    # Define random phase (epsilon)
    np.random.seed(41) # used to reproduce the same results and create comparable simulations
    epsilon = 2 * np.pi * np.random.rand(len(f_p))

    # Wind spectrum (S_wind)
    S_wind = 4 * I**2 * V_10 * l / (1 + 6 * f_p * l / V_10)**(5/3)

    # Calculate b_p and omega_p
    b_p = np.sqrt(2 * S_wind * df)
    omega_p = 2 * np.pi * f_p

    # Initialize wind speed time series (V_hub)
    V_hub = np.zeros(len(time))
    
    # Compute wind speed at each time step
    for i in range(len(time)):
        V_hub[i] = V_10 + np.sum(b_p * np.cos(omega_p * time[i] + epsilon))

    return V_hub