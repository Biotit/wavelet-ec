
# built-in modules

# 3rd party modules
import numpy as np
import pandas as pd
import logging
from scipy.optimize import curve_fit

# project modules
from .commons import j2sj

logger = logging.getLogger('corrections')

def density_correction(data, correction_density, pTq_cols,
                         average_period):
    """
    function: Applies a density correction to with open-path analyzers measured scalars according to Detto & Katul, 2007 (https://doi.org/10.1007/s10546-006-9105-1. In contrast to all other parts of the code, in here the units are important!)
        The code is modified based upon the code by Zahn, 2022 (PartitioningMethods) and Skaggs, 2018 (Fluxpart)
    call: density_correction()
    Input:
        * data (pandas.DataFrame): Data to be corrected.
        * vars_correct (list): List of variables (coloum names of data) to be density corrected, e.g. ["h2o", "co2"]. Unit has to be mmol/m3.
        * pTq_cols (list): List of column names from data with pressure, sonic temperature, and water vapor, e.g. ["Pressure", "Ts", "h2o"]. Temperature in °C, pressure in kPa, water vapor in mmol/m3.
        * average_period (str, default '30min'). The averaging period for which to calculate the mean densities and temperature for the formulas.
    Return:
        A new class object named var_ with class attributes data and saved. Data includes the averaged wavelet transformed, cross calculated variables. saved_files contains strings with paths to where the saved files are placed. If save return as test = main(), access data via test.data or test.saved.
    """
    
    # Starting logger --------------------------------------------------------
    logger = logging.getLogger('waveletec.corrections.density_correction')
    logger.info("In contrast to all other parts of the code, in here the units are important!")
    logger.info("Please make sure your input data for h2o and the scalars to correct is in mmol/m3.")
    logger.info("Please make sure your Temperature is in °C, and your pressure in kPa")
    #logger.debug(f"correction_density: {correction_density}")
    #logger.debug(f"pTq_cols: {pTq_cols}")
    #logger.debug(f"data.columns: {data.columns}")
    
    # Constants ---------------------------------------------------------------
    Ra = 287 # J/kg/K gas constant for dry air
    #R = 8.314462618 # J/mol/K universal gas constant
    Mq = 0.018016 # kg/mol Molar weight of water (vapor)
    Ma = 0.0289645 # kg/mol Molar weight of dry air
    mu = Ma / Mq # ratio of dry air mass to water vapor mass
    
    # Assertions --------------------------------------------------------------
    # Input columns with Pressure in kPa, sonic Temperature in °C,  water vapor in mmol/m3
    p_col, T_col, q_col = pTq_cols
    assert all(col in data.columns for col in pTq_cols), f"Column for pressure, temperature or water vapor missing in input data. Column names given as pTq_cols: {pTq_cols}, column names of the data: {data.columns}."
    
    # Main data loop over averaging periods -----------------------------------
    # Need to compute fluctuations and mean per average_period
    data["TIMESTAMP_av"] = data["TIMESTAMP"].dt.floor(average_period)
    corrected_chunks = []
    for av_period, data_a in data.groupby("TIMESTAMP_av"):
        # logger.debug(f"Correcting for {av_period}.")
        data_a = data_a.copy()
        if data_a.empty or data_a[correction_density].isnull().all().all():
            # In case only NA in averaged time
            corrected_chunks.append(data_a)
            continue
        # Create temporary dataframe to collect temporary variables
        data_temp = pd.DataFrame(index=data_a.index)
    
        # Find mean and perturbations / fluctuations --------------------------
        # vars_mean: mean of the vars_correct in mol/m3,
        # all vars_correct need to be in mmol/m3
        vars_mean = data_a[correction_density].mean() * 10**-3
        # Mean of water vapor
        q_mean = data_a[q_col].mean() * 10**-3
        data_temp['q_pert'] = data_a[q_col] * 10**-3 - q_mean
        
        # Calculate air density in kg/m3 --------------------------------------
        data_temp['rho_moist_air_mass'] = 1000 * data_a[p_col] / (Ra * (273.15 + data_a[T_col]))
        # Calculate dry air density in mol/m3 for further calculation
        data_temp['rho_dry_air_mass'] = ( 
            data_temp["rho_moist_air_mass"] - (data_a[q_col] * 10**-3 * Mq)
            )
        # Calculate air density as mol/m3 for formula further down
        data_temp['rho_dry_air'] = data_temp['rho_dry_air_mass'] / Ma
        rho_dry_air_mean = data_temp["rho_dry_air"].mean()
        
        # Thermodynamic temperture and its perturbations ----------------------
        # Calculate specific humidity kg/kg
        q = data_a[q_col] * 10**-3 * Mq / (data_temp["rho_moist_air_mass"])
        # Calculate thermodynamic (= real) temperature from sonic temperature [C]
        data_temp["T"] = (data_a[T_col] + 273.15) / (1.0 + 0.51 * q) - 273.15  
        # Get mean temperature
        T_mean = data_temp["T"].mean() 
        data_temp["T_pert"] = data_temp["T"] - T_mean
        # Calculate "mixing ratio" using moles [mol_q/mol_air], see Paper for Reference
        sigmaq = q_mean / rho_dry_air_mean
            
            
        # Apply corrections Eq. 5 in paper -----------------------------------
        # Loop through the variables to correct
        for var in correction_density:
            # logger.debug(f"Correcting {var}")
            # By adding left and right of the equal the mean,
            # we can get the total time series corrected not only the perturbation
            # (as in the formula)
            # See Fluxpart code, they did the same.
            data_a[var] = data_a[var] + ( 
                mu * (vars_mean[var]/rho_dry_air_mean) * data_temp['q_pert']
                + vars_mean[var] * (1.0 + mu * sigmaq) * (data_temp["T_pert"] / (T_mean + 273.15))
                ) * 10**3 # mmol/m3
            
        del data_temp
        corrected_chunks.append(data_a)
    
    return pd.concat(corrected_chunks).sort_values("TIMESTAMP")


def mauder2013(x, q=7):
    # it does not do the check for n consecutive spikes 
    x = np.array(x)
    x_med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - x_med))
    bounds = (x_med - (q * mad) / 0.6745, x_med + (q * mad) / 0.6745)
    #print("median", x_med, "mad", mad, "bounds", bounds)
    x[x < min(bounds)] = np.nan
    x[x > max(bounds)] = np.nan

    #if fill is not None:
    #    x = fill(pd.Series(x) if fill in (pd.Series.ffill, pd.Series.interpolate) else x)
    return x


def __despike__(X, method=mauder2013):
    N = len(X)
    X = method(X)
    Xna = np.isnan(X)
    try:
        X = np.interp(np.linspace(0, 1, N), 
                        np.linspace(0, 1, N)[Xna == False],
                X[Xna==False])
    except Exception as e:
        logger.error(f"UserWarning: {str(e)}")
    return X


def __fit_whitenoise__(spec, freqs=None, fmax=5, a=1):
    """Fit a white noise spectrum to the data
    
    Parameters
    ----------
    spec : array
        The spectrum to fit
    fmax : int
        The maximum frequency to fit
    a : float
        The exponent of the spectrum ("color" of the noise)
    Returns
    -------
    fit : array
        The fitted spectrum
    """
    if freqs is None:
        dt = 20
        freqs = [1/j2sj(j, 1/dt) for j in np.arange(1, len(spec)+1)]
    curve_0 = lambda f, b: np.log((f**a)*b)
    specna = np.where(np.isnan(spec[:fmax]) | (
        spec[:fmax] <= 0), False, True)
    try:
        freqs = np.array(freqs)
        spec = np.array(spec)
        params_0, _ = curve_fit(curve_0, freqs[:fmax][specna], np.log(
            spec[:fmax][specna]), bounds=(0, np.inf))
    except Exception as e:
        logger.error(f"UserWarning: {str(e)}")
        params_0 = [0]
    fit = np.array([(f**a)*params_0[0] for f in freqs])
    return type('var_', (object,), {'curve_0': curve_0, 'params_0': params_0, 'fit': fit})


def smooth_2d_data(data, **kwargs):
    if data.shape[0] == 1:
        return smooth_data(data, **kwargs)
    else:
        for i in range(data.shape[0]):
            data[i] = smooth_data(data[i], **kwargs)
        return data


def smooth_data(data, method='convolve', smoothing=3, **kwargs):
    if smoothing == 0:
        return data
    if method == 'convolve':
        return np.convolve(data, np.ones((smoothing,)) /
                           smoothing, mode='same')
    # if method == 'gaussian':
    #     return gaussian_filter1d(data, sigma=smoothing, **kwargs)
    # if method == 'savgol':
    #     return savgol_filter(data, window_length=smoothing, polyorder=2, **kwargs)
    if method == 'moving':
        stat = kwargs.get('stat', 'mean')
        if stat == 'mean':
            return data.rolling(smoothing, min_periods=1).mean()
        if stat == 'std':
            return data.rolling(smoothing, min_periods=1).std()
        if stat == 'max':
            return data.rolling(smoothing, min_periods=1).max()
        if stat == 'quantile':
            return data.rolling(smoothing, min_periods=1).quantile(0.5)
    if method == 'repeat':
        # Create array of indices to group by
        chunk_indices = np.floor(np.arange(len(data)) / smoothing).astype(int)
        # Group by those indices and apply the function
        return pd.Series(data).groupby(
            chunk_indices).transform(np.nanmean)

    return data
    for i in range(len(φcs)):
        φcs[i] = np.convolve(φcs[i], np.ones(
            (smoothing,))/smoothing, mode='same')
