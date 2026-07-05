# built-in modules
import logging

# 3rd party modules
import pandas as pd
import numpy as np

# project modules
from .._core.commons import __input_to_series__

def ETpartition_DWCS(data=None, ET='ET', T='T', E='E', H2O='wh2o', DownwardH2O='DownwardH2O',
                  CS_H2O_H2Opos_CO2neg='wh2o+wco2-', CS_H2O_H2Oneg_CO2neg='wh2o-wco2-', 
                  CS_H2O_H2Oneg_CO2pos='wh2o-wco2+', CS_H2O_H2Opos_CO2pos='wh2o+wco2+'): #NIGHT=None):
    logger = logging.getLogger('wvlt.partition.ETpartition_DWCS')
    logger.debug('Running ETpartition_DWCS, partitioning ET.')
    
    if isinstance(data, str): data = pd.read_file(data)
    else: data = data.copy()
    
    H2O = __input_to_series__(data, H2O)
    CS_H2O_H2Opos_CO2neg = __input_to_series__(data, CS_H2O_H2Opos_CO2neg)
    CS_H2O_H2Oneg_CO2neg = __input_to_series__(data, CS_H2O_H2Oneg_CO2neg)
    CS_H2O_H2Oneg_CO2pos = __input_to_series__(data, CS_H2O_H2Oneg_CO2pos)
    CS_H2O_H2Opos_CO2pos = __input_to_series__(data, CS_H2O_H2Opos_CO2pos)
    
    if data is None: data = pd.DataFrame()
    #if NIGHT is not None:
    #    islight = np.where((np.isnan(data[NIGHT]) == False) * (data[NIGHT]), 0, 1)
    #else:
    #    islight = np.array([1] * len(data))
    
    data[ET] = H2O
    data[T] = CS_H2O_H2Opos_CO2neg
    #logger.debug(f'T is {data[T]}')
    data[E] = CS_H2O_H2Opos_CO2pos
    data[DownwardH2O] = CS_H2O_H2Oneg_CO2neg + CS_H2O_H2Oneg_CO2pos
    logger.debug('Finished ETpartition_DWCS, partitioned ET.')
    return data


def _method_statistics_(data, average_period, cols_t_stat, dt, 
                        cols_corr, t_scale_thres):
    """
    function: Calculates the statistics of the conditional sampling process for provdided data and columns.
    call: _method_statistics_()
    Input:
        * data (pandas.DataFrame): Data to be processed.
        * average_period (str): Averaging period for which to calculate the statistics. Format: pandas time string, e.g. "30min". Possible specifications are s, min, h, d.
        * cols_t_stat (list): The columns in the data from which the statistics are being calculated. It makes sense, that its only the conditionally sampled columns.
        * dt (float): Sampling interval of the data (1/sampling frequency). Necessary to recalculate from the number of contigous data for an event to the event time scale. 
        * cols_corr (list): The columns in the data from which the correlations are being calculated.
    Return:
        pandas.DataFrame in long format with TIMESTAMP, variable and value. The variables are the method statistics. _t_fract denotes the time fraction of sampled events and _t_scale the average (mean) time scale of sampled events.
    """
    
    logger = logging.getLogger('wvlt.partition._method_statistics_')
    logger.debug('Calculating method statistics ')
    
    group_cols = ["TIMESTAMP", "natural_frequency"]
    
    # Reduce data to data of interest: the columns to do the statistics and the grouping columns  
    logger.debug(f'Calculating time fraction and scale for the columns {cols_t_stat}')
    cols_t_stat_t = cols_t_stat + group_cols
    data_time = data[cols_t_stat_t].copy()
    data_time["TIMESTAMP"] = data_time["TIMESTAMP"].dt.floor(average_period)
    
    # Time fraction of sampled events per Timestamp and Frequency
    time_frac = time_fraction(data_time, 
                              group_cols)
    
    # Mean time scale of sampled events per Timestamp and Frequency
    time_scal = time_scales(data_time, 
                group_cols,
                dt=dt,
                threshold=t_scale_thres)
    
    # Correlation coefficients
    logger.debug(f'Calculating correlation for the columns {cols_corr}')
    cols_corr_t = cols_corr + group_cols
    data_corr = data[cols_corr_t].copy()
    data_corr["TIMESTAMP"] = data_corr["TIMESTAMP"].dt.floor(average_period)
    corr_df = correlation_stats(data_corr,
                                group_cols)
    
    # Combine the calculated statistics
    data_stat = pd.concat([time_frac, time_scal, corr_df])
    logger.debug(f"Method statistics in data_stat: {data_stat}")
    
    return data_stat

def correlation_stats(data, group_cols):
    # Calculate correlations per partitioned flux/quadrant and per averaging time
    
    logger = logging.getLogger('wvlt.partition.correlation_stats')
    logger.debug(f'Calculating correlations, data is {data}')
    
    # 1. Correlate
    corr_m = (data.sort_values(group_cols)
              .groupby(group_cols)
              .corr()
              )
    
    # # OLD WAY USING MASK:
    # # 2. Masking to remove doubled values and diagonal from correlation matrix
    # # Create a mask that matches the SMALL correlation matrix
    # # for each combination of the grouping variables (e.g., 2x2)
    # # Get the number of variables (e.g., 2 if co2 and h2o)
    # n_vars = len(data.columns.drop(group_cols))
    # small_mask = np.triu(np.ones((n_vars, n_vars)), k=1).astype(bool)
    
    # # Tile the mask to match the full length of corr_m
    # # This repeats the 2x2 mask for every single group/frequency
    # mask = np.tile(small_mask, (len(corr_m) // n_vars, 1))

    # # 3. Filter results
    # corr_m = (corr_m
    #         .where(mask)
    #         .stack()
    #         .rename_axis(group_cols + ['corr_1', 'corr_2'])
    #         .reset_index(name='value')
    #         )
    
    # NEW WAY USING INDICES, better if missing values 
    # -- which are almost impossible, after gapfilling etc.
    # It produces the same output as the mask before.
    corr_m = corr_m.stack()
    corr_m.index.names = group_cols + ['corr_1', 'corr_2']
    corr_m = corr_m[corr_m.index.get_level_values('corr_1') < corr_m.index.get_level_values('corr_2')]
    corr_m = corr_m.reset_index(name='value')
    
    # 4. Combine 
    corr_m['variable'] = corr_m["corr_1"] + '-' + corr_m["corr_2"] + "_r"
    
    # 5. Delete old level indices and sort dataframe
    corr_m = (corr_m.drop(["corr_1", "corr_2"], axis="columns")
              .sort_values(by=["variable"] + group_cols))
    
    logger.debug(f"Calculated correlations, corr_m is {corr_m}")
    return corr_m
    
    

def time_fraction(data, group_cols):
    logger = logging.getLogger('wvlt.partition.time_fraction')
    logger.debug(f'Calculating time fraction, data is {data}')
    
    time_frac = (data.groupby(group_cols)
                 .agg(lambda x: x.dropna().astype(bool).mean())
                 .add_suffix("_t_fract")
                 .reset_index()
                 #.sort_values(by=["TIMESTAMP_av", "natural_frequency"])
                 .melt(id_vars=group_cols))
    
    logger.debug(f'Calculated time fraction, time_frac is {time_frac}')
    
    return time_frac


def time_scales(df, group_cols, dt, threshold=10):
    logger = logging.getLogger('wvlt.partition.time_scales')
    logger.debug(f"df columns: {df.columns}")
    logger.debug(f"df data: {df.head()}")
    
    # Ensure the dataframe is sorted so "consecutive" actually means time-consecutive
    df = df.sort_values(group_cols).reset_index(drop=True)
    
    # Create a 'Local Index' (0, 1, 2, 3...) for each unique group
    # This acts as a counter: "This is the 1st sample for 5Hz, this is the 2nd..."
    df['local_idx'] = df.groupby(group_cols).cumcount()
    
    results = {}
    
    # remove TIMESTAMP and natural_frequency from target columns
    target_cols = [c for c in df.columns if c not in group_cols + ['local_idx']]
    for col in target_cols:
        # Identify where non-zero values are that are not NaN
        is_nz = (df[col] != 0) & df[col].notna()
        
        # Get the index of every non-zero row
        nz_indices = df.index[is_nz]
        
        if len(nz_indices) == 0:            
            # If all were NaN, set to NaN, if there were data but just no events take 0.0
            valid_counts = df.groupby(group_cols)[col].count()
            results[col + "_t_scale"] = pd.Series(np.where(valid_counts == 0, np.nan, 0.0), index=valid_counts.index)
            continue
        
        # Identify the "starts" of events (where it goes from 0 to non-zero)
        # We also treat it as a new event if the gap since the last non-zero is > threshold
        
        # We group by and then calculate the index difference
        def get_event_ids(group):
            # Difference between row indices within THIS group
            
            # Calculate the gap between current non-zero and the previous non-zero
            # gaps = group.index.to_series().diff()
            gaps = group['local_idx'].diff()
            
            # If gap > (threshold + 1), it's a brand new "buffered" event
            new_event = gaps >= (threshold + 1)
            
            # Create the Group IDs for these events
            event_ids = new_event.cumsum()
            return event_ids
        
        # Create a temporary dataframe of just the non-zero hits to aggregate
        temp_nz = df.loc[nz_indices, group_cols + ['local_idx']].copy()
        temp_nz['event_id'] = temp_nz.groupby(group_cols, group_keys=False).apply(get_event_ids, include_groups=False)
        
        # Count the size of each event, then mean per group
        streak_counts = temp_nz.groupby(group_cols + ['event_id']).size()
        calculated_means = (streak_counts
                            .groupby(level=list(range(len(group_cols))))
                            .mean()) # level = list(range(len(group_cols))) to get rid of ['event_id'] as grouping because we want mean per grouping (TIMESTAMP + Frequency), not per event
        # Build smart defaults: NaN if the group was all-NaN, 0.0 if it had actual numbers (zeros)
        valid_counts = df.groupby(group_cols)[col].count()
        group_defaults = pd.Series(np.where(valid_counts == 0, np.nan, 0.0), index=valid_counts.index)
        # Combine them: calculated means take priority; missing groups fall back to defaults
        results[col + "_t_scale"] = calculated_means.combine_first(group_defaults)
    
    # logger.debug(f"results: {results}")
    results = pd.DataFrame(results).reset_index().melt(group_cols)
    logger.debug(f"Number of consecutive samples as event scale before calculation to time scale: {results}")
    
    # Convert to mean amount of consecetively sampled events 
    # to time scale using sampling interval dt = 1/fs (sampling frequency)
    if dt is not None:
        results["value"] = results["value"] * dt
    
    logger.debug(f"Calculated time scale {results.head()}.")
    logger.debug(f"Column names of time_scal {results.columns}")
    
    return results



def int_tests(data, f0, n, f_low, variables_not, calc_na=False):
    """
    function: Calculates the Stationarity (STA), and Ogive Test (OG) and assigns stepped QA/QC flags. It checks for discrete wavelet frequency band overlaps and logs safety warnings accordingly.
    call: int_tests()
    Input:
        * data (pandas.DataFrame): Long-format wavelet-decomposed data containing 'natural_frequency', 'variable', 'TIMESTAMP', and 'value' columns.
        * f0 (float): Baseline integration frequency threshold (Hz).
        * n (int or float): For the stationarity test, the ratio of the period duration of the lower integration period (T* < T) to the normal period duration (T = 1/f0). Typically 6, following the traditional Stationarity test (5min/30min).
        * f_low (float): The lower frequency for the ogive test
        * variables_not (str or tuple of str): Suffix or tuple of suffixes passed to .str.endswith() to exclude specific variables from calculation.
        * calc_na (bool, default False): if False, if any of the frequencies has NaN values, the integrated flux is set NA instead of integrating over only the remaining frequencies (0 if all frequencies have NaN values).
    Return:
        The processed Pandas DataFrame containing integrated band sums, STA, OG and QAQC step flags (0 <= 30%, 1 <= 100%, 2 > 100%).
    """
    
    logger = logging.getLogger('wvlt.partition.int_tests')
    logger.debug("Calculating Stationarity and Integration Scale Test")
    
    # PRE-CALCULATION for Stationarity test
    f_high = 1/((1/f0)/n) # higher frequency  
    
    # FILTER DATA TO RELEVANT VARIABLES
    dataf = data[~data['variable'].str.endswith(variables_not)]
    
    # WARNINGS
    # output a warning for f_high and f_low which dicrete wavelet level is taken for integration
    freq = data['natural_frequency'].unique()
    
    if f_low <= freq.min():
        logger.warning(
            f"Specified lower frequency {f_low} Hz is below or equal to the lowest available "
            f"wavelet frequency ({freq.min()} Hz). OG will integrate the entire spectrum. "
            f"This can be fine, but please know what you're doing since it includes the approximation coefficient."
        )
    if f_low > freq.min():
        # Filter arrays
        higher_vals = freq[freq >= f_low]
        lower_vals = freq[freq < f_low]
        
        lower_valsf0 = freq[freq < f0]
        
        # Safety check for edge cases
        if higher_vals.size > 0 and lower_vals.size > 0:
            max_band_high_freq = higher_vals.min() # highest frequency band: max frequency border
            max_band_low_freq = lower_vals.max() # highest frequency band: min frequency border
            
            max_band_low_freqf0 = lower_valsf0.max()
            
            if max_band_low_freq == max_band_low_freqf0:
                logger.warning(f"Specified lower integration frequency for OG was {f_low} Hz ({1/f_low} s as period duration). From available frequencies of the wavelet transform integrate from highest frequency up to INCLUDING the frequency band from {max_band_high_freq} Hz ({1/max_band_high_freq} s) to {max_band_low_freq} Hz ({1/max_band_low_freq} s). This corresponds to the same frequency band as the normal integration frequency f0 {f0}. The test will by definition equal to 0. Please decrease lower integration frequency f_low for OG to be corrsponding to a lower frequency band than f0. Lower than {max_band_low_freqf0}.")
            
            if (1/f_low) <= (((1/max_band_high_freq)+(1/max_band_low_freq))/2):
                logger.warning(f"Specified lower integration frequency for OG was {f_low} Hz ({1/f_low} s as period duration). From available frequencies of the wavelet transform integrate from highest frequency up to INCLUDING the frequency band from {max_band_high_freq} Hz ({1/max_band_high_freq} s) to {max_band_low_freq} Hz ({1/max_band_low_freq} s). Please make sure you know what you're doing - potentially changing the integration frequency slightly below {max_band_high_freq} Hz ({1/max_band_high_freq} s) or {max_band_low_freq} Hz ({1/max_band_low_freq} s).")
            else:
                logger.info(f"Specified lower integration frequency for OG was {f_low} Hz ({1/f_low} s as period duration). From available frequencies of the wavelet transform integrate from highest frequency up to INCLUDING the frequency band from {max_band_high_freq} Hz ({1/max_band_high_freq} s) to {max_band_low_freq} Hz ({1/max_band_low_freq} s).")
        else:
            logger.warning(f"Specified lower integration frequency for OG was {f_low} Hz ({1/f_low} s as period duration). This is already the highest frequency of either the lowest or the highest available frequency band of the wavelet transform. Please make sure you know what you're doing.")
    

    # STATIONARITY TEST
    # normal integration
    data0 = (dataf[(np.isnan(dataf['natural_frequency']) == False) * (dataf['natural_frequency'] >= f0)]
        .groupby(['variable', 'TIMESTAMP'])['value']
        .agg(lambda x: x.sum(skipna=calc_na))
        .reset_index()
        .rename(columns={'value': 'value_f0'})
        .set_index(['variable', 'TIMESTAMP'])
        )
    
    # integration only up to higher frequency f_high
    data1 = (dataf[(np.isnan(dataf['natural_frequency']) == False) * (dataf['natural_frequency'] >= f_high)]
        .groupby(['variable', 'TIMESTAMP'])['value']
        .agg(lambda x: x.sum(skipna=calc_na))
        .reset_index()
        .rename(columns={'value': 'value_f_high'})
        .set_index(['variable', 'TIMESTAMP'])
        )
    
    # INTEGRATION SCALE TEST (adjusted ogive test with absolute values)
    # normal integration with ABSOLUTE values
    # data0_abs = (dataf[(np.isnan(dataf['natural_frequency']) == False) * (dataf['natural_frequency'] >= f0)]
    #     .groupby(['variable', 'TIMESTAMP'])['value']
    #     .agg(lambda x: abs(x).sum(skipna=calc_na))
    #     .reset_index()
    #     .rename(columns={'value': 'value_f0_abs'})
    #     .set_index(['variable', 'TIMESTAMP'])
    #     )
    
    # integrate up to smaller frequency f_low, but with ABSOLUTE values
    # data2 = (dataf[(np.isnan(dataf['natural_frequency']) == False) * (dataf['natural_frequency'] >= f_low)]
    #     .groupby(['variable', 'TIMESTAMP'])['value']
    #     .agg(lambda x: abs(x).sum(skipna=calc_na))
    #     .reset_index()
    #     .rename(columns={'value': 'value_f_lowabs'})
    #     .set_index(['variable', 'TIMESTAMP'])
    #     )
    
    # OGIVE TEST
    # Ogive test, adjusted from Charuchittipan 2014, Foken 2006
    # Filter for the low frequency band to integrate: f_low <= frequency < f0
    data_ot_filtered = dataf[(np.isnan(dataf['natural_frequency']) == False) & 
                             (dataf['natural_frequency'] >= f_low) & 
                             (dataf['natural_frequency'] < f0)]
    data_ot_sorted = data_ot_filtered.sort_values(by='natural_frequency', ascending=False)
    
    # Helper function to track cumulative max safely matching calc_na setup
    def calc_ot_max(x):
        if not calc_na and x.isna().any():
            # If not calc_na and there is a NaN in there, all just NaN
            return np.nan
        # If calc_na, calculate with remaining values, if not calc_na and were here, just take all values
        cleaned = x.dropna() if calc_na else x
        if cleaned.empty:
            # If no values in the cleaned, just take 0 if calc_na, otherwise NaN
            return 0.0 if calc_na else np.nan
        # Calculate ogive using cumsum, take the absolute of it for max calculation
        return cleaned.cumsum().abs().max()
    
    data_ot_res = (data_ot_sorted.groupby(['variable', 'TIMESTAMP'])['value']
        .agg(calc_ot_max)
        .reset_index()
        .rename(columns={'value': 'max_abs_cumsum'})
        .set_index(['variable', 'TIMESTAMP'])
        )
    
    # COMBINE DATA
    # datanew = data0.join([data1, data2, data0_abs], how='inner').reset_index()
    datanew = data0.join([data1], how='inner')
    datanew = datanew.join(data_ot_res, how='left').reset_index()
    
    # CALCULATE STATISTICS
    datanew['STA'] = abs((datanew['value_f_high'] - datanew['value_f0']) / datanew['value_f0']) * 100
    # datanew['IST'] = (datanew['value_f_lowabs'] - datanew['value_f0_abs']) / datanew['value_f0_abs'] * 100
    # Force any tiny negative floating-point noise (e.g. -0.0000000000000119009068481075) to be exactly 0
    # datanew['IST'] = datanew['IST'].clip(lower=0)
    datanew['OG'] = abs(datanew['max_abs_cumsum'] / datanew['value_f0']) * 100
    
    # QUALITY FLAGS
    sta_conditions = [
        datanew['STA'] <= 30,   # Condition 1 -> Flag 0
        datanew['STA'] <= 100,  # Condition 2 -> Flag 1
        datanew['STA'].notna()  # Condition 3 (Anything > 100 and not NaN) -> Flag 2
    ]
    sta_choices = [0, 1, 2]
    datanew['QAQC_STA'] = np.select(sta_conditions, sta_choices, default=np.nan)
    
    # ist_conditions = [
    #     datanew['IST'] <= 30,      # Condition 1 -> Flag 0
    #     datanew['IST'] <= 100,     # Condition 2 -> Flag 1
    #     datanew['IST'].notna()     # Condition 3 (Anything > 100 and not NaN) -> Flag 2
    # ]
    # ist_choices = [0, 1, 2]
    # datanew['QAQC_IST'] = np.select(ist_conditions, ist_choices, default=np.nan)
    
    ot_conditions = [
        datanew['OG'] <= 30,       # Condition 1 -> Flag 0
        datanew['OG'] <= 100,      # Condition 2 -> Flag 1
        datanew['OG'].notna()      # Condition 3 (Anything > 100 and not NaN) -> Flag 2
    ]
    ot_choices = [0, 1, 2]
    datanew['QAQC_OG'] = np.select(ot_conditions, ot_choices, default=np.nan)
    
    logger.debug(f"datanew: {datanew}")
    return datanew 
    
    
    
    
    
    
    
    
    
    
    
    
