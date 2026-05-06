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
                 .agg(lambda x: x.astype(bool).mean())
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
        # Identify where non-zero values are
        is_nz = df[col] != 0
        
        # Get the index of every non-zero row
        nz_indices = df.index[is_nz]
        
        if len(nz_indices) == 0:
            results[col + "_t_scale"] = df.groupby(group_cols).size().astype(float) * 0
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
        results[col + "_t_scale"] = (streak_counts
                                     .groupby(level=list(range(len(group_cols))))
                                     .mean()) # level = list(range(len(group_cols))) to get rid of ['event_id'] as grouping because we want mean per grouping (TIMESTAMP + Frequency), not per event
    
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



