# built-in modules
import logging

# 3rd party modules
import pandas as pd

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


def _method_statistics_(data, average_period, cols_for_stat):
    logger = logging.getLogger('wvlt.partition._method_statistics_')
    logger.debug('Calculating method statistics ')
    
    cols_for_stat.extend(["TIMESTAMP", "natural_frequency"])
    data = data[cols_for_stat].copy()
    #metstat_data_l = []
    data["TIMESTAMP"] = data["TIMESTAMP"].dt.floor(average_period)
    
    # Time fraction of sampled events per Timestamp and Frequency
    logger.debug(f'Calculating time fraction, data is {data}')
    time_frac = (data.groupby(["TIMESTAMP", "natural_frequency"])
                 .agg(lambda x: x.astype(bool).mean())
                 .add_suffix("_t_fract")
                 .reset_index()
                 #.sort_values(by=["TIMESTAMP_av", "natural_frequency"])
                 .melt(id_vars=["TIMESTAMP", "natural_frequency"]))
    logger.debug(f'Calculated time fraction, time_frac is {time_frac}')
    
    return time_frac



