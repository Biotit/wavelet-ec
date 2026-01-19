"""
This script is a key part of the following publications:
    - Herig Coimbra, Pedro Henrique and Loubet, Benjamin and Laurent, Olivier and Mauder, Matthias and Heinesch, Bernard and 
    Bitton, Jonathan and Delpierre, Nicolas and Depuydt, Jérémie and Buysse, Pauline, Improvement of Co2 Flux Quality Through 
    Wavelet-Based Eddy Covariance: A New Method for Partitioning Respiration and Photosynthesis. 
    Available at SSRN: https://ssrn.com/abstract=4642939 or http://dx.doi.org/10.2139/ssrn.4642939

The main function is:  
- process
- main

"""

# built-in modules
import os
import re
import warnings
import logging
import copy
import time
import datetime
import glob

# 3rd party modules
from functools import reduce
import numpy as np
import pandas as pd
import itertools
from scipy.optimize import curve_fit
import pywt
import yaml
import tqdm # progress bar
try: 
    import pycwt
except ImportError as e:
    pycwt = None
    pass
try: 
    import fcwt
except ImportError as e:
    fcwt = None
    pass

# project modules
from . import commons as hc24
from .read_data import loaddatawithbuffer
from .wavelet_functions import universal_wt, formula_to_vars, prepare_signal, bufferforfrequency_dwt, bufferforfrequency
from ..partitioning import coimbra_et_al_2025 as ptt
from ..partitioning import schoendorf as pttET
from .._extra import eddypro_tools as eddypro

logger = logging.getLogger('wvlt.pipeline')


# def sample_raw_data(input_path, datetimerange, acquisition_frequency=20, fileduration=30, **kwargs):
#     raw_kwargs = {'path': input_path, 'fkwargs': {'dt': 1/acquisition_frequency}}
#     kwargs['fmt'] = kwargs.get('fmt', {})
#     if 'gas4_name' in kwargs.keys(): kwargs['fmt'].update({kwargs.pop('gas4_name'): '4th gas'})
#     raw_kwargs.update({k: v for k, v in kwargs.items() if k in ['fmt']})

#     ymd = [datetimerange.split('-')[0], datetimerange.split('-')[1], f'{fileduration}min']
#     _, _, _f = ymd
#     ymd = hc24.list_time_in_period(*ymd, '1D', include='both')

#     for ymd_i, yl in enumerate(ymd):
#         data = hc24.loaddatawithbuffer(
#             yl, d1=None, freq=_f, buffer=0, f_freq=_f, **raw_kwargs)
#         break
#     return data


# def sample_raw_data(input_path, datetimerange, acquisition_frequency=20, fileduration=30, processduration='1D'):
#     ymd = [datetimerange.split(
#         '-')[0], datetimerange.split('-')[1], f'{fileduration}min']
#     _, _, _f = ymd
#     ymd = waveletec._core.list_time_in_period(
#         *ymd, processduration, include='both')

#     for ymd_i, yl in enumerate(ymd):
#         data = waveletec._core.loaddatawithbuffer(
#             yl, d1=None, freq=_f, buffer=0, f_freq=_f, **{'path': input_path, 'fkwargs': {'dt': 1/acquisition_frequency}})
#         break
#     return data

# raw_kwargs = {'path': input_path, 'fkwargs': {'dt': 1/acquisition_frequency}}
# kwargs['fmt'] = kwargs.get('fmt', {})
# if 'gas4_name' in kwargs.keys(): kwargs['fmt'].update({kwargs.pop('gas4_name'): '4th gas'})
# raw_kwargs.update({k: v for k, v in kwargs.items() if k in ['fmt']})

# ymd, raw_kwargs, output_folderpath = None, verbosity = 1,
# overwrite = False, processing_time_duration = "1D",
# internal_averaging = None, dt = 0.05, wt_kwargs = {},
# method = "dwt", averaging = 30, **kwargs)

    

# def eddypro_wavelet_run(site_name, input_path, outputpath, datetimerange, acquisition_frequency=20, fileduration=30, 
#          processduration='1D', integration_period=None, preaverage=None,
#          covariance = None, variables_available=['u', 'v', 'w', 'ts', 'co2', 'h2o'], denoise=0, deadband=[], 
#          method = 'dwt', wave_mother='db6', **kwargs):
#     local_args = locals()

#     if outputpath is not None:
#         hc24.start_logging(outputpath)

#         # Select output file path
#         if method == 'cov':
#             outputpath = str(os.path.join(outputpath, str(site_name)+'{}_{}.csv'))
#         else:
#             outputpath = str(os.path.join(outputpath, 'wavelet_full_cospectra', str(site_name)+'_CDWT{}_{}.csv'))

#         # Save args for run
#         hc24.mkdirs(outputpath)
#         with open(os.path.join(os.path.dirname(os.path.dirname(outputpath)), f'log/setup_{datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")}.yml'), 'w+') as stp:
#             yaml.safe_dump(local_args, stp)

#     # Select covariances
#     # x*y → Cov(x, y)
#     # x*y|x*z|x*... → Cov(x, y)|Cov(x, z),Cov(x, ...)
#     if covariance is None:
#         covariance = hc24.available_combinations(
#             hc24.DEFAULT_COVARIANCE, variables_available)

#     # RUN WAVELET FLUX PROCESSING
#     # ymd = [START_DATE, END_DATE, FILE_FREQUENCY]
#     raw_kwargs = {'path': input_path, 'fkwargs': {'dt': 1/acquisition_frequency}}
#     kwargs['fmt'] = kwargs.get('fmt', {})
#     if 'gas4_name' in kwargs.keys(): kwargs['fmt'].update({kwargs.pop('gas4_name'): '4th gas'})
#     raw_kwargs.update({k: v for k, v in kwargs.items() if k in ['fmt']})
#     data = wavelet_functions.load_data_and_loop(ymd = [datetimerange.split('-')[0], datetimerange.split('-')[1], f'{fileduration}min'],
#                                          output_path = outputpath,
#                                          varstorun = covariance,
#                                          averaging = [fileduration],
#                                          processing_time_duration = processduration,
#                                          method = method,
#                                          wt_kwargs = {'fs': acquisition_frequency, 'wavelet': wave_mother},
#                                          raw_kwargs = raw_kwargs,
#                                          verbosity=5)
#     return data


def condition_sampling_partition(site_name, output_folderpath, 
                                 integration_period=None, 
                                 variables_available=['w', 'co2', 'h2o'], 
                                 newlog=False, **kwargs):
    """
    function: Read avaraged and integrated cospectra file, apply conditional sampling, i.e. different partitioning algorithms, and save result.
    call: condition_sampling_partition()
    Input: 
        * site_name (str): Site name of the data to be loaded in. Nessessary to construct file names to be loaded. See variable output_folderpath for more information.
        * output_folderpath (str): Path to folder where the input and output files files are saved. Inside this folder there has to be a file with the pattern os.path.join(output_folderpath, f"{site_name}_CDWT_fulldata_integrated_*min.csv"). Usually produced by integrate_full_spectra_into_file() or by process(). 
        * integration_period (int, default None): If different files with different integration_period inside the output_folderpath, this helps to find the correct file for conditional sampling. In those functions it is the integration period of the wavelength signal in s. Works as a high-pass filter for the wavelet cospectra (as f0 = 1/integration_period) inside integrate_cospectra(). Also relevant for the filename of saved data. It gets constructed similar to os.path.join(output_folderpath, str(site_name)+f'_CDWT_partitioning_H2O.csv' dependent on the used partitioning algorithm.
        * variables_available (list, default ['w', 'co2', 'h2o']): From which variables are data available. Necessary to know, which partitioning algorithms can be run.
        * newlog (bool, default False): if new log file in the subfolder log inside the output_folderpath is created using start_logging(). Useful if the function condition_sampling_partition() is called on its own, e.g. outside of eddypro_wavelet_run() or with time delay after other functions.
        **kwargs
    Return:
        No return.
    """
    # RUN PARTITIONING
    #dst_path = os.path.join(output_folderpath, str(
    #    site_name)+f'_CDWT_full_cospectra.csv')
    # variables_available=['u', 'v', 'w', 'ts', 'co2', 'h2o']
    
    # activate new logging file? Useful if function is called on its own, e.g. outside of eddypro_wavelet_run and with time delay after the process().
    if (output_folderpath is not None) and newlog:
        hc24.start_logging(output_folderpath)
        
    logger = logging.getLogger('wvlt.handler.condition_sampling_partition')
    
    # to be able to have different integration_period = 1/f0, hence different high pass filters in the folder
    # search for the pattern with variable minutes
    if not integration_period:
        pattern = os.path.join(output_folderpath, f"{site_name}_CDWT_fulldata_integrated_*min.csv")
        matches = glob.glob(pattern)
        if not matches:
            raise FileNotFoundError(f"No file found matching {pattern}")
        if len(matches) > 1:
            logger.warning(f"Multiple files of integrated cospectra found: {matches}, taking the first one for partitioning! Set integration_period if want to use another one.")
        dst_path = matches[0]
        logger.info(f'Taking the file {dst_path} for partitioning.')
    else:
        dst_path = os.path.join(output_folderpath + str(site_name) + f"_CDWT_fulldata_integrated_{integration_period//60}min" + ".csv")
        logger.debug(f"Specified integration_period={integration_period}. Hence taking the file {dst_path}")
        

    h2o_dw_required_variables = ['w','co2','h2o']
    is_lacking_variable = sum([v not in variables_available for v in h2o_dw_required_variables])
    if is_lacking_variable:
        logger.debug(f'For partition_DWCS_H2O with {h2o_dw_required_variables} lacking variables.')
    if not is_lacking_variable:
        logger.debug(f'For partition_DWCS_H2O no lacking variables. Necessary variables were {h2o_dw_required_variables}.')
        try:
            ptt.partition_DWCS_H2O(str(dst_path), 
                                        NEE='NEE', GPP='GPP', Reco='Reco', CO2='wco2', 
                                        CO2neg_H2Opos='wco2-wh2o+', 
                                        CO2neg_H2Oneg='wco2-wh2o-', NIGHT=None)\
                                    .filter(['TIMESTAMP', 'NEE', 'GPP', 'Reco'])\
                                    .to_file(os.path.join(output_folderpath, str(site_name)+f'_CDWT_partitioning_H2O.csv'), index=False)
        except Exception as e:
            logger.warning(str(e))
    
    h2o_co_dw_required_variables = ['w','co2','h2o','co']
    is_lacking_variable = sum([v not in variables_available for v in h2o_co_dw_required_variables])
    if is_lacking_variable:
        logger.debug(f'For partition_DWCS_CO with {h2o_co_dw_required_variables} lacking variables.')
    if not is_lacking_variable:
        logger.debug(f'For partition_DWCS_CO no lacking variables. Necessary variables were {h2o_co_dw_required_variables}.')
        try:
            ptt.partition_DWCS_CO(str(dst_path), 
                                        NEE='NEE', GPP='GPP', Reco='Reco', ffCO2='ffCO2',
                                        CO2='wco2', 
                                        CO2neg_H2Opos='wco2-wh2o+', 
                                        CO2neg_H2Oneg='wco2-wh2o-', 
                                        CO2pos_COpos='wco2+wco+', 
                                        CO2pos_COneg='wco2+wco-',
                                        NIGHT=None)\
                                        .filter(['TIMESTAMP', 'NEE', 'GPP', 'Reco', 'ffCO2'])\
                                        .to_file(os.path.join(output_folderpath, str(site_name)+f'_CDWT_partitioning_H2O_CO.csv'), index=False)
        except Exception as e:
            logger.warning(str(e))
    
    co_dw_required_variables = ['w','co2','co']
    is_lacking_variable = sum([v not in variables_available for v in co_dw_required_variables])
    if is_lacking_variable:
        logger.debug(f'For partition_DWCS_CO with {co_dw_required_variables} lacking variables.')
    if not is_lacking_variable:
        logger.debug(f'For partition_DWCS_CO no lacking variables. Necessary variables were {co_dw_required_variables}.')
        try:
            ptt.partition_DWCS_CO(str(dst_path), 
                                        NEE='NEE', GPP='GPP', Reco='Reco', ffCO2='ffCO2',
                                        CO2='wco2', 
                                        CO2neg_H2Opos=['wco2-wco+', 'wco2-wco-'], 
                                        CO2neg_H2Oneg=None, 
                                        CO2pos_COpos='wco2+wco+', 
                                        CO2pos_COneg='wco2+wco-',
                                        NIGHT=None)\
                                        .filter(['TIMESTAMP', 'NEE', 'GPP', 'Reco', 'ffCO2'])\
                                        .to_file(os.path.join(output_folderpath, str(site_name)+f'_CDWT_partitioning_CO.csv'), index=False)
        except Exception as e:
            logger.warning(str(e))
        
    ch4_dw_required_variables = ['w','co2','ch4']
    is_lacking_variable = sum([v not in variables_available for v in ch4_dw_required_variables])
    if is_lacking_variable:
        logger.debug(f'For partition_DWCS_CO with {ch4_dw_required_variables} lacking variables.')
    if not is_lacking_variable:
        logger.debug(f'For partition_DWCS_CO no lacking variables. Necessary variables were {ch4_dw_required_variables}.')
        try:
            ptt.partition_DWCS_CO(str(dst_path), 
                                        NEE='NEE', GPP='GPP', Reco='Reco', ffCO2='ffCO2',
                                        CO2='wco2', 
                                        CO2neg_H2Opos=['wco2-wch4+', 'wco2-wch4-'], 
                                        CO2neg_H2Oneg=None, 
                                        CO2pos_COpos='wco2+wch4+', 
                                        CO2pos_COneg='wco2+wch4-',
                                        NIGHT=None)\
                                        .filter(['TIMESTAMP', 'NEE', 'GPP', 'Reco', 'ffCO2'])\
                                        .to_file(os.path.join(output_folderpath, str(site_name)+f'_CDWT_partitioning_CH4.csv'), index=False)
        except Exception as e:
            logger.warning(str(e))


def integrate_cospectra(data, f0, dst_path=None):
    logger = logging.getLogger('wvlt.pipeline.integrate_cospectra')
    logger.debug(f"Integrate cospectra with f0 = {f0}")
    
    data0 = data[(np.isnan(data['natural_frequency']) == False) * (data['natural_frequency'] >= f0)
                 ].groupby(['variable', 'TIMESTAMP'])['value'].agg("sum").reset_index(drop=False)
    data1 = data[np.isnan(data['natural_frequency'])].drop(
        'natural_frequency', axis=1)

    datai = pd.concat([data1[np.isin(
        data1['variable'], data0['variable'].unique()) == False], data0]).drop_duplicates()
    datai = datai.pivot_table('value', 'TIMESTAMP',
                              'variable').reset_index(drop=False)

    if dst_path:
        logger.debug(f"Writing cospectra with f0 = {f0} to file {dst_path}")
        datai.to_file(dst_path, index=False)
    return datai

def integrate_cospectra_from_file(root, f0, pattern='_full_cospectra_([0-9]+)_', 
                                  dst_path=None, newlog=False):
    """
    function: integrate cospectra from output files of process() (or main()) into a file.
    call: integrate_cospectra_from_file()
    Input:
        * root (str): Path to the folder with the files to be loaded. Usually the folder is named wavelet_full_cospectra. 
        * pattern (str, default '_full_cospectra_([0-9]+)_'): Pattern to be searched for in the files inside the folder. Usually they contain the pattern '_CDWT_full_cospectra_([0-9]{12})_'.
        * f0 (int, default None): Works as a high-pass filter for the wavelet cospectra (see similar process function f0 = 1/integration_period) inside integrate_cospectra().
        * newlog (bool, default False): if new log file in the subfolder log inside the output_folderpath is created using start_logging(). Useful if the function integrate_full_spectra_into_file() is called on its own, e.g. outside of eddypro_wavelet_run or with time delay after the function process().
        **kwargs
    Return:
        The integrated cospectrum. Also file saved accordingly.
    """
    
    # use glob.glob to find files matching the pattern
    logger = logging.getLogger('wvlt.pipeline.integrate_cospectra_from_file')
    logger.debug(f"Try to integrate cospectra from file in folder {root} with f0 = {f0}.")
    if isinstance(root, str):
        saved_files = {}
        for name in os.listdir(root):
            # logger.debug(f"Looking through {name}")
            dateparts = re.findall(pattern, name, flags=re.IGNORECASE)
            if len(dateparts) == 1:
                saved_files[dateparts[0]] = os.path.join(root, name)
        # logger.debug(f"Found {saved_files} that match the provided pattern: {pattern}.")

        def __read__(date, path):
            r = pd.read_csv(path, skiprows=11, sep=',')
            if 'natural_frequency' not in r.columns: 
                logger.warning(f'Skipping spectral file. Natural frequency column not found ({path}).')
                return pd.DataFrame()
            if r.natural_frequency.dtype != float: print(date, r.natural_frequency.dtype)
            r['TIMESTAMP'] = pd.to_datetime(date, format='%Y%m%d%H%M')
            return r

        data = pd.concat([__read__(k, v) for k, v in saved_files.items()])
    else:
        data = root
    
    return integrate_cospectra(data, f0, dst_path=dst_path)



def integrate_full_spectra_into_file(site_name, output_folderpath, 
                                     integration_period=60*30, newlog=False, **kwargs):
    """
    function: integrate cospectra from output files of process() (or main()) into a file.
    status: not necessary, please JUST USE integrate_cospectra_from_file instead.
    call: integrate_full_spectra_into_file()
    Input:
        * site_name: Site name of the data to be loaded in. Nessessary to construct file names to be loaded. File names have to correspond to output file names of the function main(), i.e. they need to contain the pattern '_CDWT_full_cospectra_([0-9]{12})_' and have to be saved inside the subfolder wavelet_full_cospectra of output_folderpath.
        * output_folderpath: Path of the folder where the output file gets saved.
        * integration_period (int, default None): integration period of the wavelength signal in s.
            Works as a high-pass filter for the wavelet cospectra (as f0 = 1/integration_period) inside integrate_cospectra().
        * newlog (bool, default False): if new log file in the subfolder log inside the output_folderpath is created using start_logging(). Useful if the function integrate_full_spectra_into_file() is called on its own, e.g. outside of eddypro_wavelet_run or with time delay after the function process().
        **kwargs
    Return:
        No return.
    """
    # CONCAT INTO SINGLE FILE
    
    # activate new logging file? Useful if function is called on its own, e.g. outside of eddypro_wavelet_run and with time delay after the process().
    if (output_folderpath is not None) and newlog:
        hc24.start_logging(output_folderpath)
    logger = logging.getLogger('wvlt.handler.integrate_full_spectra_into_file')
    
    dst_path = os.path.join(output_folderpath + str(site_name) + f"_CDWT_fulldata_integrated_{integration_period//60}min" + ".csv")
    
    logger.info(f'Integrate files with _CDWT_full_cospectra_([0-9]{12})_ inside {output_folderpath} into {dst_path}.')
    
    #dst_path = os.path.join(output_folderpath, str(
    #    site_name)+f'_CDWT_full_cospectra.csv')
    integrate_cospectra_from_file(os.path.join(output_folderpath, 'wavelet_full_cospectra'),
                                          1/integration_period, '_CDWT_full_cospectra_([0-9]{12})_', dst_path)


def decompose_variables(data, variables=['w', 'co2'], method='dwt', 
                     nan_tolerance=.3, identifier='0000', 
                     **kwargs):
    """    Calculate data decomposed with wavelet transform
    """
    # dictionary of wavelet transforms
    φ = {'info_names': []}
    sj = None

    # run by couple of variables (e.g.: co2*w -> mean(co2'w'))
    try:
        logger.debug(f"Debug variables are: {variables} and data columns are: {data.columns.tolist()}")
        assert len(
            variables), 'Empty list of covariances to run. Check available variables and covariances to be performed.'
        
        for var in variables:
            if var not in φ.keys():
                ready_signal = prepare_signal(
                    data[var], nan_tolerance=nan_tolerance, identifier=identifier)
                logger.debug(f"signal is ready: {ready_signal.signal.shape}")
                wt_signal = universal_wt(
                    signal=ready_signal.signal, method=method, **kwargs)
                logger.debug(
                    f"wt_signal is ready: {wt_signal.wave.shape}, {wt_signal.sj}")
                # φ[var], sj
                φ[var] = wt_signal.wave
                # COI
                coi = wt_signal.coi
                if coi is None: # coi is None for dwt
                    coi = np.ones_like(φ[var]) # produce same shape array as if there would be a coi, all with values 1. --> for dwt for quality control only nan amount.
                φ[f'{var}_qc'] = np.where(ready_signal.signan, 0, coi)
                #φ[f'{var}_qc'] = np.where(ready_signal.signan, 0, wt_signal.coi)
                φ['info_names'] += [var, f'{var}_qc']
                logger.debug(f"wt_signal is done.")
        φ.update({'sj': wt_signal.sj})
        φ.update({'coi': wt_signal.coi})
    except Exception as e:
        logger.critical(e)
        logger.error(f"Error in decompose_variables: {e}")
    # return φ
    return type('var_', (object,), φ)


def decompose_data(data, variables=['w', 'co2'], dt=0.05, method='dwt', 
                   nan_tolerance=.3, memory_eff=True, verbosity=1, **kwargs):
    """
    function: Calculate data decomposed with wavelet transform
    call: decompose_data()
    
    Input:
        * data (pandas.DataFrame): The data to be decomposed. Column TIMESTAMP necessary and variables as column names.
        * variables (list, default ['w', 'co2']): variables to decompose. Need to correspond to column names in data.
        * dt (float, default 0.05): Time step between the observations/measurements. Necessary to calculate frequency from wavelet scales, and fs, f1.
        * method (str, default 'dwt'): Wavelet decomposition method. Possible methods are 'dwt', 'cwt', 'fcwt'.
        * nan_tolerance (float, default .3): Passed to decompose_variables() and then prepare_signal(). Proportion or absolute value of NAN values allowed within the data. Passed to prepare_signal(). If too much NAN, warning gets called. Otherwise NAN get linear interpolated.
        * memory_eff (bool, default True): if False, fast but memory-heavy algorithm is used to combine all decomposed data. Otherwise memory-light but slow algorithm is used.
        **kwargs
    
    Return:
        Pandas DataFrame with columns: TIMESTAMP, natural_frequency, the decomposed variables, and their quality control (_qc).
    """
    logger = logging.getLogger('wvlt.pipeline.decompose_data')
    assert method in [
        'dwt', 'cwt', 'fcwt'], "Method not found. Available methods are: dwt, cwt, fcwt"
    
    # φ = {}
    # run by couple of variables (e.g.: co2*w -> mean(co2'w'))
    info_t_startvarloop = time.time()
    
    kwargs['fs'] = kwargs.get('fs', 1/dt) # if fs specified in wt_kwargs in process() its used here, otherwise taking 1/dt as fs.
    logger.debug(f'fs set to {kwargs['fs']}.')
    kwargs['f1'] = kwargs.get('f1', (1/2)*(1/dt)) # f1: lowest scale (2x sampling rate), if acquisition_frequency = 20 --> we can highest capture the 10 Hz frequencies because of sampling theorem.
    logger.debug(f'f1 set to {kwargs['f1']}')
    
    
    # run wavelet transform
    φ = decompose_variables(data, variables=variables, method=method,
                              nan_tolerance=nan_tolerance, **kwargs)
    # sj = φ.pop('sj', None)
    # for v in variables:
    #     if v not in φ.keys():
    #         signal = np.array(data[v])
    #         signan = np.isnan(signal)
    #         N = len(signal)
    #         Nnan = np.sum(signan)
    #         if Nnan:
    #             if (nan_tolerance > 1 and Nnan > nan_tolerance) or (Nnan > nan_tolerance * N):
    #                 logger.warning(
    #                     f"UserWarning: Too much nans ({np.sum(signan)}, {np.round(100*np.sum(signan)/len(signal), 1)}%) in {data.TIMESTAMP.head(-1)[0]}.")
    #         if Nnan and Nnan < N:
    #             signal = np.interp(np.linspace(0, 1, N), 
    #                     np.linspace(0, 1, N)[signan == False],
    #                     signal[signan==False])
    #         φ[v], sj = universal_wt(signal, method, **wt_kwargs)
    #         N = len(signal)
        

    logger.debug(f'\t\tDecompose all variables took {round(time.time() - info_t_startvarloop)} s.')
    
    φs_names = []
    # φs_names = [f'valid_{l}' if l else 'valid' for l in φ.sj]
    logger.debug(f'\t\tφ.info_names ({φ.info_names}).')
    for n in φ.info_names:
        if vars(φ)[n].shape[0] > 1:
            for l in φ.sj:  # use '+ ['']' if __integrate_decomposedarrays_in_dictionary__
                if l: φs_names += [f'{n}_{l}'] 
                else: φs_names += [n] 
        else: φs_names += [n]

    logger.debug(f"\t\tvars(φ).keys(): {vars(φ).keys()}")
    logger.debug(f'\t\tφs_names: {φs_names} ({len(φs_names)}).')

    # transform 2D arrays to DataFrame
    values = [vars(φ)[n] for n in φ.info_names]
    logger.debug(f'\t\tTransform 2D arrays to DataFrame with columns `{"`; `".join(φs_names)}`.')
    logger.debug(f'\t\t{[np.array(v).shape for v in values]}.')
    
    info_t_mattime = time.time()
    __temp__ = hc24.matrixtotimetable(np.array(data.TIMESTAMP),
                                      np.concatenate(values, axis=0), columns=φs_names)
    logger.debug(f'\t\t\tMatrix to timetable took {round(time.time() - info_t_mattime)} s.')
    
    __temp__.set_index('TIMESTAMP', inplace=True)
    if verbosity > 1: logger.debug(f'\t\tpd.MultiIndex.from_tuples: {[tuple(c.split("_")) for c in __temp__.columns]}.')
    
    if verbosity > 1: logger.debug(f'Calculating natural_frequency from scales for {__temp__.head()}.')
    __temp__.columns = pd.MultiIndex.from_tuples([tuple(c.rsplit('_', 1)) for c in __temp__.columns])
    if verbosity > 1: logger.debug(f'Calculating natural_frequency from scales for MultiIndex {__temp__.head()}.')
    
    if memory_eff==False: # old processing: memory heavy, but fast
        t0 = time.perf_counter()
        __temp_l__ = __temp__.stack(1).reset_index(1).rename(columns={"level_1": "natural_frequency"}).reset_index(drop=False) # very memory heavy. Try to replace
        elapsed = time.perf_counter() - t0
        logger.debug(f"Fast, but memory heavy processing took {elapsed:.3f} seconds")
    else: 
    # new processing: memory efficient but slow, enables using longer processing_time_durations or lower f0 --> more levels for dwt --> far more data loaded in, especially using necessary buffer to avoid edge effects or coi.
        t0 = time.perf_counter()
        # Initialize empty DataFrame, gets filled during loop step by step for each variable
        __temp_l__ = pd.DataFrame(index=pd.MultiIndex.from_product([__temp__.index, __temp__.columns.get_level_values(1).unique()], names=["TIMESTAMP", "natural_frequency"]))
        for var in __temp__.columns.get_level_values(0).unique():
            logger.debug(f'Looping over {var}.')
            df_var = __temp__[var]
            __temp__.drop(columns=var, level=0, inplace=True) # delete directly after using it to save RAM space
            s = df_var.stack()
            del df_var # delete directly after using it to save RAM space
            s.index = s.index.set_names(["TIMESTAMP", "natural_frequency"])
            __temp_l__[var] = s
            del s
        logger.debug(f'Looping worked, column names are {__temp_l__.columns}, now resetting index.')
        __temp_l__.reset_index(inplace=True)
        elapsed = time.perf_counter() - t0
        logger.debug(f"Slow, but memory efficient processing took {elapsed:.3f} seconds")
    
    # Making sure the new method produces the same output as the old:
    #__temp_l__.to_csv("../test_outputs/Newmethod.csv")
    #__temp2__.to_csv('../test_outputs/Oldmethod.csv')
    # --> the data is the same but ordered differently.
    # --> when the output of the process function is compared (already averaged), the output is exactly the same between old and new method.

    #pattern = re.compile(r"^(?P<variable>.+?)_?(?P<natural_frequency>(?<=_)\d+)?$")
    #__temp__ = __temp__.melt(['TIMESTAMP'] + φ.keys())
    #__temp__ = pd.concat([__temp__.pop('variable').str.extract(pattern, expand=True), __temp__], axis=1)
    __temp_l__['natural_frequency'] = __temp_l__['natural_frequency'].apply(lambda j: 1/hc24.j2sj(j, 1/dt) if j else np.nan)
    if verbosity > 1: logger.debug(f'Calculating natural_frequency finished. Returning DataFrame {__temp_l__.head()}')
    logger.debug(f'\t\tDecompose data (full) took {round(time.time() - info_t_startvarloop)} s.')
    return __temp_l__


def _calculate_product_from_formula_(data, formula='w*co2|w*h2o'):
    logger = logging.getLogger('wvlt.pipeline._calculate_product_from_formula_')

    φs = {}
    logger.debug(
        f"\t\tformula: {formula} ({type(formula)}).")

    formulavar = formula_to_vars(formula) if isinstance(
        formula, str) else formula
    xy_name = ''.join(formulavar.xy)
    if xy_name not in data.columns:
        for ci, c in enumerate(formulavar.xy):
            XY = XY * np.array(data[c]).conjugate() if ci else data[c]
        φs[xy_name] = XY
        logger.debug(
            f"\t\tDecomposed covariance shape: {XY.shape}, named: {xy_name}.")
    # ({round(XY.shape[1] / (24*60*60/20), 2)} days, for dt=20Hz)
                        
    # logger.debug(f"\t\tformulavar.condsamp_pair: {formulavar.condsamp_pair}.")
    # logger.debug(f"\t\tφs: {φs}.")
    # φ = pd.DataFrame(φs)
    # logger.debug(f"\t\tφ: {φ.head()}.")
    for cs in formulavar.condsamp_pair:
        cs_name = ''.join(cs)
        # logger.debug(f"\t\tcs_name: {cs_name} (from {cs}) | φs.keys(): {φs.keys()}.")
        if (cs_name not in φs.keys()) and (cs_name not in data.columns):
            for ci, c in enumerate(cs):
                # logger.debug(f"\t\tCurrent c: {c}.")
                # logger.debug(f"\t\tCurrent data: {data.head()}.")
                # logger.debug(f"\t\tDecomposed covariance shape: {data[c].shape}.")
                CS = CS * np.array(data[c]).conjugate() if ci else data[c]
            φs[cs_name] = CS

    # data = pd.concat([data, pd.DataFrame(φs)], axis=1)
    return pd.DataFrame(φs)


def _calculate_conditional_sampling_from_formula_(data, formula='w*co2|w*h2o', cond_samp_both=False):
    """
    function: calculate from formula what arguments passed to conditional_sampling().
    
    call: _calculate_conditional_sampling_from_formula_()
    
    Input:
        * data (pandas.DataFrame): Data to be conditionally sampled. Need to contain columns fitting to the variables in the formula. If the formula contains 'w*co2|w*h2o', using the function formula_to_vars() it gets to the names ['wh2o', 'wco2']. Those names need to be in the column names of the dataframe.
        * formula (str, default 'w*co2|w*h2o'): Which variables are how to be conditionally sampled. "w*co2|w*h2o" means: conditionally sample w*co2 depending on w*co2 and w*h2o.
        * cond_samp_both (bool, default True): If True both parts of the formula are conditionally sampled. If False, only the leading part of the formula is sampled. E.g. if False in case of 'w*co2|w*h2o', we get the output columns wco2+wh2o+,wco2-wh2o+,wco2+wh2o-,wco2-wh2o-, stating wco2 being conditionally sampled e.g. for wco2+wh2o+ when wco2 is positiv AND wh2o is positive. If True, we get the output columns wco2+wh2o+,wco2+wh2o-,wco2-wh2o+,wco2-wh2o-,wh2o+wco2+,wh2o+wco2-,wh2o-wco2+,wh2o-wco2-, hence, we get both, wco2 and wh2o conditionally sampled.
    
    Return:
        pd.DataFrame(φc) (pandas.DataFrame): Dataframe containing the conditionally sampled data.
    """
    logger = logging.getLogger('wvlt.pipeline._calculate_conditional_sampling_from_formula_')

    logger.debug(f"\t\tformula: {formula}.")
    formulavar = formula_to_vars(formula) if isinstance(formula, str) else formula
    
    logger.debug(f"\t\tformulavar: {formulavar}.")
    names = [''.join(formulavar.xy)] + [''.join(cs)
                                        for cs in formulavar.condsamp_pair]

    logger.debug(f"\t\tnames: {names} ({data.columns.to_list()}).")

    if cond_samp_both==True:
        logger.debug(f"Conditional sampling for all parts of formula: {names}.")
        φc = {}
        for i, name in enumerate(names):
            logger.debug(f"Run conditional sampling for {name}.")
            names_i = names.copy()
            # shift the name that is being conditionally sampled to the front
            # that the column names created by conditional_sampling() fit.
            names_i.insert(0, names_i.pop(i))
            φc_a = ptt.conditional_sampling(
                np.array(data[name]), *[np.array(data[n]) for n in names_i], names=names_i, label={1: "+", -1: "-"}) if names_i else {}
            φc = {**φc, **φc_a}
            #logger.debug(f"φc is within the loop is {φc}")
    else:
        logger.debug(f"Conditional sampling only for the first part of formula: {names[0]}.")
        φc = ptt.conditional_sampling(
            np.array(data[names[0]]), *[np.array(data[n]) for n in names], names=names, label={1: "+", -1: "-"}) if names else {}

    # data = pd.concat([data, pd.DataFrame(φs)], axis=1)
    return pd.DataFrame(φc)


def __save_cospectra__(data, dst_path, overwrite=False, **meta):
    logger = logging.getLogger('wvlt.pipeline.__save_cospectra__')

    # saved_files = []

    info_t_startsaveloop = time.time()

    # for __datea__, __tempa__ in data.groupby(data.TIMESTAMP):
    # dst_path = output_path.format(pd.to_datetime(__datea__).strftime('%Y%m%d%H%M'))
    logger.debug(f'\t\tSaving {dst_path} with shape {data.shape}.')
    # if os.path.exists(dst_path): continue
    use_header = False

    if overwrite or (not os.path.exists(dst_path)):
        use_header = True
        header  = "wavelet_based_(co)spectra\n"
        header += f"--------------------------------------------------------------\n"
        header += f"TIMESTAMP_START = {meta.get('TIMESTAMP_START', min(data.TIMESTAMP))}\n"
        header += f"TIMESTAMP_END = {meta.get('TIMESTAMP_END', max(data.TIMESTAMP))}\n"
        header += f"N: {meta.get('N', len(data.TIMESTAMP))}\n"
        header += f"TIME_BUFFER [min] = {meta.get('buffer', np.nan)/60}\n"
        header += f"frequency [Hz]\n"
        header += f"y-axis -> wavelet_reconstructed\n"
        header += f"mother_wavelet -> {meta.get('method', '')}\n"
        header += f"acquisition_frequency [Hz] = {1/meta.get('dt', np.nan)}\n"
        header += f"averaging_interval = {meta.get('averaging', '')}\n"
        hc24.mkdirs(dst_path)
        with open(dst_path, 'w+') as part: part.write(header)
        # legitimate_to_write = 1
        logger.debug(f'\t\tSaving header of DataFrame took {round(time.time() - info_t_startsaveloop)} s.')
        # saved_files.append(dst_path)
    
    # if not legitimate_to_write: continue
    
    data.drop('TIMESTAMP', axis=1, inplace=True)
    with open(dst_path, 'a+', newline='') as part:
        data.to_file(part, header=use_header, chunksize=500, index=False)
    logger.debug(f'\t\tSaving DataFrame took {round(time.time() - info_t_startsaveloop)} s.')
    
    # del data
        
    #arr_slice = np.unique(data.TIMESTAMP, return_index=True)
    #for __datea__ in arr_slice[0]:
    #    dst_path = output_path.format(suffix, pd.to_datetime(__datea__).strftime('%Y%m%d%H%M'))
    #    if os.path.exists(dst_path+'.part'): os.rename(dst_path+'.part', dst_path)
    
    # return saved_files
    return

def cs_partition_NEE_ET(site_name, output_folderpath, NEE=True, ET=True, 
    integration_period=None, 
    variables_available=['h2o', 'wh2o+wco2-', 'wh2o-wco2-', 'wh2o-wco2+', 'wh2o+wco2+', 'co2', 'wco2-wh2o+', 'wco2-wh2o-'], 
        newlog=False):
    """
    function: Read avaraged and integrated cospectra file, with conditionally sampled fluxes, apply conditional sampling for NEE and ET.
    Input: 
        * site_name (str): Site name of the data to be loaded in. Nessessary to construct file names to be loaded. See variable output_folderpath for more information.
        * output_folderpath (str): Path to folder where the input and output files files are saved. Inside this folder there has to be a file with the pattern os.path.join(output_folderpath, f"{site_name}_CDWT_fulldata_integrated_*min.csv"). Usually produced by integrate_full_spectra_into_file() or by process().
        * NEE (bool, default True): If True, NEE is partitioned.
        * ET (bool, default True): If True, ET is partitioned.
        * integration_period (int, default None): If different files with different integration_period inside the output_folderpath, this helps to find the correct file for conditional sampling. In those functions it is the integration period of the wavelength signal in s. Works as a high-pass filter for the wavelet cospectra (as f0 = 1/integration_period) inside integrate_cospectra(). Also relevant for the filename of saved data. It gets constructed similar to os.path.join(output_folderpath, str(site_name)+f'_CDWT_partitioning_H2O.csv' dependent on the used partitioning algorithm.
        * variables_available (list, default ['h2o', 'wh2o+wco2-', 'wh2o-wco2-', 'wh2o-wco2+', 'wh2o+wco2+', 'co2', 'wco2-wh2o+', 'wco2-wh2o-']): From which variables are data available. Necessary to test, if partitioning algorithms can be run.
        * newlog (bool, default False): if new log file in the subfolder log inside the output_folderpath is created using start_logging(). Useful if the function condition_sampling_partition() is called on its own, e.g. outside of eddypro_wavelet_run() or with time delay after other functions.
    Return:
        dat (Pandas.DataFrame): Dataframe with the partitioned fluxes.
    """
    # function to load output file from integrate_full_spectra_into_file() or process()
    # and, partitiones the data and save its results as well in a new file.
    # inside process file or inside handler?
    logger = logging.getLogger('wvlt.pipeline.cs_partition_NEE_ET')
    # activate new logging file? Useful if function is called on its own, 
    # e.g. outside of eddypro_wavelet_run and with time delay after the process().
    if (output_folderpath is not None) and newlog:
        hc24.start_logging(output_folderpath)
        
    # to be able to have different integration_period = 1/f0, hence different high pass filters in the folder
    # search for the pattern with variable minutes
    if not integration_period:
        pattern = os.path.join(output_folderpath, f"{site_name}_CDWT_fulldata_integrated_*min.csv")
        matches = glob.glob(pattern)
        if not matches:
            raise FileNotFoundError(f"No file found matching {pattern}")
        if len(matches) > 1:
            logger.warning(f"Multiple files of integrated cospectra found: {matches}, taking the first one for partitioning! Set integration_period if want to use another one.")
        dst_path = matches[0]
        logger.info(f'Taking the file {dst_path} for partitioning.')
    else:
        dst_path = os.path.join(output_folderpath + str(site_name) + f"_CDWT_fulldata_integrated_{integration_period//60}min" + ".csv")
        logger.debug(f"Specified integration_period={integration_period}. Hence taking the file {dst_path}")
    
    list_dat = []
    
    # Partitioning ET
    if ET:
        logger.debug("Trying to partition ET")
        ETpartition_DWCS_required_variables = ['h2o', 'wh2o+wco2-', 'wh2o-wco2-', 'wh2o-wco2+', 'wh2o+wco2+']
        is_lacking_variable = sum([v not in variables_available for v in ETpartition_DWCS_required_variables])
        if is_lacking_variable:
            logger.warning(f'For ETpartition_DWCS with {ETpartition_DWCS_required_variables} lacking variables.')
        if not is_lacking_variable:
            logger.debug(f'For ETpartition_DWCS no lacking variables. Necessary variables were {ETpartition_DWCS_required_variables}.')
            dat_part = pttET.ETpartition_DWCS(str(dst_path))\
                        .filter(['TIMESTAMP', 'ET', 'T', 'E', 'Dew'])
            list_dat.append(dat_part)
    
    # Partitioning NEE
    if NEE:
        logger.debug("Starting to partition NEE")
        partition_DWCS_H2O_required_variables = ['co2', 'wco2-wh2o+', 'wco2-wh2o-']
        is_lacking_variable = sum([v not in variables_available for v in partition_DWCS_H2O_required_variables])
        if is_lacking_variable:
           logger.warning(f'For partition_DWCS_H2O with {partition_DWCS_H2O_required_variables} lacking variables.')
        if not is_lacking_variable:
            logger.debug(f'For partition_DWCS_H2O no lacking variables. Necessary variables were {partition_DWCS_H2O_required_variables}.')
            dat_part = ptt.partition_DWCS_H2O(str(dst_path), 
                                        NEE='NEE', GPP='GPP', Reco='Reco', CO2='wco2', 
                                        CO2neg_H2Opos='wco2-wh2o+', 
                                        CO2neg_H2Oneg='wco2-wh2o-', NIGHT=None)\
                                    .filter(['TIMESTAMP', 'NEE', 'GPP', 'Reco'])
            list_dat.append(dat_part)
            
    if list_dat:
        dat = list_dat[0]
    for df in list_dat[1:]:
        dat = dat.merge(df, on="TIMESTAMP", how="outer")
    dat.to_file(os.path.join(output_folderpath, str(site_name)+'_CDWT_partitioning_NEE_ET.csv'), index=False)
    
    return dat

def process(datetimerange, fileduration, input_path, acquisition_frequency,
            covariance=None, cond_samp_both=False, output_folderpath=None,
            overwrite=False, processing_time_duration="1d",
            integration_period=None, partition=None,
            method="dwt", average_period="30min", sitename="00000",
            wt_kwargs={}, meta={}, **kwargs):
    """
    function: process data. (1) gets data, (2) performs wavelet transform, (3) cross calculate variables using conditional_sampling, (4) averages, (5) saves. Implemented as loops to prevent RAM overflow.
    
    call: process()
    
    Input:
        * datetimerange (str): date time range from which the data is processed. Format: YYYYMMDDHHMM-YYYYMMDDHHMM
        * fileduration (int): time range that the input files cover in minutes (e.g. 30). Needed inside bufferforfrequency_dwt() for only taking steps in file-size to calculate buffer size. Also necessary to calculate correct file names from bmmflux inside universal_reader().
        * input_path (str): path to the folder where the input files are located.
        * acquisition_frequency (int): frequency of the data in Hz. Used as dt = 1/acquisition_frequency. Used to calculate the sampling frequency fs for wavelet decomposition inside of decompose_data() and then passed to universal_wt().
        * covariance (list, default: None): variables to be considered in the calculations as strings in a list. * denotes the covariance. | denotes conditional sampling. Format: e.g. ["w*co2|w*h2o"]. In this example, "w*co2|w*h2o" means: conditionally sample w*co2 depending on w*co2 and w*h2o. This produces the columns wco2+wh2o+,wco2-wh2o+,wco2+wh2o-,wco2-wh2o-, which mean e.g. in the case wco2+wh2o+ that wco2 is sampled when wco2 is positive AND wh2o is positive.
        * cond_samp_both (bool, default True): If True both parts of the formula are conditionally sampled. If False, only the leading part of the formula is sampled. E.g. if False in case of 'w*co2|w*h2o', we get the output columns wco2+wh2o+,wco2-wh2o+,wco2+wh2o-,wco2-wh2o-, stating wco2 being conditionally sampled e.g. for wco2+wh2o+ when wco2 is positiv AND wh2o is positive. If True, we get the output columns wco2+wh2o+,wco2+wh2o-,wco2-wh2o+,wco2-wh2o-,wh2o+wco2+,wh2o+wco2-,wh2o-wco2+,wh2o-wco2-, hence, we get both, wco2 and wh2o conditionally sampled.
        * output_folderpath (str, default: None): path to the folder where the output is saved.
        * overwrite (bool, default False): if files can be overriden. If True, output files not get overriden and no calculation is performed for these data.
        * processing_time_duration (str, default "1d"): Time duration over which the calculation is perfomed in a loop. Important setting to prevent overflowing of RAM. Format: pandas time offset string, e.g. "3h". Possible specifications are s, min, h, d.
        * integration_period (int, default None): integration period of the wavelength signal in s. Works as a high-pass filter for the wavelet cospectra (as f0 = 1/integration_period) inside integrate_cospectra().
        * partition (list, default None): Gives if ET or NEE should be partitioned. Set as strings in a list, e.g. ["ET", "NEE"], or in case only NEE: ["NEE"]. Necessary to set an integration_period for this.
        * method (str, default "dwt"): One of 'dwt', 'cwt', 'fcwt', passed as kwargs to the functions main() and decompose_data().
        * average_period (str, default '30min'): Averaging period for averaging the wavelet decompositioned values. Format: pandas time string, e.g. "30min". Possible specifications are s, min, h, d. Passed to the main function.
        * sitename (str, default "00000"): Sitename, files get named accordingly.
        * wt_kwargs (dict, default {}): **kwargs passed to the wavelet tranformation itself. Can e.g. include wavelet specification. See wavelet_function.py for more details. Important setting include f0, which is the lowest frequency for wavelet decomposition. Because adding buffer is necessary to prevent edge effects, a lower f0 drastically increases the amount of data loaded in and the memory usage.
        * meta (dict, default {}): Header lines in the output files. Get filled successively during the code run.
        **kwargs
        
    Return:
        fulldata (pandas.DataFrame): Containing all processed data. If integration_period is specified already integrated.
    
    """
    logger = logging.getLogger('wvlt.pipeline.process')
    local_args = locals()
    
    info_t_start = time.time()

    def _date_from_yl(date):
        date = re.sub('[-: ]', '', str(date))
        if processing_time_duration.endswith("d"):
            date = date[:8]
        if processing_time_duration.endswith("h") or processing_time_duration.endswith("min"):
            date = date[:12]
        return date
    
    def _validate_run(date, yl, compare_start=True, compare_end=False):
        # recheck if files exist and overwrite option
        # doesn't save time (maybe only save 5min)
        file_name = os.path.basename(output_pathmodel.format(date))
        part_name0 = file_name.rsplit('_', 1)[0] + '_' if compare_start else ''
        part_name1 = file_name.rsplit('.', 1)[-1] if compare_end else ''
        current_files = [p for p in os.listdir(os.path.dirname(
            output_pathmodel)) if p.startswith(part_name0) and p.endswith(part_name1)]

        if not overwrite and current_files: #file_name in os.path.exists(output_pathmodel.format(date)):
            logger.warning(
                "UserWarning: Skipping, file already exists ({}).".format(date))
            return False

        if hc24.checkifinprogress(curoutpath_inprog):
            return False
        return True
        
    def _load_data():
        start_time = time.time()
        data = loaddatawithbuffer(yl, d1=None, freq=_f, buffer=buffer, f_freq=_f, **load_kwargs)
        if data.empty:
            logger.warning(f"UserWarning: No file found ({date}, path: {load_kwargs.get('path', 'default')}).")
            return None
        logger.debug(f'\tLoading data took {round(time.time() - start_time)} s.')
        return data
    
    def _exit():
        if os.path.exists(curoutpath_inprog):
            os.remove(curoutpath_inprog)

    # Group parameters for each function
    load_kwargs = {
        # 'datetimerange': datetimerange,
        'fileduration': fileduration,
        'path': input_path,
        # 'acquisition_frequency': acquisition_frequency,
        'fkwargs': {'dt': 1/acquisition_frequency},
        'fmt': kwargs.get("fmt", {}),  # allow user to override or extend
        **kwargs.get("load_kwargs", {}),  # allow user to override or extend,
    }

    if 'gas4_name' in kwargs.keys():
        load_kwargs['fmt'].update({kwargs.pop('gas4_name'): '4th gas'})

    transform_kwargs = {
        'dt': 1/acquisition_frequency,
        'method': method,
        'average_period': average_period,
        'varstorun': covariance or hc24.available_combinations(hc24.DEFAULT_COVARIANCE),
        **kwargs.get("transform_kwargs", {}),
        **wt_kwargs
    }

    output_kwargs = {
        'output_folderpath': output_folderpath,
        'overwrite': overwrite,
        'meta': {'acquisition_frequency': acquisition_frequency},
        **kwargs.get("output_kwargs", {})
    }
    
    logging_kwargs = {
        **kwargs.get("logging_kwargs", {})
    }

    run_time = datetime.datetime.now().strftime("%Y%m%d%H%M")
    output_path = ""
    output_pathmodel = ""
    curoutpath_inprog = ""

    # kwargs.update({'covariance': covariance,
    #           'processduration': processduration, })
    
    if output_folderpath is not None:
        hc24.start_logging(output_folderpath, **logging_kwargs)
        if method == 'cov':
            output_pathmodel = str(os.path.join(output_folderpath, str(
                sitename)+'_cov_{}_'+run_time+'.csv'))
        else:
            output_pathmodel = str(os.path.join(output_folderpath, 'wavelet_full_cospectra', str(
                sitename)+'_CDWT_full_cospectra_{}_'+run_time+'.csv'))
        hc24.mkdirs(output_pathmodel)
        try:
            hc24.save_locals(local_args, os.path.join(output_folderpath, f'log/setup_{run_time}.yml'))
        except Exception as e:
            logger.warning(f"Could not save local arguments to file: {e}")
        # with open(, 'w+') as stp:
        #     yaml.safe_dump(local_args, stp)
    else:
        output_pathmodel = ""
    
    logger.debug(f'Output path: {output_pathmodel}.')
    
    ymd = [datetimerange.split(
        '-')[0], datetimerange.split('-')[1], f'{fileduration}min']
    
    fulldata = pd.DataFrame()
    info_t_start = time.time()
    logger.info(f'In load_main.')

    print(f'\nRUNNING WAVELET TRANSFORM ({method})')
    
    _, _, _f = ymd
    ymd = hc24.list_time_in_period(*ymd, processing_time_duration, include='both')

    buffer = 0 if method == 'cov' else (
        bufferforfrequency_dwt(n_=_f, **wt_kwargs) / 2 if method == 'dwt'
        else bufferforfrequency(wt_kwargs.get("f0", 1/(3*60*60))) / 2
    )
    meta.update({'buffer': buffer})
    logger.debug(f"Buffer (s): {buffer}.")
    logger.info(f'Start date loop at {round(time.time() - info_t_start)} s (load_main).')

    # Skip two line
    prev_print = '\n'
    print(f"{prev_print} Wavelet decompositon seperated into {processing_time_duration} chunks of data.{prev_print}")
    
    prev_print = '\x1B[1A\r'
    for yl in tqdm.tqdm(ymd): # for yl in ymd:
        info_t_yl_ymd = time.time()
        date = _date_from_yl(yl[0])
        
        print(prev_print, 'Current chunk from', yl[0], 'to', yl[-1])
        
        # print(prev_print, date, 'reading', ' '*10, sep=' ', end='\n')
        #prev_print = '\x1B[1A\r'

        if output_pathmodel:
            curoutpath_inprog = output_pathmodel.format(date).rsplit(".", 1)[
                0] + ".inprogress"
            logger.debug(f'In progress file: {curoutpath_inprog}.')
            if not _validate_run(date, yl):
                continue
        
        data = _load_data()
        if data is None:
            _exit()
            continue
        

        try:
            # main run
            output_path = output_pathmodel.format("{}")
            # # run directly
            # output_kwargs.update({'output_path': output_path})
            # fulldata = main(data, period=[min(yl), max(yl)],
            #                 output_kwargs=output_kwargs, **transform_kwargs)
            
            # run by varstorun
            output_kwargs.update({'output_path': output_path + '.part'})
            allvars = transform_kwargs['varstorun']
            logger.debug(f"Allvars that get looped through: {allvars}")
            saved_files = []
            for f in allvars:
                f_transform_kwargs = transform_kwargs.copy()
                f_transform_kwargs['varstorun'] = [f]
                logger.debug(f"Available data columns: {data.columns.tolist()}")
                logger.debug(f"Processing formula: {f}, varstorun passed to main(): {transform_kwargs['varstorun']}")
                output = main(data, period=[min(yl), max(yl)], 
                              cond_samp_both=cond_samp_both,
                              meta=meta,
                              output_kwargs=output_kwargs, **f_transform_kwargs)
                saved_files.append(output.saved)
                fulldata = pd.concat([fulldata, output.data], axis=0)
                
            for f in [s for s_ in saved_files for s in s_]:
                if os.path.exists(f):
                    os.rename(f, f.replace('.part', ''))
        except Exception as e:
            logger.critical(e)
            print(str(e))
        logger.debug(f'Date loop ({yl}) took {round(time.time() - info_t_yl_ymd)} s.')

        _exit()
    
    logger.debug(f'End date loop at {round(time.time() - info_t_start)} s.')
    logger.debug(f"integration_period: {integration_period}.")
    logger.debug(f"output_pathmodel: {output_pathmodel}.")
    logger.debug(f"fulldata: {fulldata.head()}.")
    if output_pathmodel and not fulldata.empty:
        
        # Integrating
        #dst_path = os.path.join(output_folderpath, os.path.basename(
        #    output_pathmodel.format(run_time)))
        dst_path = os.path.join(output_folderpath + str(sitename) + f"_CDWT_fulldata" + ".csv")
        if integration_period:
            print(f'\n \n Integration of all frequencies until a period of {integration_period} s.')
            logger.debug(f'Running Integration of wavelet flux with {integration_period}')
            dst_path = os.path.join(output_folderpath + str(sitename) + f"_CDWT_fulldata_integrated_{integration_period//60}min" + ".csv")
            fulldata = integrate_cospectra(fulldata, 1/integration_period, dst_path=None)
        fulldata.to_csv(dst_path, index=False)
        
        # Partitioning
        if partition:
            if integration_period:
                print(f'\n Partitioning of integrated wavelet flux.')
                if "NEE" in partition:
                    NEE = True
                if "ET" in partition:
                    ET = True
                av_var = fulldata.columns
                if not NEE and not ET:
                    logger.warning("Neiter ET nor NEE set to partition but wanted to partion. This is not possible.")
                else:
                    logger.debug(f'Running partitioning of ET:{ET} and NEE:{NEE} with available variables {av_var}')
                    fulldata = cs_partition_NEE_ET(site_name=sitename, 
                                        output_folderpath=output_folderpath,
                                        NEE=NEE,
                                        ET=ET,
                                        integration_period=integration_period,
                                        variables_available=av_var)
            else:
                logger.warning("Integration period not set but wanted to partion. This is not possible.")
    logger.debug(f'\t\tFull process took {round(time.time() - info_t_start)} s (run_wt).')
    return fulldata


def main(data, varstorun, period=None, average_period='30min', 
         cond_samp_both=False, output_kwargs={}, meta={}, **kwargs):
    """
    function: Performs wavelet transform for specified variables, cross calculate variables using conditional_sampling and averages. Can save the data in files.
    call: main()
    Input:
        * data (pandas.DataFrame): Data to be processed.
        * varstorun (list): variables to be considered in the calculations as strings in a list. * denotes the covariance. | denotes conditional sampling. Format: e.g. ["w*co2|w*h2o"]
        * period (list, default None): List with two entries. Decomposed signal only used for data['TIMESTAMP'] > period[0]) & data['TIMESTAMP'] < period[1].
        * average_period (str, default '30min'): Averaging period for averaging the wavelet decompositioned values. Format: pandas time string, e.g. "30min". Possible specifications are s, min, h, d.
        * output_kwargs (dict, default {}): Specify output variables. For saving the data, output_path needs to be set as string containing an element {0} to paste the data in, e.g. output_kwargs={'output_path':'../test_outputs/test_{0}.csv'}. Possible further specification is overwrite (bool) specifiying if files can get overwritten.
        * meta (dict, default {}): Header lines in the output files. Get filled successively during the code run.
        **kwargs
    Return:
        A new class object named var_ with class attributes data and saved. Data includes the averaged wavelet transformed, cross calculated variables. saved_files contains strings with paths to where the saved files are placed. If save return as test = main(), access data via test.data or test.saved.
    """
    logger = logging.getLogger('wvlt.pipeline.main')
    logger.info('In main.')
    logger.debug(f'Input data shape: {data.shape}.')
    info_t_main = time.time()

    vars_unique = list(set([var for f in varstorun for var in formula_to_vars(f).uniquevars]))

    logger.debug(f'Input data is ready, data shape is {data.shape}, with unique vars: {vars_unique}, based on {"; ".join(varstorun)}.')
    
    
    # decompose all required variables
    wvvar = decompose_data(data, vars_unique,
                           nan_tolerance=.3,
                           identifier='0000',
                           **kwargs)
    meta.update({'averaging': average_period,
                 'method': f"{kwargs.get('method', '')} ~{kwargs.get('mother_wavelet', '')}",
                 'dt': kwargs.get('dt', np.nan)})
    
    logger.debug(f'Decompose data is over, data shape is {wvvar.shape}.')
    logger.debug(f'\n{wvvar.head()}\n')
    
    # select valid dates
    logger.debug(f'Period of interest is from {period[0]} to {period[1]}')
    if period: wvvar = wvvar[(wvvar['TIMESTAMP'] > period[0]) & (wvvar['TIMESTAMP'] < period[1])]
    wvvar = wvvar.reset_index(drop=True)

    logger.debug(f'Screen data over period of interest yielded data shape {wvvar.shape}.')

    # calculate covariance
    info_t_calc_product = time.time()
    # wvout = _calculate_product_from_formula_(wvvar, varstorun)
    logger.debug(f'varstorun. {varstorun}')
    uniquecovs = list(set(
        [c for f in varstorun for c in formula_to_vars(f).combinations]))
    logger.debug(f'uniquecovs. {uniquecovs}')

    wvout = pd.concat(
        # [wvvar[['TIMESTAMP', 'natural_frequency']]] +
        [(_calculate_product_from_formula_(wvvar, formula=f)
        #   .drop(columns=formula_to_vars(f).uniquevars if i == 0 else ['TIMESTAMP'] + formula_to_vars(f).uniquevars)
          )
         for i, f in enumerate(uniquecovs)], axis=1)
    
    logger.debug(
        f'\tCalculate product from formula took {round(time.time() - info_t_calc_product)} s.')

    growingdata = pd.concat([wvvar, wvout], axis=1)
    logger.debug(f'Growing data shape {growingdata.shape}.')

    # calculate conditional sampling
    info_t_calc_cond_samp = time.time()
    logger.debug(f'Starting _calculate_conditional_sampling_from_formula_.')
    wvcsp = pd.concat(
        # [wvvar[['TIMESTAMP', 'natural_frequency']]] +
        [_calculate_conditional_sampling_from_formula_(growingdata, f, cond_samp_both)
         for f in varstorun], axis=1)
         
    logger.debug(f'\tCalculate conditional sampling took {round(time.time() - info_t_calc_cond_samp)} s.')

    # despike
    # denoise
    # smoothing
    # Y12 = smooth_2d_data(Y12, method='repeat', smoothing=smoothing)
    # for i in range(len(φcs)):
    #     φcs[i] = smooth_2d_data(
    #         φcs[i], method='repeat', smoothing=smoothing)

    # assemble data
    info_t_assemble_data = time.time()
    growingdata = pd.concat([growingdata, wvcsp], axis=1)
    logger.debug(f'\tAssemble data took {round(time.time() - info_t_assemble_data)} s.')

    # average
    info_t_average = time.time()
    for thisdate, thisdata in growingdata.groupby(growingdata['TIMESTAMP'].dt.floor(
            average_period)):
        thisdate_ = thisdate.strftime('%Y%m%d%H%M')
        meta.update({thisdate_: meta})
        meta[thisdate_].update({
                    'TIMESTAMP_START': min(thisdata['TIMESTAMP']),
                    'TIMESTAMP_END': max(thisdata['TIMESTAMP']),
                    'N': len(thisdata['TIMESTAMP'])})

    growingdata['TIMESTAMP'] = growingdata['TIMESTAMP'].dt.floor(
        average_period)
    __ID_COLS__ = list(
        set(['TIMESTAMP', 'natural_frequency']) & set(growingdata.columns))
    growingdata = growingdata.groupby(__ID_COLS__).agg(
        "mean").reset_index(drop=False)
    logger.debug(f'\tAveraging data took {round(time.time() - info_t_average)} s.')


    # save in dataframe and .csv
    growingdata = (growingdata.sort_values(by=__ID_COLS__)
                .melt(__ID_COLS__))
    
    saved_files = []
    if output_kwargs.get('output_path', None):
        logger.debug(f"\tSaving data in {output_kwargs['output_path']}.")
        info_t_save_cospectra = time.time()
        for thisdate, thisdata in growingdata.groupby(growingdata.TIMESTAMP):
            thisdate_ = thisdate.strftime('%Y%m%d%H%M')
            dst_path = output_kwargs.get('output_path').format(thisdate_)
            overwrite = output_kwargs.get('overwrite', False)
            logger.debug(f"\t\t\t... in {dst_path}.")
            __save_cospectra__(thisdata, dst_path, overwrite, **meta[thisdate_])
            saved_files += [dst_path]
        logger.debug(f'\tSaving data took {round(time.time() - info_t_save_cospectra)} s.')

    
    # rename saved files when done
    # os.rename(output_kwargs['output_path'].format(suffix, pd.datetime.now().strftime('%Y%m%dT%H%M%S_%f')),
    #           )
    
    logger.debug(f'\tMain took {round(time.time() - info_t_main)} s.')
    # save in .nc
    return type('var_', (object,), {'data': growingdata, 'saved': saved_files})



def run_from_eddypro(path,
                     # ="input/EP/FR-Gri_sample.eddypro",
                     #  covariance=["w*co2|w|co2|h2o", "w*co2|w*h2o", "w*h2o",],
                     #  processduration='6H',
                     **kwargs):
    c = eddypro.extract_info_from_eddypro_setup(eddypro=path)
    c.update(**kwargs)

    for path in ['input_path', 'output_folderpath']:
        if c.get(path, None) is not None:
            c[path] = os.path.abspath(c[path])

    return process(**c)


