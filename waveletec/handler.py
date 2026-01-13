# standard modules
import os
import sys
import re
import logging
import datetime
import glob
# 3rd party modules
import yaml
# Project modules
from . import _core as hc24
from .partitioning import coimbra_et_al_2025 as ptt
#from .partitioning import schoendorf as pttET
from ._core import wavelet_functions as wavelet_functions
from ._extra import eddypro_tools as eddypro
from ._core import pipeline


logger = logging.getLogger('wvlt.handler')

def sample_raw_data(input_path, datetimerange, acquisition_frequency=20, fileduration=30, **kwargs):
    raw_kwargs = {'path': input_path, 'fkwargs': {'dt': 1/acquisition_frequency}}
    kwargs['fmt'] = kwargs.get('fmt', {})
    if 'gas4_name' in kwargs.keys(): kwargs['fmt'].update({kwargs.pop('gas4_name'): '4th gas'})
    raw_kwargs.update({k: v for k, v in kwargs.items() if k in ['fmt']})

    ymd = [datetimerange.split('-')[0], datetimerange.split('-')[1], f'{fileduration}min']
    _, _, _f = ymd
    ymd = hc24.list_time_in_period(*ymd, '1D', include='both')

    for ymd_i, yl in enumerate(ymd):
        data = hc24.loaddatawithbuffer(
            yl, d1=None, freq=_f, buffer=0, f_freq=_f, **raw_kwargs)
        break
    return data


def sample_raw_data(input_path, datetimerange, acquisition_frequency=20, fileduration=30, processduration='1D'):
    ymd = [datetimerange.split(
        '-')[0], datetimerange.split('-')[1], f'{fileduration}min']
    _, _, _f = ymd
    ymd = waveletec._core.list_time_in_period(
        *ymd, processduration, include='both')

    for ymd_i, yl in enumerate(ymd):
        data = waveletec._core.loaddatawithbuffer(
            yl, d1=None, freq=_f, buffer=0, f_freq=_f, **{'path': input_path, 'fkwargs': {'dt': 1/acquisition_frequency}})
        break
    return data

def run_from_eddypro(path,
                     #="input/EP/FR-Gri_sample.eddypro",
                    #  covariance=["w*co2|w|co2|h2o", "w*co2|w*h2o", "w*h2o",],
                    #  processduration='6H',
                    **kwargs):
    c = eddypro.extract_info_from_eddypro_setup(eddypro=path)
    c.update(**kwargs)

    for path in ['input_path', 'output_folderpath']:
        if c.get(path, None) is not None:
            c[path] = os.path.abspath(c[path])

    return pipeline.process(**c)

# raw_kwargs = {'path': input_path, 'fkwargs': {'dt': 1/acquisition_frequency}}
# kwargs['fmt'] = kwargs.get('fmt', {})
# if 'gas4_name' in kwargs.keys(): kwargs['fmt'].update({kwargs.pop('gas4_name'): '4th gas'})
# raw_kwargs.update({k: v for k, v in kwargs.items() if k in ['fmt']})

# ymd, raw_kwargs, output_folderpath = None, verbosity = 1,
# overwrite = False, processing_time_duration = "1D",
# internal_averaging = None, dt = 0.05, wt_kwargs = {},
# method = "dwt", averaging = 30, **kwargs)

    

def eddypro_wavelet_run(site_name, input_path, outputpath, datetimerange, acquisition_frequency=20, fileduration=30, 
         processduration='1D', integration_period=None, preaverage=None,
         covariance = None, variables_available=['u', 'v', 'w', 'ts', 'co2', 'h2o'], denoise=0, deadband=[], 
         method = 'dwt', wave_mother='db6', **kwargs):
    local_args = locals()

    if outputpath is not None:
        hc24.start_logging(outputpath)

        # Select output file path
        if method == 'cov':
            outputpath = str(os.path.join(outputpath, str(site_name)+'{}_{}.csv'))
        else:
            outputpath = str(os.path.join(outputpath, 'wavelet_full_cospectra', str(site_name)+'_CDWT{}_{}.csv'))

        # Save args for run
        hc24.mkdirs(outputpath)
        with open(os.path.join(os.path.dirname(os.path.dirname(outputpath)), f'log/setup_{datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")}.yml'), 'w+') as stp:
            yaml.safe_dump(local_args, stp)

    # Select covariances
    # x*y → Cov(x, y)
    # x*y|x*z|x*... → Cov(x, y)|Cov(x, z),Cov(x, ...)
    if covariance is None:
        covariance = hc24.available_combinations(
            hc24.DEFAULT_COVARIANCE, variables_available)

    # RUN WAVELET FLUX PROCESSING
    # ymd = [START_DATE, END_DATE, FILE_FREQUENCY]
    raw_kwargs = {'path': input_path, 'fkwargs': {'dt': 1/acquisition_frequency}}
    kwargs['fmt'] = kwargs.get('fmt', {})
    if 'gas4_name' in kwargs.keys(): kwargs['fmt'].update({kwargs.pop('gas4_name'): '4th gas'})
    raw_kwargs.update({k: v for k, v in kwargs.items() if k in ['fmt']})
    data = wavelet_functions.load_data_and_loop(ymd = [datetimerange.split('-')[0], datetimerange.split('-')[1], f'{fileduration}min'],
                                         output_path = outputpath,
                                         varstorun = covariance,
                                         averaging = [fileduration],
                                         processing_time_duration = processduration,
                                         method = method,
                                         wt_kwargs = {'fs': acquisition_frequency, 'wavelet': wave_mother},
                                         raw_kwargs = raw_kwargs,
                                         verbosity=5)
    return data


def integrate_full_spectra_into_file(site_name, output_folderpath, integration_period=60*30, newlog=False, **kwargs):
    """
    function: integrate cospectra from output files of process() (or main()) into a file.
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
    pipeline.integrate_cospectra_from_file(os.path.join(output_folderpath, 'wavelet_full_cospectra'),
                                          1/integration_period, '_CDWT_full_cospectra_([0-9]{12})_', dst_path)
    
    #pipeline.integrate_cospectra_from_file(os.path.join(output_folderpath, 'wavelet_full_cospectra'),
    #                                      1/integration_period, '_CDWT_full_cospectra_([0-9]{12}).csv$', dst_path)
    #hc24.concat_into_single_file(
    #    os.path.join(outputpath, 'wavelet_full_cospectra'), str(site_name)+f'_CDWT_full_cospectra.+.{fileduration}mn.csv', 
    #    output_path=dst_path, skiprows=10)


def condition_sampling_partition(site_name, output_folderpath, integration_period=None, variables_available=['w', 'co2', 'h2o'], newlog=False, **kwargs):
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

def hndl_cs_partition_NEE_ET(site_name, output_folderpath, NEE=True, ET=True, 
    integration_period=None, 
    variables_available=['h2o', 'wh2o+wco2-', 'wh2o-wco2-', 'wh2o-wco2+', 'wh2o+wco2+', 'co2', 'wco2-wh2o+', 'wco2-wh2o-'], 
        newlog=False):
    pipeline.cs_partition_NEE_ET(site_name=site_name, output_folderpath=output_folderpath,
                        NEE=NEE, ET=ET, integration_period=integration_period,
                        variables_available=variables_available, newlog=newlog)
    

def set_cwd(workingdir=None):
    """
    Set the current working directory to the specified path.
    
    Parameters:
    - workingdir: Path to the directory to set as current working directory.
    """
    if workingdir is not None:
        os.chdir(workingdir)
        print(f"Current working directory set to: {os.getcwd()}")
        logger.info(f"Current working directory set to: {os.getcwd()}")
    else:
        logger.debug("No working directory specified, using current directory.")


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-ep', '--eddypro',   type=str,
                        help='Path to EddyPro setup file')
    parser.add_argument('-m', '--metadata',   type=str,
                        help='Path to EddyPro metadata file')
    parser.add_argument('-s', '--site_name',   type=str,
                        help='Site name for the processing')
    parser.add_argument('-i', '--input_path',  type=str,
                        help='Path to the input data folder')
    parser.add_argument('-o', '--output_folderpath', type=str,
                        help='Path to the output folder where results will be saved')
    parser.add_argument('-d', '--datetimerange', type=str,
                        help='Date and time range for processing, format: YYYYMMDD-HHMM-YYYYMMDD-HHMM')
    parser.add_argument('-af', '--acquisition_frequency', type=int, default=20,
                        help='Acquisition frequency in Hz (default: 20 Hz)')
    parser.add_argument('-fd', '--fileduration', type=int, default=30,
                        help='File duration in minutes (default: 30 minutes)')
    parser.add_argument('-ip', '--integration_period', type=int, default=None,
                        help='Integration period in seconds (default: None)')
    parser.add_argument('-v', '--variables_available', type=str, nargs='+',
                        help='List of available variables (default: ["u", "v", "w", "ts", "co2", "h2o"])')
    # parser.add_argument('-dk', '--despike', type=int)  # , nargs=1)
    # parser.add_argument('-dn', '--denoise', type=int)  # , nargs=1)
    # parser.add_argument('-db', '--deadband', type=str, nargs='+')
    parser.add_argument('-cov', '--covariance', type=str, nargs='+',
                        help='List of covariances to compute (default: None)')
    parser.add_argument('--method', type=str, default='dwt', choices=['cov', 'dwt', 'cwt'],
                        help='Method to use for wavelet processing (default: "dwt")')
    parser.add_argument('--wave_mother', type=str, default='db6',
                        help='Mother wavelet to use for wavelet processing (default: "db6")')
    parser.add_argument('--run', type=int, default=1,
                        help='Run the wavelet processing (default: 1)')
    parser.add_argument('--concat', type=int, default=1,
                        help='Concatenate results into a single file (default: 1)')
    parser.add_argument('--partition', type=int, default=1,
                        help='Run partitioning on the results (default: 1)')
    parser.add_argument('--processing_time_duration', type=str, default='1D',
                        help='Processing time duration for the wavelet processing (default: "1D")')
    # parser.add_argument('--preaverage', type=str, default=None)
    parser.add_argument('-cwd', '-workingdir', type=str, default=None,
                        help='Set the current working directory (default: None)')
    parser.add_argument('--log_level', type=str, default='INFO', choices=[
                        'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help='Set the logging level')

    args = parser.parse_args()
    return vars(args)


def prepare_args(args):
    # Set the current working directory if specified
    if args.get('cwd'):
        set_cwd(os.chdir(args['cwd']))

    # Set logging level
    log_level = getattr(logging, args.pop(
        'log_level', 'INFO').upper(), logging.INFO)
    args['logging_kwargs'] = {'level': log_level}

    # Extract and pop run flags
    run = args.pop('run', 1)
    concat = args.pop('concat', 1)
    partition = args.pop('partition', 1)
    ep_setup = args.pop('eddypro')
    ep_meta = args.pop('metadata')

    # Retrieve eddypro setup
    exta = eddypro.extract_info_from_eddypro_setup(ep_setup, ep_meta)
    args = eddypro.update_args_with_extracted_info(args, exta)

    # Default integration period
    logger.debug(
        f"Integration period {args['integration_period']}, setting to {args['fileduration']} minutes.")
    if args['integration_period'] is None:
        args['integration_period'] = args['fileduration'] * 60

    # Sanitize site name
    args['site_name'] = args['site_name'].replace('/', '_').replace('\\', '_')

    # Convert paths
    for path in ['input_path', 'output_folderpath', 'eddypro', 'metadata']:
        if args.get(path):
            args[path] = os.path.abspath(args[path])

    # Return args and control flags
    return args, run, concat, partition, ep_setup


def validate_args(args):
    required = ['site_name', 'input_path', 'output_folderpath',
                'datetimerange', 'acquisition_frequency', 'fileduration']
    missing = [f'`{k}`' for k in required if args.get(k, None) is None]
    if missing:
        raise ValueError(f'Missing argument(s): {", ".join(missing)}')


def log_args(args):
    # print('\n'.join(
    #     [f'{k}:\t{v[:5] + "~" + v[-25:] if isinstance(v, str) and len(v) > 30 else v}' for k, v in args.items()]), end='\n\n')

    print('Start run with:')
    for k, v in args.items():
        if isinstance(v, str) and len(v) > 30:
            v = f"{v[:5]}~{v[-25:]}"
        print(f"{k}:\t{v}")
    print()


def main(**args):
    args, run, concat, partition, ep_setup = prepare_args(args)
    validate_args(args)
    log_args(args)

    if args['method'] == 'cov':
        concat = partition = False

    if run:
        run_from_eddypro(ep_setup, **args)
    # if concat:
    #     integrate_full_spectra_into_file(**args)
    if partition:
        condition_sampling_partition(**args)


if __name__ == '__main__':
    args = parse_args()
    main(**args)
