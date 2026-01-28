"""
This script is a key part of the following publications:
    - Herig Coimbra, Pedro Henrique and Loubet, Benjamin and Laurent, Olivier and Mauder, Matthias and Heinesch, Bernard and 
    Bitton, Jonathan and Delpierre, Nicolas and Depuydt, Jérémie and Buysse, Pauline, Improvement of Co2 Flux Quality Through 
    Wavelet-Based Eddy Covariance: A New Method for Partitioning Respiration and Photosynthesis. 
    Available at SSRN: https://ssrn.com/abstract=4642939 or http://dx.doi.org/10.2139/ssrn.4642939

Modified by Daniel Schöndorf, University of Bayreuth.
"""

##########################################
###     IMPORTS                           
##########################################

# standard modules
import copy
import os
import re
import warnings
import logging
import time
import datetime
from functools import reduce

# 3rd party modules
import yaml
import numpy as np
import pandas as pd
from itertools import chain
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression, RANSACRegressor, HuberRegressor, QuantileRegressor
import zipfile
from io import StringIO
from pandas.api.types import is_numeric_dtype, is_object_dtype

# project modules
from .commons import update_nested_dicts, structuredData, yaml_to_dict
from .addons import *

##########################################
###     PROJECT CHOICES                           
##########################################

SITES_TO_STUDY = ['SAC']

DEFAULT_FILE_RAW = {
    'file_pattern': '([0-9]{8}-[0-9]{4})_raw_dataset_.*.txt', 
    'date_format': '%Y%m%d-%H%M', 
    'dt': 0.05, 
    'tname': "TIMESTAMP", 
    'id': None,
    'datefomatfrom': '%Y%m%d%H%M%S.%f', 
    'datefomatto': '%Y-%m-%dT%H:%M:%S.%f'
}

DEFAULT_READ_CSV = {
    'handle_eddypro_raw_dataset': True,
    'skiprows': 8,
    'sep': "\\s+",
    'na_values': ['NaN', 'nan', -9999],
}

DEFAULT_READ_GHG = {
    'skiprows': 7,
    'sep': r"\t",
    'engine': 'python'
}

DEFAULT_FMT_DATA = {
    '4thgas': '4th gas',
}


class structuredDataFrame:
    def __init__(self, data=None, dt=None, **kwargs):
        if data is None:
            loopvar = kwargs.pop('lookup', [])
            loopvar = [l.to_list() if isinstance(l, list)==False else l for l in loopvar]
            for l in loopvar:
                result = universal_reader(lookup=l, **kwargs, verbosity=0)
                self.__dict__.update(result.__dict__)
        
        else:
            assert dt is not None, 'Missing measurement frequency (dt).'
            self.data = data
            self.dt = dt
            self.__dict__.update(**kwargs)
    
    def filter(self, items: dict):
        for k, v in items.items():
            if isinstance(v, tuple):
                self.data = self.data.loc[(self.data[k] > v[0])
                                          & (self.data[k] < v[1])].copy()
            else:
                self.data = self.data[self.data[k].isin(v)].copy()
        return self

    def rename(self, names: dict):
        self.data = self.data.rename(columns=names)
        return self

    def modify(self, items: dict):
        for k, v in items.items():
            self.data[k] = v
        return self

    def format(self, 
               cols={'t':'ts'}, 
               keepcols=['u', 'v', 'w', 'ts', 'co2', 'co2_dry', 'h2o', 'h2o_dry', 'ch4', 'n2o'],
               addkeep=[],
               colsfunc=str.lower, cut=False, **kwargs):
        
        if isinstance(self, pd.DataFrame):
            formated = self
        else:
            fmt_clas = structuredDataFrame(**self.__dict__)
            formated = fmt_clas.data

        if colsfunc is not None:
            formated.columns = map(colsfunc, formated.columns)
        #cols.update(kwargs)
        cols.update({v.lower(): k.lower() for k, v in kwargs.items() if isinstance(v, list)==False})
        cols = {v: k for k, v in {v: k for k, v in cols.items()}.items()}
        cols.update({'timestamp': 'TIMESTAMP'})
        #formated.TIMESTAMP = formated.TIMESTAMP.apply(np.datetime64)
        if cut:
            #formated = formated[[
            #    c for c in formated.columns if c in cols.keys()]]
            formated = formated.loc[:, np.isin(formated.columns, keepcols+addkeep+list(cols.keys()))]
        
        formated = formated.rename(columns=cols)

        if isinstance(self, pd.DataFrame):
            return formated
        else:
            fmt_clas.data = formated
            return fmt_clas
    
    def interpolate(self, cols=["co2", "w"], qcname="qc"):
        interpolated = structuredDataFrame(**self.__dict__)
        interpolated.data[qcname] = 0
        for c_ in list(cols):
            interpolated.data[qcname] = interpolated.data[qcname] + 0 * \
                np.array(interpolated.data[c_])
            interpolated.data.loc[np.isnan(interpolated.data[qcname]), qcname] = 1
            interpolated.data[qcname] = interpolated.data[qcname].astype(int)
            interpolated.data[c_] = interpolated.data[c_].interpolate(method='pad')
            
        return interpolated


def __universal_reader__(path, **kw_csv):
    logger = logging.getLogger('waveletec.read_data.__universal_reader__')
    
    handle_eddypro_raw_dataset = kw_csv.pop('handle_eddypro_raw_dataset', False)
    handle_bmmflux_raw_dataset = kw_csv.pop('handle_bmmflux_raw_dataset', False)
    if handle_eddypro_raw_dataset and handle_bmmflux_raw_dataset:
        logger.critical("Specified handle_eddypro_raw_dataset and handle_bmmflux_raw_dataset both True. This should not be possible!")
    
    if path.endswith('.gz'): kw_csv.update(**{'compression': 'gzip'})
    elif path.endswith('.csv'): kw_csv.pop('compression', None)
    if path.endswith('.ghg'):
        with zipfile.ZipFile(path, 'r') as zip_ref:
            datafile = [zip_ref.read(name) for name in zip_ref.namelist() if name.endswith(".data")][0]
        datafile = str(datafile, 'utf-8')
        path = StringIO(datafile)
        # DEFAULT_READ_GHG
        kw_csv.update(DEFAULT_READ_GHG)
        
    try:
        if handle_eddypro_raw_dataset:
            with open(path, 'r') as file: header = [c.replace('\n', '').strip() for c in file.readlines()[9].split('  ') if c]
            kw_csv.update({'skiprows': 10, 'sep': '\\s+', 'na_values': ['NaN', 'nan', -9999], 'names': header})
        elif handle_bmmflux_raw_dataset:
            #logger.debug('Inside handle_bmmflux_raw_dataset')
            with open(path, 'r') as file: header = file.readlines()[0].strip().split(',')
            kw_csv.update({'skiprows': 2, 'sep': ',', 'na_values': ['NaN', 'nan', -9999], 'names': header})
        df_td = pd.read_csv(path, **kw_csv)
    except Exception as e:
        # (EOFError, pd.errors.ParserError, pd.errors.EmptyDataError):
        try:
            #if verbosity>1: warnings.warn(f'{e}, when opening {path}, using {kw_csv}. Re-trying using python as engine and ignoring bad lines.')
            df_td = pd.read_csv(path, on_bad_lines='warn', engine='python', **kw_csv)
        except Exception as ee:
            warnings.warn(f'{ee}, when opening {str(path)}, using {kw_csv}')
            return None
    return df_td

def universal_reader(path, lookup=[], fill=False, fmt={}, onlynumeric=True, verbosity=1, fkwargs={}, tipfile="readme.txt", **kwargs):
    """
    path: str or dict 
    e.g. (str): path = str(os.path.join(inputpath, 'eddypro_raw_datasets/level_6'))
    e.g. (dict): path = {k: str(os.path.join(v, 'eddypro_raw_datasets/level_6')) for k, v in inputpath.items()}
        where k will become the suffix in case of duplicates
    """
    logger = logging.getLogger('waveletec.read_data.universal_reader')
    
    
    if isinstance(path, dict):
        dataframes = {}
        for k, v in path.items():
            assert isinstance(v, str), f'Path ({v}) is not string.'
            path_ = universal_reader(v, lookup, fill, fmt, onlynumeric, verbosity, fkwargs, tipfile, **kwargs)
            dt = path_.dt
            dataframes[k] = path_.data
        #df_site = dataframes.pop(k, pd.DataFrame())
        dup_keys = {}
        for k in dataframes.keys():
            dup_keys[k] = [set(dataframes[k].columns).intersection(set(v_.columns)) - set(['TIMESTAMP']) for k_, v_ in dataframes.items() if k_!=k]
            dup_keys[k] = list(chain.from_iterable(dup_keys[k]))
            dataframes[k].rename(columns={c: c + k for c in dup_keys[k]}, inplace=True)
        
        try:
            df_site = reduce(lambda left, right: pd.merge(left, right, on='TIMESTAMP', how='outer', suffixes=('', '_DUP')), 
                            list(dataframes.values()))
        except Exception as e:
            logger.error(str(e))
            logger.debug('Concatenating instead of merging.')
            #df_site = reduce(lambda left, right: pd.concat([left, right], axis=1), dataframes)
            return structuredDataFrame(pd.DataFrame(), dt=dt)
        #for k, v in dataframes.items():
        #    df_site = df_site.merge(v, on='TIMESTAMP', how='outer', suffixes=('', k))
        return structuredDataFrame(df_site, dt=dt)
    
    df_site = pd.DataFrame()
    
    folders = [path + p + '/' for p in os.listdir(path) if os.path.isdir(path + p)]
    folders = folders if folders else [path]
    
    if verbosity > 1: logger.debug(f"Check readable files in {path}")
    #print("Check readable files in", path if len(path)<40 else f'{path[:5]}...{path[-30:]}')#, fmt, fkwargs, kwargs)

    for path_ in folders:
        if verbosity > 1: logger.debug(f"Read files from {path_}")
        df_td = pd.DataFrame()

        # read tips file
        if verbosity > 1: logger.debug(f"Trying to read tipfile from {os.path.join(path, tipfile)}")
        kw_ = update_nested_dicts({"FILE_RAW": DEFAULT_FILE_RAW, "READ_CSV": DEFAULT_READ_CSV, "FMT_DATA": DEFAULT_FMT_DATA}, 
                                  os.path.join(path, tipfile), os.path.join(path_, tipfile),
                                  {"FILE_RAW": fkwargs, "READ_CSV": kwargs, "FMT_DATA": fmt},
                                  fstr=lambda d: yaml_to_dict(d))
        if verbosity > 1: logger.debug(f"kw_ is {kw_}")
        kw = structuredData(**kw_['FILE_RAW'])
        
        kw_csv = kw_['READ_CSV']
        if kw_csv != DEFAULT_READ_CSV: kw_csv['handle_eddypro_raw_dataset'] = False
        
        if verbosity > 1: logger.debug(f"lookup before handle_bmmflux_raw_dataset is {lookup}")
        if kw_csv.get('handle_bmmflux_raw_dataset'):
            logger.info('handle_bmmflux_raw_dataset was True, therefore resetting the matching file_pattern and the matching date format.')
            kw.file_pattern = "([0-9]{8}_[0-9]{6}).*?_3Drot_frc_tempplausadj\\.csv"
            kw.date_format = '%Y%m%d_%H%M%S'
            # need to redefine lookup, because bmmflux files always have the
            # mid of the file period in their file name not the beginning
            fileduration = kw_csv.get('fileduration')
            lookup = [f + datetime.timedelta(minutes=(fileduration/2)) for f in lookup]
        
        try:
            if ('header_file' in kw_csv.keys()) and (os.path.exists(kw_csv['header_file'])):
                kw_csv['header_file'] = "[" + open(kw_csv['header_file']).readlines()[0].replace("\n", "") + "]"
        except:
            None
        
        if verbosity > 1: logger.debug(f"lookup is {lookup}")
        lookup_ = list(set([f.strftime(kw.date_format) for f in lookup]))
        if verbosity > 1:
            logger.debug(f"lookup_ is {lookup_}")
            
        files_list = {}

        for root, directories, files in os.walk(path_):
            for name in files:
                dateparts = re.findall(kw.file_pattern, name, flags=re.IGNORECASE)
                if len(dateparts) == 1:
                    files_list[dateparts[0]] = os.path.join(root, name)
        
        if verbosity > 1: logger.debug(f"Five entries to files_list.keys() with the found dateparts are {list(files_list.keys())[0:5]}")
        
        #logger.debug(f"Found {len(files_list.keys())} files. One of them is {list(files_list.keys())[1]}.")
        
        # removing fileduration from the dict if it exists in there, because not needed inside __universal_reader__
        kw_csv.pop('fileduration', None) 
        # if lookup input existed (list with datetimes), 
        # we only load those files that their name pattern has a matching
        # time interval 
        # --> see the finding of dateparts in the file names above.
        # These dateparts make up the dictionary files_list
        # and hence, they need to match the lookup_, which is lookup but converted
        # to strings following the specified kw.date_format.
        for td in set(lookup_) & files_list.keys() if lookup_ != [] else files_list.keys():
            path_to_tdfile = files_list[td]
            if os.path.exists(path_to_tdfile):
                if verbosity > 1: logger.debug(f"Reading file {path_to_tdfile}")
                df_td = __universal_reader__(path_to_tdfile, **kw_csv)
                if verbosity > 1: logger.debug(f'File read, header is {df_td.head()}.')
                
                if kw.datefomatfrom == 'drop':
                    df_td = df_td.rename({kw.tname: kw.tname+'_orig'})
                
                if kw.tname not in df_td.columns or kw.datefomatfrom == 'drop':
                    if verbosity > 1: logger.debug(f"Not found your Timestamp column name {kw.tname} or datefomatfrom=='drop'. Creating new column {kw.tname}.")
                    if "date" in df_td.columns and "time" in df_td.columns:
                        df_td[kw.tname] = pd.to_datetime(
                            df_td.date + " " + df_td.time, format='%Y-%m-%d %H:%M')
                    elif any(col in df_td.columns for col in ['Year', 'Month', 'Day', 'Hour', 'Minute', 'Second', 'Millisecond']):
                        # for bmmflux files
                        present_col = [col for col in ['Year', 'Month', 'Day', 'Hour', 'Minute', 'Second', 'Millisecond'] if col in df_td.columns]
                        if verbosity > 1:
                            logger.debug(f"Creating new column {kw.tname} from the columns {present_col}.")
                        df_td[kw.tname] = pd.to_datetime(df_td[present_col])
                    else:
                        if verbosity > 1: logger.debug("Did not found data or time columns. Creating Timestamp using specified dt.")
                        df_td[kw.tname] = pd.to_datetime(
                            td, format=kw.date_format) - datetime.timedelta(seconds=kw.dt) * (len(df_td)-1 + -1*df_td.index)
                            #td, format=kw.date_format) + datetime.timedelta(seconds=kw.dt) * (df_td.index)
                        df_td[kw.tname] = df_td[kw.tname].dt.strftime(
                            kw.datefomatto)
                else:
                    if verbosity > 1: logger.debug("Found your specified timestamp column {kw.tname}. Now testing if numeric or if necessary to transform.")
                    try:
                        if is_numeric_dtype(df_td[kw.tname]):
                            df_td.loc[:, kw.tname] = df_td.loc[:, kw.tname].apply(lambda e: pd.to_datetime('%.2f' % e, format=kw.datefomatfrom).strftime(kw.datefomatto))
                        elif is_object_dtype(df_td[kw.tname]):
                            df_td.loc[:, kw.tname] = df_td.loc[:, kw.tname].apply(lambda e: pd.to_datetime(e).strftime(kw.datefomatto))
                        else:
                            df_td.loc[:, kw.tname] = pd.to_datetime(df_td[kw.tname], format=kw.datefomatfrom).strftime(kw.datefomatto)
                    except:
                        warnings.warn(f'error when converting {kw.tname} from {kw.datefomatfrom} to {kw.datefomatto}.')
                        continue
                
                df_td['file'] = td
                #df_site = df_site.append(df_td)
                df_site = pd.concat([df_site, df_td], ignore_index=True).reset_index(drop=True)
        
        if df_td.empty == False:
            break
        
    #print('df_td.empty ', df_td.empty)
    if onlynumeric:
        valcols = [i for i in df_site.columns if i.lower() not in [kw.tname.lower(), 'file']]
        _bf = df_site.dtypes
        #df_site.loc[:, valcols] = df_site.loc[:, valcols].apply(pd.to_numeric, errors='coerce')
        df_site[valcols] = df_site[valcols].apply(pd.to_numeric, errors='coerce')
        _af = df_site.dtypes
        if verbosity>1:
            _bfaf = []
            for (k, b) in _bf.items():
                if b!=_af[k]:
                    _nonnum = [s for s in np.unique(df_site[k].apply(lambda s: str(s) if re.findall('[A-z/]+', str(s)) else '')) if s]
                    _bfaf += ['{}, changed from {} to {}. ({})'.format(k, b, _af[k], ', '.join(_nonnum) if _nonnum else 'All numeric')]
            if _bfaf:
                warnings.warn(', '.join(_bfaf))
    
    #if kw_fmt:
    df_site = structuredDataFrame.format(df_site, **kw_['FMT_DATA'])

    if fill:
        if lookup:
            minmax = [min(lookup), max(lookup)]
        else:
            minmax = [np.nanmin(df_site[kw.tname]),
                      np.nanmax(df_site[kw.tname])]
        df_site = df_site.set_index(kw.tname).join(pd.DataFrame({kw.tname: pd.date_range(*minmax, freq=str(kw.dt*1000) + 'ms')}).set_index(kw.tname),
                how='outer').ffill().reset_index()
        #if 'co2' in df_site.columns and (abs(np.max(df_site.co2)) < 1000) and (abs(np.min(df_site.co2)) < 1000):
        #    df_site.loc[:, "co2"] = df_site.loc[:, "co2"] * 1000  # mmol/m3 -> μmol/m3
    
        
    if kw.id is not None:
        return {kw.id: structuredDataFrame(df_site, dt=kw.dt)}
    else:
        return structuredDataFrame(df_site, dt=kw.dt)


def loaddatawithbuffer(d0, d1=None, freq=None, buffer=None, 
                       tname="TIMESTAMP", f_freq=30, sort=True, **kwargs):
    """
    Load data with buffer.
    d0: start date
    d1: end date
    freq: frequency of the data
    buffer: buffer in seconds
    tname: name of the timestamp column
    f_freq: frequency of the data in minutes
    kwargs: other arguments to pass to the function"""
    if isinstance(d0, (pd.DatetimeIndex, list, set, tuple)) and d1 is None:
        d0, d1 = [np.nanmin(d0), np.nanmax(d0)]
    
    if buffer == None:
        datarange = [pd.date_range(start=d0, end=d1, freq=f"{f_freq}min")[:-1] + pd.Timedelta(freq)]
    else:
        # buffer align with file frequency (e.g. 30 min)
        freqno = int(re.match(r"\d*", f"{f_freq}min")[0])
        bufi = np.ceil(buffer / (freqno*60)) * freqno
        datarange = [
            pd.date_range(
                start=pd.to_datetime(d0) - pd.Timedelta(bufi, unit='min'),
                end=pd.to_datetime(d1) + pd.Timedelta(bufi, unit='min'),
                freq=freq)[:-1] + pd.Timedelta(freq)]
                
    if not datarange:
        return pd.DataFrame()
    
    data = structuredDataFrame(lookup=datarange, **kwargs)
    if data == None or data.data.empty:
        return data.data
    data.data[tname] = pd.to_datetime(data.data[tname])

    if sort:
        data.data = data.data.sort_values(by=tname).reset_index(drop=True)
    
    if buffer:
        d0 = pd.to_datetime(d0) - pd.Timedelta(buffer*1.1, unit='s')
        d1 = pd.to_datetime(d1) + pd.Timedelta(buffer*1.1, unit='s')
        data.filter({tname: (d0, d1)})

    # garantee all data points, if any valid time, else empty dataframe
    if np.sum(np.isnat(data.data.TIMESTAMP)==False):
        #data.data = pd.merge(pd.DataFrame({tname: pd.date_range(*nanminmax(data.data.TIMESTAMP), freq="0.05S")}),
        #                    data.data,
        #                    on=tname, how='outer').reset_index(drop=True)
        return data.data
    else:
        pd.DataFrame()

