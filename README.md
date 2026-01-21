# WaveletEC Module - Version NEE_ET

WaveletEC is a Python module designed for processing and analyzing Eddy Covariance data using wavelet transforms. This module provides functionalities to handle EddyPro setup files, perform wavelet analysis, and integrate results.

The module is forked from [Pedro H H Coimbra](https://github.com/pedrohenriquecoimbra/wavelet-ec) which was published in

[![DOI](https://zenodo.org/badge/DOI/10.1016/j.agrformet.2025.110684.svg)](https://doi.org/10.1016/j.agrformet.2025.110684)

[![DOI](https://zenodo.org/badge/786866970.svg)](https://zenodo.org/doi/10.5281/zenodo.11071327)

Pedro H H Coimbra, Benjamin Loubet, Olivier Laurent, Matthias Mauder, Bernard Heinesch, Jonathan Bitton, Nicolas Delpierre, Daniel Berveiller, Jérémie Depuydt, Pauline Buysse. Evaluation of a novel approach to partitioning respiration and photosynthesis using eddy covariance, wavelets and conditional sampling. https://doi.org/10.1016/j.agrformet.2025.110684

\* corresponding author: pedro-henrique.herig-coimbra@inrae.fr

## Version NEE_ET
This fork, developed by Daniel Schöndorf, contains some additions:
- the possibility to partition ET in addition to NEE (in development)
- reading bmmflux high-frequency corrected output files
- the option to run a more memory-efficient but slower algorithm (in the code at the moment just acivated. Need to be changed in the code in ```decompose_data```, variable ```memory_eff=True```)
- some minor bug fixes
- additional documentation about several functions

## Installation

### From git

To install the WaveletEC module, clone the repository and install the required dependencies:

```bash
git clone https://github.com/Biotit/wavelet-ec
cd waveletec
```
(Optional) Using conda to run savely in environment use:
 ```
 conda create -f conda/environment.yml
 conda activate wavec2
 ```
 
Inside the folder setup.py:
```bash
pip install
```


## Usage

### Running _waveletec_ in general

This requires ```pip install waveletec```.
```python
import waveletec
```
Then
```python
waveletec.process(datetimerange, fileduration, input_path, acquisition_frequency,
            covariance=None, cond_samp_both=False, output_folderpath=None,
            overwrite=False, processing_time_duration="1d",
            integration_period=None, partition=None,
            method="dwt", average_period="30min", sitename="00000",
            wt_kwargs={}, meta={}, **kwargs)
```
The documentation of process is:
``` python
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
        * partition (list, default None): Gives if ET and/or NEE should be partitioned. Set as strings in a list, e.g. ["ET", "NEE"], or in case only NEE: ["NEE"]. Necessary to set an integration_period for this.
        * method (str, default "dwt"): One of 'dwt', 'cwt', 'fcwt', passed as kwargs to the functions main() and decompose_data().
        * average_period (str, default '30min'): Averaging period for averaging the wavelet decompositioned values. Format: pandas time string, e.g. "30min". Possible specifications are s, min, h, d. Passed to the main function.
        * sitename (str, default "00000"): Sitename, files get named accordingly.
        * wt_kwargs (dict, default {}): **kwargs passed to the wavelet tranformation itself. Can e.g. include wavelet specification. See wavelet_function.py for more details. Important setting include f0, which is the lowest frequency for wavelet decomposition. Because adding buffer is necessary to prevent edge effects, a lower f0 drastically increases the amount of data loaded in and the memory usage.
        * meta (dict, default {}): Header lines in the output files. Get filled successively during the code run.
        **kwargs
        
    Return:
        fulldata (pandas.DataFrame): Containing all processed data. If integration_period is specified already integrated.
    
"""
```


or directly run from a DataFrame or dictionary:
```python
import waveletec
data = ...
waveletec.main(data, varstorun, period=None, average_period='30min', 
         cond_samp_both=False, output_kwargs={}, meta={}, **kwargs):
```
With corresponding documentation:
```python
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
```

Further important functions are:
If e.g. inside process no integration_period was specified:
```python
import waveletec

waveletec.integrate_cospectra_from_file(root, f0, pattern='_full_cospectra_([0-9]+)_', 
                                  dst_path=None, newlog=False)
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
```

If no inside process no partition was specified:
```python
import waveletec

waveletec.cs_partition_NEE_ET(site_name, output_folderpath, NEE=True, ET=True, 
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

```



#### Saving Results using the default functions

If you're using ```waveletec.process```, just pass the ```output_folderpath``` parameter like this: ```waveletec.process(..., output_folderpath='PATH/TO/OUTPUT/FOLDER/')```. This will create a new subfolder called 'wavelet_full_cospectra' for you, and you'll find the sum of all your files neatly stored in the output_folderpath you specified.

Prefer using ```waveletec.main```? No problem! Just pass the ```output_kwargs``` parameter like this: ```waveletec.main(data, ..., output_kwargs={'output_path': 'PATH/TO/OUTPUT/FOLDER'})```.


### (Recommended) Using EddyPro

1. Make sure [EddyPro®](https://www.licor.com/support/EddyPro/software.html) is installed.
2. Run EddyPro, saving level 6 raw data.
   - go to Advanced Settings (top menu) > Output Files (left menu) > Processed raw data (bottom);
   - select Time series on "level 6 (after time lag compensation)";
   - select all variables;
   - proceed as usual running on "Advanced Mode".
3. Run _waveletec_:\
    This requires ```pip install waveletec```.
    Inside a script: See example file [example.ipynb](https://github.com/pedrohenriquecoimbra/wavelet-ec/blob/main/sample/FR-Gri_20220514/example.ipynb).
    ```python
    import waveletec as wEC

    wEC.run_from_eddypro(
        path="input/EP/FR-Gri_sample.eddypro", 
        output_folderpath="output/wavelet_flux",
        datetimerange="20220513T0000-20220516T0000",
        processing_time_duration="3h",
        covariance=["w*co2|w*h2o"],
        )

    wEC.integrate_cospectra_from_file(
        "output/wavelet_flux/wavelet_full_cospectra",
        f0=1/30,
        dst_path="output/wavelet_flux/FR-Gri_sample_CDWT_full_cospectra.csv",
    )

    wEC.condition_sampling_partition(
        folder="output/wavelet_flux",
        output_name="FR-Gri_sample_CDWT_partitioning",
    )
    ```

### Using bmmflux files

1. Activate the high-frequency output inside bmmflux. These files will contain all necessary corrections.
2. Run _waveletec_:\
    Inside a script using the ```waveletec.process()``` function, with the argument ```load_kwargs = {'handle_bmmflux_raw_dataset':True}```.
    Example (example data not provided yet):
    ```python
    data = waveletec.process(datetimerange = '20250407T0805-20250410T0805',
                  fileduration = 5,
                  input_path = "../test_input_bmmflux/",
                  acquisition_frequency = 10,
                  covariance=["w*h2o|w*co2"],
                  cond_samp_both=True,
                  output_folderpath='../test_outputs/',
                  processing_time_duration="6h",
                  partition=["ET", "NEE"],
                  integration_period=30*60, 
                  average_period = "30min",
                  sitename='EC_North',
                  wt_kwargs = {'f0':(1/(1*60*60))},
                  #transform_kwargs = {'memory_eff':True}, # also default at the moment True
                  load_kwargs = {'handle_bmmflux_raw_dataset':True}
                  )
    ```
    This will by default run wavelet decomposition, conditional sampling, integrating using a high-pass filter and partition NEE and ET.


### Using the command line / terminal

After installing the package using e.g ```pip install waveletec```, 
following commands from the command line are possible:

```bash
waveletEC-eddypro_run
waveletEC-run
waveletEC-partition
waveletEC-integrate

```
using ```--help``` it is possible to see an overview of all possible options to
specifiy the commands. Those are similar to the ones when running 
manually in a script the functions
```run_from_eddypro```, ```process```, ```integrate_cospectra_from_file```, 
```cs_partition_NEE_ET```, respectively.



### Output Format cospectra

The output file for the wavelet-based (co)spectra analysis is structured as follows:

```cs
1   wavelet_based_(co)spectra
2   --------------------------------------------------------------
3   TIMESTAMP_START = 2022-05-13 00:00:00
4   TIMESTAMP_END = 2022-05-13 00:30:00
5   N: 133
6   TIME_BUFFER [min] = nan
7   frequency [Hz]
8   y-axis -> nan
9   mother_wavelet -> dwt
10  acquisition_frequency [Hz] = 20.0
11  averaging_interval [Min] = 30min
12  natural_frequency,variable,value
13  3.814697265625e-05,co2,417.0002460141172
.   ...,...,...
.   5.0,co2,1.708199383374261e-05
.   10.0,co2,1.7947312058017124e-07
.   ...,...,...
```

#### Explanation of Fields:

- **`wavelet_based_(co)spectra`**: Indicates that the data pertains to wavelet-based (co)spectra analysis.

- **`TIMESTAMP_START` and `TIMESTAMP_END`**: These fields specify the start and end times of the data collection period.

- **`N`**: Represents the number of data points or samples included in the analysis.

- **`TIME_BUFFER [min]`**: Indicates the time buffer in minutes. A value of 0 means no buffer was applied.

- **`frequency [Hz]`**: Specifies the frequency of the data is in Hertz.

- **`y-axis_->_wavelet_coefficient_*_`**: Indicates that the y-axis represents the wavelet coefficients.

- **`mother_wavelet -> dwt`**: Specifies the type of mother wavelet used in the analysis, in this case, `dwt` (Discrete Wavelet Transform).

- **`acquisition_frequency [Hz]`**: The frequency at which data was acquired, given in Hertz.

- **`averaging_interval [Min]`**: The interval over which the data was averaged, specified in minutes.

- **`natural_frequency,variable,value`**: This line contains the actual data points:
  - **`natural_frequency`**: The frequency at which the data point was recorded.
  - **`variable`**: The variable being measured, such as `co2`.
  - **`value`**: The value of the variable at the specified natural frequency.


### Output Format Integrated Spectra

The integrated spectra is the sum of all frequencies up to the frequency ```f0```
or corresponding to the specified ```integration_period``` (```f0 = 1/integration_period```).
Hence, ```f0``` works as a high-pass filter for the wavelet cospectrum.

In general, in the cospectra and the integrated spectra, the variables containing
```_qc``` specify quality control. In case of discrete wavelet transform (dwt) these
show the amount of NaN-values (in the cospectra the ratio, in the integrated spectra summed up over the frequencies).
In case of continous wavelet transform, they also consider the amount of values outside of the cone of influence.

Values such as ```wco2+wh2o+``` are conditionally sampled fluxes. 
In this example it is the ```wco2``` flux, sampled when the ```wco2``` gust flux is positive AND
the ```wh2o``` gust flux is positive (i.e. on the high frequency data).
Hence, the order of the fluxes ```wco2+wh2o+``` vs ```wh2o+wco2+``` is important!
The second case is the ```wh2o``` flux, when the ```wh2o``` gust flux is positive AND
the ```wco2``` gust flux is positive.

```cs
TIMESTAMP,co2,co2_qc,h2o,h2o_qc,w,w_qc,wco2,wco2+wh2o+,wco2+wh2o-,wco2-wh2o+,wco2-wh2o-,wh2o
2025-04-07 08:00:00,-0.10021955266763895,6.998833138856476,22.33837837212945,6.998833138856476,-0.03433317151915598,6.991831971995333,0.020236662912589777,0.009884740387139105,0.017499776507520038,-0.00022041844773451258,-0.006927435534334856,-1.6948474791162642
2025-04-07 12:00:00,0.2965510005139298,0.21466666666666664,51.18079698756063,0.21466666666666664,-0.06683468485321119,2.291333333333333,0.0681909688582618,0.02922524693866301,0.044957744497812195,
```

### Output Format Partitioned

```cs
TIMESTAMP,ET,T,E,Dew,NEE,GPP,Reco
2025-04-07 08:00:00,-1.6948474791162642,0.000914847844653,0.0583496289049483,-1.7541119558658658,0.0202366629125897,-0.0036841362149019003,0.0239207991274916
```

## Example

For an example follow the [demonstration.ipynb](https://github.com/pedrohenriquecoimbra/wavelet-ec/blob/main/sample/FR-Gri_20220514/demonstration.ipynb) and [example.ipynb](https://github.com/pedrohenriquecoimbra/wavelet-ec/blob/main/sample/FR-Gri_20220514/example.ipynb) files.

## License

This project is licensed under the EU EUPL License. See the [LICENSE](LICENSE) file for details.
