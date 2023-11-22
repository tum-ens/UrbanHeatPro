# UrbanHeatPro

[//]: # (marker_text_start_description)

A Bottom-up model for the simulation of heat demand profiles of urban areas

----------------------------------------------------------------------------

## Features

  - UrbanHeatPro is a Python-based bottom-up model for the simulation of heat demand profiles of urban areas.
  - It considers both the space heating demand and hot water demand. So far, the hot water demand is calculated only for residential buildings.
  - Characteristic values for the building stock, building thermal properties, building set-temperature and annual hot water consumption are based on statistics for Germany.
  - DSM strategies for the reduction of heat demand such as building renovation, heat load reduction, night-set back operation, etc, are easily implemented.
  - The size of the study area can start from one building. Buildings to simulate should be included in the input (csv) file.
  - By default, the model operates on an hourly time steps. However, the temporal resolution is configurable.

[//]: # (marker_text_end_description)

## Requirements & Installation

### System requirements

The project uses POSIX paths as is standard in UNIX systems such as Linux or macOS. On Windows machines using the WSL is might mitigate possible issues.


### Python requirements & dependencies
Python 3.9 or higher (installation with Anaconda recommended).

The latest version was tested with Python 3.10.

Use the provided requirements.txt or the setup.py to install the dependencies via pip.

### Installation

There are two ways to install UrbanHeatPro:

#### 1. As a package:

Ideally install the project as an editable package.
Inside the root directory of the repository for UrbanHeatPro run:
    
```bash
$ pip install -e .
```

All dependencies are installed and UrbanHeatPro can be used as the package `UrbanHeatPro`.
This allows to use the folder structure expected by UrbanHeatPro and use it from other packages.

#### 2. As a standalone project with an entry script

To install only the dependencies via pip, run inside the root directory of the repository for UrbanHeatPro:

```bash
$ pip install -r requirements.txt
```


## Folder structure & settings

[//]: # (marker_text_start_folders_files)

The expects the input files in the `input/ directory`.
Per default the results are written to the `results/` directory. A different result directory can be defined in the 
settings file or as a parameter of `run_uhp()`.

The example and default settings files are in the `settings/` directory.
Settings are provided to the module as a configuration file in the yaml format. An example of the expected structure 
can be seen in the example settings file `settings/uhp_settings_example.yaml`. All possible settings are also shortly 
described in this file.

**Note:** Do not move or modify the default settings file `settings/uhp_default_settings.yaml`. This file is needed in case no other settings are given.

## Input files (csv)
An example for the expected input data is given for the region Unterhaching. 

### Buildings
Each building is described with the following information:
1.  Area (required) \
  Building ground floor area in m²
2.  Use (required) \
    Building use as integer:
    *   0 Commercial
    *   1 Industrial
    *   2 Public
    *   3 Residential

3.  Bid (optional) \
  Building identification number as integer
4.  Free_walls (required) \
  Number of walls in contact with ambient temperature. The building is assumed to be a rectangular box with four walls.
5.  Construction year class (optional) \
  Construction year class from TABULA Typology as integer:

    *   0 <1859
    *   1 1860 - 1918
    *   2 1919 - 1948
    *   3 1949 - 1957
    *   4 1958 - 1968
    *   5 1969 - 1978
    *   6 1979 - 1983
    *   7 1984 - 1994
    *   8 1995 - 2001
    *   9 2002 - 2009
    *   10 >2009

6.  Building type (optional) \
  Building type from TABULA Typology as integer:

    *   0 Single-Family House (SFH)
    *   1 Terraced House (TH)
    *   2 Multi-family House (MFH)
    *   3 Apartment Block (AB)

7.  Refurbishment level (optional) \
  Refurfishment level for all building elements from TABULA Typology as integer:

    *   1 No refurbishment (Existing state)
    *   2 Usual refurbishment
    *   3 Advanced refurbishment

8.  Number of occupants (optional) \
  Number of occupants living in the building

All columns are required, if no information is given, the cell value is taken as `NaN`.
This input file should be located in the corresponding folder `input/Buildings`.

### Regional Data
In addition to the buildings to be modeled, regional data is needed. This data should be placed in the folder `input/Regional Data`.
Each region has its own directory of regional data, e.g. `input/Regional Data/Unterhaching`.

[//]: # (marker_text_end_folders_files)

## Using UrbanHeatPro
### From the entry script
To run the model with the given input data, change the desired information in the runme.py file, create a settings 
file as in the given example and run this file in the command line:

```sh
$ python runme.py
```

### As a library or package
Alternatively the model can be used as a library or package. 
Note: This requires UrbanHeatPro to be installed as a package (see Installation 1).
The `run_uhp` function in the `run_uhp.py` module is the entrypoint for running the model.


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.


## Copyright
Copyright (C) 2018-2021 Anahi Molar-Cruz

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details: http://www.gnu.org/licenses/.
