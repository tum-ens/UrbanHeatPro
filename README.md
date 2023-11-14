# UrbanHeatPro
A Bottom-up model for the simulation of heat demand profiles of urban areas

----------------------------------------------------------------------------

## Features

  - UrbanHeatPro is a python-based bottom-up model for the simulation of heat demand profiles of urban areas.
  - It considers both the space heating demand and hot water demand. So far, the hot water demand is calculated only for residential buildings.
  - Characteristic values for the building stock, building thermal properties, building set-temperature and annual hot water consumption are based on statistics for Germany.
  - DSM strategies for the reduction of heat demand such as building renovation, heat load reduction, night-set back operation, etc, are easily implemented.
  - The size of the study area can start from one building. Buildings to simulate should be included in the input (csv) file.
  - By default, the model operates on an hourly time steps. However, the temporal resolution is configurable.


## Requirements & Installation

### System requirements

The project uses POSIX paths as is standard in UNIX systems such as Linux or MacOS. On Windows machines using the WSL is might mitigate possible issues.


### Python requirements & dependencies
Python 3.9 or higher (installation with Anaconda recommended).

The latest version was tested with Python 3.10.

Use the provided requirements.txt or the setup.py to install the dependencies via pip.

### Installation

Ideally install the project as editable.
Inside the directory of the repository for UrbanHeatPro run:
> pip install -e .

This allows to use the folder structure expected by UrbanHeatPro.

## Folder structure & settings
The expects the input files in the input/ directory.
Per default the results are written to the results/ directory. A different result directory can be defined in the 
settings file or as a parameter of `run_uhp()`.

The example and default settings files are in the settings/ directory.
Settings are provided to the module as a configuration file in the yaml format. An example of the expected structure 
can be seen in the example settings file settings/uhp_settings_example.yaml. All possible settings are also shortly 
described in this file.

Note: The default settings file should not be moved or modified as it is needed in case no other settings are given.

## Input files (csv)
An example for the expected input data is given for the region Unterhaching. 

### Buildings
Each building is described with the following information:
<ol>
  <li>Area (required)</li>
      Building ground floor area in m²
  <li>Use (required)</li>
      Building use as integer:
        <ul>
          <li>0	Commercial</li>
          <li>1	Industrial</li>
          <li>2	Public</li>
          <li>3	Residential</li>
        </ul>
  <li>Bid (optional)</li>
      Building identification number as integer
  <li>Free_walls (required)</li>
      Number of walls in contact with ambient temperature. The building is assumed to be a rectangular box with four walls.
  <li>Construction year class (optional)</li>
      Construction year class from TABULA Typology as integer:
        <ul>
          <li>0	<1859</li>
          <li>1	1860 - 1918</li>
          <li>2	1919 - 1948</li>
          <li>3	1949 - 1957</li>
          <li>4	1958 - 1968</li>
          <li>5	1969 - 1978</li>
          <li>6	1979 - 1983</li>
          <li>7	1984 - 1994</li>
          <li>8	1995 - 2001</li>
          <li>9	2002 - 2009</li>
          <li>10 >2009</li>
        </ul>
  <li>Building type (optional)</li>
      Building type from TABULA Typology as integer:
        0	Single-Family House (SFH)
        1	Terraced House (TH)
        2	Multi-family House (MFH)
        3	Apartment Block (AB)
  <li>Refurbishment level (optional)</li>
      Refurfishment level for all building elements from TABULA Typology as integer:
        1	No refurbishment (Existing state)
        2	Usual refurbishment
        3	Advanced refurbishment
  <li>Number of occupants (optional)</li>
      Number of occupants living in the building
        
</ol>

All columns are required, if no information is given, the cell value is taken as NaN.
This input file should be located in the corresponding folder `input/Buildings`.

### Regional Data
In addition to the buildings to be modeled, regional data is needed. This data should be placed in the folder `input/Regional Data`.
Each region has its own directory of regional data, e.g. `input/Regional Data/Unterhaching`.

## Using UrbanHeatPro
### From the entry script
To run the model with the given input data, change the desired information in the runme.py file, create a settings 
file as in the given example and run this file in the command line:
> python runme.py

### As a library
Alternatively the model can be used as a library.
The `run_uhp` function in the `run_uhp.py` module is the entrypoint for running the model.


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.


## Copyright
Copyright (C) 2018-2021 Anahi Molar-Cruz

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details: http://www.gnu.org/licenses/.
