# `raw-data/`

* Catalogue and documentation of raw data compiled for biological indicators for steelhead starting in 2022 

* Includes data files in original form (as-received) and source PDFs from which raw data have been extracted

* Intended as read-only

# Data types

* spawner surveys counts and abundance indexes/estimates
* juvenile surveys
* hatchery releases
* test fishery data
* age and life history data
* catch data 

# Important files

`1-raw-data-catalogue` is an inventory of raw data files with region/stream, source, and data type for each data file 

# Processing Notes: Data Extraction From PDFs

Process for extracting data from PDF reports:

* Move PDF to `raw-population-indicator-data` 

* Extract data from tables in PDF (Adobe Acrobat Pro: File > Export to... > Spreadsheet > Excel) into spreadsheet and save as CSV in `raw-population-indicator-data`. If this function does not work for a given PDF, data are extracted by hand or from figures using [Plot Digitizer](https://plotdigitizer.com/app).

* QAQC formatting and data checks (check row and column names and alignment; spot check data values)

* Add to as new row to `1-raw-data-catalogue`. 


