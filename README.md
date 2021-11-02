# A space-time-age model for subnational child mortality estimation in India

This repository contains code used to produce the results described in Chapter 4 of the doctoral thesis, "Assessing local health outcomes using spatially-resolved health surveillance data."

## Code structure

This repository is organized into the following folders:
- `inputs/`: Summarized input data from the Sample Registration System of India (SRS)
- `data_prep/`: Scripts used to prepare raw data for use in the geostatistical model and secondary analyses
- `model/`: Functions and scripts used to run the space-time-age model based on survey data. The main execution script, `run_u5m_model.R`, calls functions defined in a [public repository](https://github.com/ihmeuw/lbd/tree/u5m-lmic-2019) based on settings defined in `config.csv`, to estimate neonatal, infant, and under-5 mortality between 2000 and 2017. The script `project_d_and_q.R` projects estimates forward to 2025 and 2030 to compare with goals set in the 2017 National Health Strategy.
- `viz/`: Scripts to analyze the raw data and model results, including all figures included in the chapters.


## Running this code

This code is designed to be executed in R version 3.6 or above. Most model configuration arguments should be set in `config.csv`; otherwise, filepaths to the input data directory and this repository should be set at the top of each script.
