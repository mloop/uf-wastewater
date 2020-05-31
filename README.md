# uf-wastewater

The purpose of this project is to estimate the concentrations of illicit compounds near Ben Hill Griffith Statium on a gameday.

## Summary of scripts

* `scripts/`
    * `01.eda`: creates some of the figures for the paper
    * `02_*`: fits the main models for the paper and performs prior and posterior checks on model fits; create summaries of the model results for the paper
    * `03_enter_metabolism_data.Rmd`: creates a dataset with information on metabolites provided by Bikram Subedi
    * can ignore other scripts in this directory for the purposes of the main paper
    * `clean_dataset.R`: takes the raw data and converts it to an analysis-ready dataset
