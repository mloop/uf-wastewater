reports/02_estimate_dosage_censored_whole_stadium_whole_game.pdf: reports/02_estimate_dosage_censored_whole_stadium_whole_game.Rmd output/02_posterior_predictive_mass_load.rds
	cd reports/ && Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/MacOS/pandoc');library(rmarkdown);render('02_estimate_dosage_censored_whole_stadium_whole_game.Rmd')"

reports/02_estimate_dosage_censored_time_location.pdf: reports/02_estimate_dosage_censored_time_location.Rmd output/02_posterior_predictive_mass_load.rds
	cd reports/ && Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/MacOS/pandoc');library(rmarkdown);render('02_estimate_dosage_censored_time_location.Rmd')"

figs/02_posterior_predictive_doses_fig.png: R/02_posterior_predictive_doses_fig.R output/02_posterior_predictive_mass_load.rds
	cd R/ && Rscript 02_posterior_predictive_doses_fig.R

figs/02_posterior_predictive_mass_load_time.png: R/02_posterior_predictive_mass_load_time_fig.R output/02_posterior_predictive_mass_load.rds
	cd R/ && Rscript 02_posterior_predictive_mass_load_time_fig.R

figs/02_posterior_predictive_mass_load_time_location.png: R/02_posterior_predictive_mass_load_time_location_fig.R output/02_posterior_predictive_mass_load.rds
	cd R/ && Rscript 02_posterior_predictive_mass_load_time_location_fig.R

output/02_posterior_predictive_mass_load.rds: R/02_posterior_predictive_distributions_doses.R data/stadium_seating.txt data/flow_rate.txt output/02_model_metabolites_censored_Amphetamine.rds data/water_cleaned.txt data/03_metabolism_data.txt
	cd R/ && Rscript 02_posterior_predictive_distributions_doses.R

output/02_model_metabolites_censored_Amphetamine.rds: R/02_model_metabolites_censored.R
	mkdir -p output/
	cd R/ && Rscript 02_model_metabolites_censored.R
