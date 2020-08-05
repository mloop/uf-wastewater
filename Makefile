output/02_posterior_predictive_mass_load_consumption_locations.rds: R/02_posterior_predictive_distributions_doses.R data/stadium_seating.txt data/flow_rate.txt output/02_model_metabolites_censored_Amphetamine.rds data/water_cleaned.txt data/03_metabolism_data.txt
	cd R/ && Rscript 02_posterior_predictive_distributions_doses.R
	
output/02_model_metabolites_censored_Amphetamine.rds: R/02_model_metabolites_censored.R
	cd R/ && Rscript 02_model_metabolites_censored.R
	
