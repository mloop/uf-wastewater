reports/02_estimate_dosage_censored.png: reports/02_estimate_dosage_censored.Rmd reports/02_model_metabolites_censored_amphetamine.RData
	cd reports/;\
	Rscript -e "rmarkdown::render('02_estimate_dosage_censored.Rmd')"


reports/02_model_metabolites_censored_amphetamine.RData: reports/02_model_metabolites_censored.Rmd data/water_cleaned.txt
	cd reports/;\
	Rscript -e "rmarkdown::render('02_model_metabolites_censored.Rmd')"
	
reports/04_model_metabolites_censored_multivariate.RData: reports/04_model_metabolites_censored_multivariate.R
	cd reports/;\
	Rscript -e 04_model_metabolites_censored_multivariate.R