reports/02_model_metabolites_censored_amphetamine.RData: reports/02_model_metabolites_censored.Rmd data/water_cleaned.txt
	cd reports/;\
	Rscript -e "rmarkdown::render('02_model_metabolites_censored.Rmd')"