
# download datasets
# need to add this script
# http://www.esapubs.org/archive/ecol/E094/195/SCBI_initial_woody_stem_census_2012.csv
# 


# filter raw data
Rscript('./scripts/empir_data_filering')

# run analysis
system('Rscript ./scripts/mobr_analysis.R', wait=F)
