#The following packages must be installed:
#xcms
#magrittr
#highcharter
#pracma
#plyr
#lsa
#reshape2
#htmlwidgets

setwd('path/to/analyse_GCMS_data/')
source('./analysis_functions.R')

data <- easy_start(data_dir = './sample_data/',
                   my_prefix = './easy_start_sample_output/sample')
