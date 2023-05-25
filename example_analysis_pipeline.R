###Set working directory
setwd('path/to/analyse_GCMS_data/')
getwd()

###Load analysis functions
source('./analysis_functions.R')


###Get library
#baseline of chromatogram, calculated from the mean of several conditioned tubes,
#then passed through a low-pass filter
#see 2022-01-04_mzR_get_background.R for details
load('smoothed.RData')

#custom library with retention times
#see 2022-11-15_build_compound_library.R for details
load('2022-11-17_custom_lib4.Rdata')
#Will add functionality later to allow adding new compounds



###Load and process data

data_dir <- './sample_data'

library(xcms) #3.6.2
library(magrittr)

data1 <- runXCMS(data_dir, rt_span = c(500, 1800), reference_sample = 'B37') %>%
  clean_sample_names() %>%
  subtract_blank(blank_IDs_file = './sample_data/blank_IDs.tab') %>% #table shows which blank corresponds to which sample
  round_mz_rt()

saveRDS(data1, './example_analysis_output/sample_raw_XCMS_data.rds') #save this intermediate step

###Some minor automatic RT adjustment
custom_lib4 <- finetune_RTs(custom_lib = custom_lib4,
                            custom_lib_spectra = custom_lib_spectra4,
                            data = data1)

#component distinctiveness, for estimating compound abundances
distinctiveness_list <- get_distinctiveness_list(custom_lib4, custom_lib_spectra4, smoothed)

###Round 1 of abundance estimation
round1_output <- estimate_abund(data1, custom_lib = custom_lib4,
                 custom_lib_spectra = custom_lib_spectra4,
                 d_list = distinctiveness_list) %>%
  filter_hits_and_subtract(data = data1,
                           custom_lib = custom_lib4,
                           custom_lib_spectra = custom_lib_spectra4)

data2 <- round1_output$chromatograms

###Round 2 of abundance estimation
round2_output <- estimate_abund(data2, custom_lib = custom_lib4,
                                custom_lib_spectra = custom_lib_spectra4,
                                d_list = distinctiveness_list) %>%
  filter_hits_and_subtract(data = data2,
                           custom_lib = custom_lib4,
                           custom_lib_spectra = custom_lib_spectra4)

###Clean-up and processing
cpd_sample_matrix <- combine_estimates(round1_output$abund_estimates, round2_output$abund_estimates) %>%
  long_to_wide() %>%
  remove_nonbio() %>%
  relative_abund()

colour_df <- get_compound_colours(cpd_sample_matrix,
                                  compound_info = './2022-12-10_custom_lib4_cpd_classes.tab',
                                  special_cpds = './special_cpd_colours.tab',
                                  simple_names = TRUE)

cpd_sample_matrix_reordered <- reorder_matrix(cpd_sample_matrix,
                                              axis = 'col',
                                              order = './sample_data/sample_order.txt') %>%
  reorder_matrix(axis = 'row',
                 order = colour_df$cpd_name)

###Plotting
static_barplot(cpd_sample_matrix_reordered,
               cols = colour_df$colours,
               legend_df = readRDS('legend_df.rds'))

pdf_barplot(cpd_sample_matrix_reordered,
            cols = colour_df$colours,
            legend_df = readRDS('legend_df.rds'),
            filename = './example_analysis_output/sample')

library(highcharter)
interactive_barplot(cpd_sample_matrix_reordered,
                    colour_df = colour_df,
                    filename = './example_analysis_output/sample',
                    use_simple_names = TRUE)

###Output data matrix
write.table(cpd_sample_matrix, file = './example_analysis_output/sample_output_matrix.tab',
            sep = '\t',
            row.names = TRUE,
            col.names=NA)
