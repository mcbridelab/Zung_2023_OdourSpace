distinctive_peaks <- function(spectrum1, spectrum2, rt1, rt2){
  #find the "distinctiveness" of ions in a given spectrum, weighted by how closely in time compounds elute
  peakwidth <- 10 #bigger value means ions will not be counted as distinctive, even for compounds that elute far away
  weight <- exp(-(1/peakwidth)*(rt1-rt2)^2)
  distinctiveness <- spectrum1 - weight*spectrum2
  distinctiveness[distinctiveness < 0] <- 0
  return(distinctiveness)
}

get_distinctiveness_list <- function(custom_lib, custom_lib_spectra, smoothed){
  median_spec <- apply(custom_lib_spectra, 1, median)
  #get "distinctiveness" of ions
  distinctiveness_list <- list()
  
  for(i in 1:dim(custom_lib)[1]){
    focal_cpd <- custom_lib$cpd_name[i]
    
    distinctiveness <- custom_lib_spectra[,i]
    
    for(j in 1:dim(custom_lib)[1]){
      if(i == j){
        next
      }
      #Start with full focal spectrum
      #Loop through all other compounds, taking away from the distinctiveness
      #of each ion in the focal spectrum according to the function above
      distinctiveness <- distinctive_peaks(distinctiveness, custom_lib_spectra[,j],
                                           custom_lib$rt[i], custom_lib$rt[j])
      
    }
    
    #Subtract the median compound spectrum (downweighted slightly by adding an artificial RT difference of 0.5s)
    distinctiveness <- distinctive_peaks(distinctiveness, median_spec,
                                         custom_lib$rt[i], 0.5+custom_lib$rt[i])
    
    #Subtract the background
    distinctiveness <- distinctiveness / smoothed[,round(custom_lib$rt[i])]
    distinctiveness[distinctiveness < 0] <- 0
    distinctiveness_list[[i]] <- distinctiveness
  }
  
  return(distinctiveness_list)
}


runXCMS <- function(data_dir,
                    rt_span = NULL,
                    reference_sample = NULL){
  #Takes in a directory of mzdata.xml files
  #Outputs a list of features (peak of a certain m/z at a certain RT) and their abundance across samples
  
  #rt_span example value: c(500, 1800). Ignore components outside this retention-time window
  
  #reference_sample example value: 'B37'. An identifying portion of the file name to use as the reference
  #sample. In the RT-alignment step, RTs will be aligned to this value
  
  setMSnbaseVerbose(opt = TRUE)
  
  files <- list.files(data_dir, full.name=TRUE, pattern = 'mzdata.xml$')
  
  raw_data <- readMSData((files = files),
                         mode = 'onDisk')
  
  if(is.null(rt_span)){
    rt_span = c(min(rtime(raw_data)), max(rtime(raw_data)))
  }
  
  raw_data <- filterRt(raw_data, rt_span)
  
  cw <- CentWaveParam(ppm = 3000,
                      prefilter = c(5,1000),
                      peakwidth = c(2,8),
                      snthresh = 4)
  
  register(SerialParam()) #next step fails for large datasets when run with parallel processing; switch to serial
  xdata <- findChromPeaks(raw_data, param = cw) #takes ~60 min on standard desktop for dataset of ~250 samples
  
  if(!is.null(reference_sample)){
    ref <- which(grepl(reference_sample, files, fixed = TRUE)) #align RTs to this sample
  } else {
    ref <- 1
  }
  
  xdata2 <- adjustRtime(xdata,
                        param = ObiwarpParam(binSize = 0.6,
                                             centerSample = ref)) 
  
  
  pdp <- PeakDensityParam(sampleGroups = rep('group1', length(files)),
                          minFraction = 0, bw = 0.3)
  
  xdata3 <- groupChromPeaks(xdata2, param = pdp)
  
  xdata4 <- fillChromPeaks(xdata3, param = FillChromPeaksParam()) 
  
  results <- featureValues(xdata4)
  results[is.na(results)] <- 0
  results <- cbind(data.frame(featureDefinitions(xdata4)@listData[c('mzmed','rtmed')]), results)
  
  
  return(results)
}


clean_sample_names <- function(data){
  #remove '.mzdata.xml' from sample names
  colnames(data) <- c('mzmed', 'rtmed',
                      substr(as.character(colnames(data)[3:dim(data)[2]]),1,
                             nchar(as.character(colnames(data)[3:dim(data)[2]]))-11))
  return(data)
}

subtract_blank <- function(data, blank_IDs_file){
  #Input: XCMS data and tab-delimited file noting which blank corresponds to which sample
  blank_IDs <- read.table(blank_IDs_file, sep='\t', header=T)
  blank_IDs$sample <- as.character(blank_IDs$sample)
  blank_IDs$blank <- as.character(blank_IDs$blank)
  
  data3 <- data
  for(i in 1:dim(blank_IDs)[1]){
    sample_ix <- which(colnames(data)==blank_IDs$sample[i])
    blank_ix <- which(colnames(data)==blank_IDs$blank[i])
    data3[,sample_ix] <- data3[,sample_ix] - data[,blank_ix]
  }
  data3[data3<0] <- 0
  data3 <- with(data3, data3[order(rtmed),])
  
  data3 <- data3[,!colnames(data3) %in% blank_IDs$blank]
  
  return(data3)
}

round_mz_rt <- function(data8){
  #round mz to nearest integer
  #round RT to nearest 0.5
  data8$mzmed <- round(data8$mzmed)
  
  data8$rtmed <- round(data8$rtmed*2)/2
  
  data8 <- aggregate(list(data8[,3:dim(data8)[2]]), by=list(mzmed=data8$mzmed, rtmed=data8$rtmed), sum)
  
  return(data8)
}


finetune_RTs <- function(custom_lib, custom_lib_spectra, data,
                         threshold = 0.7,
                         rt_wiggle_range = seq(-5,5, 0.5), #range of 'wiggle' values to try, in seconds
                         print_report = FALSE){
  #Retention times may not be exactly right; drift may occur over long periods of time
  #This function 'wiggles' the library RT slightly to see if the match score is
  #improved. New wiggled RT is accepted if match score is above threshold.
  #The function is only designed to fix RTs that are uniformly different across
  #samples. If sample RTs differ (e.g., because they were run far apart in time),
  #need to fix a different way. May be possible using XCMS.
  
  print('Adjusting RTs...')
  
  best_rt_wiggle <- data.frame(cpd = custom_lib$cpd_name,
                               match_score = NA,
                               rt_wiggle = NA,
                               improvement = NA)
  
  for(i in 1:dim(custom_lib)[1]){
    cpd1 <- custom_lib$cpd_name[i]
    #print(cpd1)
    spectrum <- custom_lib_spectra[,i]
    
    match_score_df <- data.frame(rt_wiggle = rt_wiggle_range,
                                 match_score = NA)
    
    for(j in 1:nrow(match_score_df)){
      rt_wiggle <- match_score_df$rt_wiggle[j]
      rt <- custom_lib$rt[i] + rt_wiggle
      
      ions_in_range <- subset(data, rtmed <= rt + 2 & rtmed >= rt - 2)
      if(nrow(ions_in_range)==0){
        match_score_df$match_score[j] <- 0
      } else {
        ions_in_range$mzmed <- round(ions_in_range$mzmed)
        
        #just sum across all samples
        foo <- cbind(ions_in_range$mzmed, rowSums(ions_in_range[,3:ncol(ions_in_range)]))
        colnames(foo) <- c('mzmed', 'count')
        
        foo <- aggregate(count ~ mzmed, foo, sum) #aggregated ions_in_range across duplicate mz
        filled_in <- data.frame(mz=c(40:250), count = rep(0, 211))
        filled_in[filled_in$mz %in% ions_in_range$mzmed,2] <- foo[foo$mzmed %in% filled_in$mz,2]
        
        match_score <- lsa::cosine(spectrum, filled_in$count)
        match_score_df$match_score[j] <- match_score
      }
    }
    
    ix <- which.max(match_score_df$match_score)
    best_rt_wiggle$match_score[i] <- match_score_df$match_score[ix] #the best match score achieved
    best_rt_wiggle$rt_wiggle[i] <- with(match_score_df[match_score_df$match_score > threshold,],
                                        weighted.mean(rt_wiggle, match_score)) #rt_wiggle decided on
    best_rt_wiggle$improvement[i] <- match_score_df$match_score[ix] - match_score_df$match_score[match_score_df$rt_wiggle==0]
    
  }
  
  
  
  #only accept those wiggles with final match score above threshold
  impose_wiggles <- best_rt_wiggle$rt_wiggle
  impose_wiggles[is.nan(impose_wiggles)] <- 0
  custom_lib$rt <- custom_lib$rt + impose_wiggles
  
  
  
  
  if(print_report){
    best_rt_wiggle$match_score <- round(best_rt_wiggle$match_score*100)/100
    best_rt_wiggle$improvement <- round(best_rt_wiggle$improvement*100)/100
    print(best_rt_wiggle)
  }
  
  return(custom_lib)
  
}



estimate_abund <- function(data, custom_lib, custom_lib_spectra, d_list){
  #estimate abundance of compounds based on component abundances
  #NO filter for overall spectrum match
  
  all_wtd_means = data.frame(sample=as.character(),
                             estimate=as.numeric(),
                             truth=as.numeric(),
                             cpd=as.character())
  
  nsamp = dim(data)[2]-2
  
  for(i in 1:dim(custom_lib)[1]){
    cpd <- custom_lib$cpd_name[i]
    print(cpd)
    spectrum <- custom_lib_spectra[,i]
    rt <- custom_lib$rt[i]
    
    ions_in_range <- subset(data, rtmed <= rt+2 & rtmed >= rt-2) #use data from build_compound_library.R
    ions_in_range$mzmed <- round(ions_in_range$mzmed)
    
    if(!pracma::isempty(intersect(ions_in_range$mzmed, c(40:250)[spectrum!=0]))){
      #if there are distinctive spectrum ions among detected features in range
      estimates <- data.frame(matrix(ncol=4, nrow=0))
      colnames(estimates) <- c('sample','mz','estimate','distinctiveness')
      
      for(j in 1:dim(ions_in_range)[1]){
        mz <- ions_in_range$mzmed[j] #focal mz
        mz_ix <- mz-39 #convert from mz to index of mz in sequence c(40:250)
        if(any(c(40:250)[spectrum!=0]==mz)){ #if focal mz is found in focal compound spectrum
          factor <- sum(spectrum)/spectrum[mz_ix] #scaling factor - what prop of ion counts in spectrum at this mz?
          if(is.finite(factor)){
            est <- t(ions_in_range[j,3:(nsamp+2)]) * factor
            rownames(est) <- NULL
            colnames(est) <- 'estimate'
            dist <- d_list[[i]][mz_ix] #distinctiveness of ion for cpd
            temp <- data.frame(sample=colnames(data)[3:(nsamp+2)], mz = rep(mz, nsamp), estimate = est, distinctiveness = rep(dist, nsamp))
            estimates <- rbind(estimates, temp)
          }
        }
      }
      
      wt_mean <- plyr::ddply(estimates, ~sample, plyr::summarize, weighted.mean(estimate, w=distinctiveness))
      
      wt_mean$cpd <- cpd
      colnames(wt_mean) <- c('sample','estimate','cpd')
      all_wtd_means <- rbind(all_wtd_means, wt_mean)
      
      
    }
  }
  
  
  return(all_wtd_means)
  
}



filter_hits_and_subtract <- function(all_wtd_means, data, custom_lib, custom_lib_spectra){
  #filter out hits that have <0.7 cosine similarity (set abundance to zero)
  #also subtract component abundances for hits that pass threshold
  #return (1) the estimates with filtered hits and (2) the 'residual sample'
  #left after subtracting (portions of) components that are already
  #accounted for by compounds that met threshold
  
  all_wtd_means2 <- all_wtd_means
  nsamp = dim(data)[2]-2
  data4 <- data
  
  for(i in 1:dim(custom_lib)[1]){
    cpd1 <- custom_lib$cpd_name[i]
    print(cpd1)
    spectrum <- custom_lib_spectra[,i]
    rt <- custom_lib$rt[i]
    
    ions_in_range <- subset(data, rtmed <= rt + 2 & rtmed >= rt - 2)
    ions_in_range$mzmed <- round(ions_in_range$mzmed)
    
    sample_estimates <- subset(all_wtd_means, cpd == cpd1)
    sample_estimates$factors1 <- sample_estimates$estimate / sum(spectrum) #est abundance of compounds, normalized to standard spectrum
    
    for(j in 1:nsamp){ #sample
      samp <- colnames(data)[j+2]
      foo <- ions_in_range[,c(1,j+2)]
      if(sum(foo[,2]) > 0){
        colnames(foo) <- c('mzmed', 'count')
        foo <- aggregate(count ~ mzmed, foo, sum) #aggregated ions_in_range across duplicate mz
        filled_in <- data.frame(mz=c(40:250), count = rep(0, 211))
        filled_in[filled_in$mz %in% ions_in_range$mzmed,2] <- foo[foo$mzmed %in% filled_in$mz,2]
        scaled_spectrum <- data.frame(mz=c(40:250), count=spectrum * sample_estimates$factors1[j])
        
        match_score <- lsa::cosine(spectrum, filled_in$count)
        
        if(match_score < 0.7){
          #if no match, set estimate to zero (compound not actually present)
          all_wtd_means2[all_wtd_means2$sample==samp & all_wtd_means2$cpd == cpd1,2] <- 0
          
        } else{
          #subtract scaled_spectrum from the chromatogram that will go through another iteration of analysis
          scaled_spect_no_zeroes <- scaled_spectrum[scaled_spectrum$count!=0,]
          for(k in 1:dim(scaled_spect_no_zeroes)[1]){
            mz1 <- scaled_spect_no_zeroes$mz[k]
            ion_index <- which.min((sqrt((5*(mz1 - data$mzmed))^2 + (rt - data$rtmed)^2)))
            data4[ion_index, which(colnames(data4)==samp)] <- data4[ion_index, which(colnames(data4)==samp)] - scaled_spect_no_zeroes$count[k]
            
          }
        }
      }
    }
  }
  data4[data4<0] <- 0
  
  return(list(chromatograms = data, abund_estimates = all_wtd_means2))
}


combine_estimates <- function(round1_estimates, round2_estimates){
  return(with(rbind(round1_estimates, round2_estimates),
              aggregate(list(estimate=estimate), by=list(sample=sample,cpd=cpd), FUN=sum)))
}

long_to_wide <- function(long){
  # check for duplicates, correct as needed
  x <- paste(long$cpd,long$sample)
  if(any(duplicated(x))){
    print(x[duplicated(x)])
    stop('Duplicates found (i.e., multiple estimates of the same compound for the same sample).')
  }
  #convert to wide format
  wide <- reshape2::dcast(long, cpd ~ sample, value.var='estimate')
  wide[is.na(wide)]<-0
  rownames(wide) <- wide$cpd
  wide <- wide[, -1]
}

remove_nonbio <- function(wide){
  #remove siloxanes and halogenated cpds
  wide <- wide[!grepl('sil', row.names(wide)),]
  wide <- wide[!grepl('TMS', row.names(wide)),]
  wide <- wide[!grepl('bromo', row.names(wide)),]
  wide <- wide[!grepl('chlor', row.names(wide)),]
  wide <- wide[!grepl('fluor', row.names(wide)),]
  wide <- wide[!grepl('Bromo', row.names(wide)),]
  wide <- wide[!grepl('Chlor', row.names(wide)),]
  wide <- wide[!grepl('Fluor', row.names(wide)),]
  
  return(wide)
}

relative_abund <- function(wide){
  return(apply(as.matrix(wide),2,function(x){x/sum(x)}))
}

reorder_matrix <- function(prop, axis, order = NULL){
  #reorder samples in the compounds x samples matrix
  #can pass a numeric vector (e.g., order = c(1,5,3,4,2))
  #or the path to a file with the names in the desired order
  #or a character vector with the names in the desired order
  if(!is.null(order)){
    if(length(order)==1){ #if length is one, assume this is a file path
      order <- as.character(read.table(order)$V1)
    }
    
    if((axis == 'row' & length(order) != nrow(prop))|
       (axis == 'col' & length(order) != ncol(prop))){
      warning('Mismatch between matrix size and length of order vector.')
    }
    
    if(is.character(order)){
      if(axis == 'row'){
        order <- order[order %in% rownames(prop)]
      } else if(axis == 'col'){
        order <- order[order %in% colnames(prop)]
      }
      
    }
    
    if(axis == 'row'){
      prop <- prop[order,]
    } else if(axis == 'col'){
      prop <- prop[,order]
    } else {
      stop("Axis must be 'row' or 'col'.")
    }
  }
  return(prop)
}

get_compound_colours <- function(prop, compound_info, special_cpds = NULL, simple_names = FALSE){
  #prop: the matrix of compounds x samples
  #compound_info: the path to a file containing info about what class compounds belong to
  #special_cpds: compounds with their own specially assigned colour
  
  cpds_class <- read.table(compound_info, sep='\t', header=T)
  cpds_class$cpd_name <- as.character(cpds_class$cpd_name)
  cpds_class$simple_cpd_name <- as.character(cpds_class$simple_cpd_name)
  cpds_class$class <- as.character(cpds_class$class)
  
  #check for missing compounds and subset to compounds found in dataset
  rownames(prop) <- gsub("'", "", rownames(prop))
  cpds_class <- subset(cpds_class, cpd_name %in% rownames(prop))
  
  missing_cpds <- rownames(prop)[!rownames(prop) %in% cpds_class$cpd_name]
  if(length(missing_cpds) > 0){
    print('These compounds were found in the data set but not in the reference compound list:')
    print(missing_cpds)
    stop('Compounds missing from the reference list.')
  }
  
  match_index <- match(cpds_class$cpd_name,rownames(prop))
  cpds_class <- cpds_class[order(match_index),]
  
  
  #sort compounds alphabetically within class
  
  #using "simple" name and not NIST name?
  if(simple_names){
    newOrder0 <- order(as.character(cpds_class$class), as.character(cpds_class$simple_cpd_name))
  } else {
    newOrder0 <- order(as.character(cpds_class$class), as.character(cpds_class$cpd_name))
  }
  cpds_class <- cpds_class[newOrder0,]
  
  #put 'special compounds' at top
  if(!is.null(special_cpds)){
    colours_sort <- read.table(special_cpds,
                               sep='\t',
                               comment.char='')
    colours_sort$V1 <- as.character(colours_sort$V1)
    colours_sort <- colours_sort[colours_sort$V1 %in% rownames(prop),]
    
    newOrder1 <- c(match(colours_sort$V1,cpds_class$cpd_name), which(!cpds_class$cpd_name %in% colours_sort$V1))
    cpds_class <- cpds_class[newOrder1,]
  }

  #Define compound colours##################
  
  #combine some classes
  cpds_class$class[cpds_class$class == 'alkane'] <- 'alk'
  cpds_class$class[cpds_class$class == 'alkene'] <- 'alk'
  cpds_class$class[cpds_class$class == 'ketone'] <- 'ketone_ester_ether'
  cpds_class$class[cpds_class$class == 'ester'] <- 'ketone_ester_ether'
  cpds_class$class[cpds_class$class == 'ether'] <- 'ketone_ester_ether'
  cpds_class$class[cpds_class$class == 'acid'] <- 'ketone_ester_ether'
  # cpds_class$primary_class[cpds_class$primary_class == 'nitrogenous'] <- 'NS'
  # cpds_class$primary_class[cpds_class$primary_class == 'sulfurous'] <- 'NS'
  
  
  hue_table <- data.frame(cbind(c('alcohol','aldehyde','alk','aromatic',
                                  'ketone_ester_ether','NS','terpene','other'),
                                c(0.686, 0.03, 0.85, 0.52, 0.77, 0.4, 0.2, 0)))
  names(hue_table)<-c('class','hue')
  hue_table$hue <- as.numeric(as.character(hue_table$hue))
  
  cpds_class$h <- 0
  cpds_class$s <- 0
  cpds_class$v <- 0
  
  set.seed(20)
  
  mean_sat <- 0.8
  mean_val <- 0.9
  
  for(i in 1:nrow(cpds_class)){
    hue <- as.numeric(
      as.character(hue_table[hue_table$class==cpds_class$class[i],]$hue))
    cpds_class$h[i] <- jitter(hue,3)
    cpds_class$s[i] <- jitter(mean_sat,4)
    cpds_class$v[i] <- jitter(mean_val,4)
  }
  
  hsv_vals <- cpds_class[,(ncol(cpds_class)-2):ncol(cpds_class)]
  hsv_vals[hsv_vals<0] <- 0
  hsv_vals[hsv_vals>1] <- 1
  
  #make "other" compounds dark grey
  hsv_vals$s[cpds_class$class=='other'] <- 0
  hsv_vals$v[cpds_class$class=='other'] <- hsv_vals$v[cpds_class$class=='other']/2
  
  cpds_class[,(ncol(cpds_class)-2):ncol(cpds_class)] <- hsv_vals
  
  
  cpds_class$colours <- apply(cpds_class[,(ncol(cpds_class)-2):ncol(cpds_class)],
                              MARGIN = 1,
                              FUN = function(x) hsv(x[1],x[2],x[3]))
  
  if(!is.null(special_cpds)){
    cpds_class$colours[1:dim(colours_sort)[1]] <- as.character(colours_sort$V2)
  }
  
  if(simple_names){
    output <- as.data.frame(with(cpds_class[1:nrow(prop),], cbind(cpd_name, simple_cpd_name, colours)))
  } else {
    output <- as.data.frame(with(cpds_class[1:nrow(prop),], cbind(cpd_name, colours)))
  }
  
  output[] <- lapply(output, as.character)
  return(output)
}

static_barplot <- function(prop, cols, legend_df){
  par(mar=c(1,10,1,1), xpd=T)
  barplot(prop, col=cols, horiz=T, las=2,
          xlim=c(0,1.5), xaxt='n', border=NA,cex.names=0.5)
  
  with(legend_df, legend(1,ncol(prop), # position
                         legend = cpd,
                         fill = hex,
                         border = NA,
                         bty = 'n',
                         xjust = 0,
                         yjust = 1,
                         y.intersp = 0.8,
                         cex = 0.5))
}

pdf_barplot <- function(prop, cols, legend_df, filename){
  
  pdf(paste0(filename, '_barplot.pdf'),
      width=8.8, height=ncol(prop)*0.13 + 1.3)
  
  static_barplot(prop, cols, legend_df)
  
  dev.off()
  
}

interactive_barplot <- function(prop, colour_df, filename, use_simple_names = FALSE){
  cols <- colour_df$colours
  
  #melt data frame
  long_props <- reshape2::melt(as.matrix(prop))
  colnames(long_props) <- c('cpd', 'sample', 'estimate')
  
  if(use_simple_names){
    long_props$cpd <- colour_df$simple_cpd_name[match(long_props$cpd, colour_df$cpd_name)]
    long_props$cpd <- factor(long_props$cpd, levels = as.character(colour_df$simple_cpd_name))
  }
  
  hc <- hchart(long_props, 'column',
               hcaes(x = 'sample', y = 'estimate', group = 'cpd'),
               stacking = 'normal',
               pointPadding = 0,
               groupPadding = 0.1,
               animation = F
  ) %>%
    hc_legend(enabled = T,
              itemStyle = list(fontSize = '11px', fontWeight = 'lighter'),
              layout = 'vertical',
              align = 'right',
              maxHeight = 700,
              y = 70,
              verticalAlign = 'top'
    ) %>%
    hc_plotOptions(
      series = list(states = list(hover = list(borderColor = 'gray',
                                               borderWidth = 1)),
                    point = list(events = list(mouseOver = JS("function() {
                                                              const point = this;
                                                              point.series.points.forEach(p => {
                                                               p.setState('hover');
                                                              });
                                                            }"),
                                               mouseOut = JS("function() {
                                                             const point = this;
                                                             point.series.points.forEach(p => {
                                                              p.setState('normal');
                                                             });
                                                            }")
                    ))
      ),
      column = list(borderWidth = 0)) %>%
    hc_colors(cols) %>%
    hc_xAxis(title = list(text = ''),
             #labels = list(enabled = F), #to get rid of labels entirely
             labels = list(step = 1, rotation = 45)) %>%
    hc_yAxis(
      visible = F,
      labels = list(enabled = F))
  
  htmlwidgets::saveWidget(hc, file = paste0(filename,'_barplot.html'))
}


easy_start <- function(data_dir, my_prefix, analysis_dir = '.'){
  ###Set working directory
  setwd(analysis_dir)
  
  ###Get library
  #baseline of chromatogram, calculated from the mean of several conditioned tubes,
  #then passed through a low-pass filter
  #see 2022-01-04_mzR_get_background.R for details
  load('smoothed.RData')
  
  #custom library with retention times
  #see 2022-11-15_build_compound_library.R for details
  load('2022-11-17_custom_lib4.Rdata')
  
  #component distinctiveness, for estimating compound abundances
  distinctiveness_list <- get_distinctiveness_list(custom_lib4, custom_lib_spectra4, smoothed)
  
  
  ###Load and process data
  library(xcms) #3.6.2
  library(magrittr)
  
  data1 <- runXCMS(data_dir) %>%
    clean_sample_names() %>%
    round_mz_rt()

  saveRDS(data1, paste0(my_prefix, '_raw_XCMS_data.rds'))
  
  ###Some minor automatic RT adjustment
  custom_lib4 <- finetune_RTs(custom_lib = custom_lib4,
                              custom_lib_spectra = custom_lib_spectra4,
                              data = data1)
  
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
  colour_df <- readRDS('colour_df.rds')
  
  cpd_sample_matrix <- combine_estimates(round1_output$abund_estimates, round2_output$abund_estimates) %>%
    long_to_wide() %>%
    remove_nonbio() %>%
    relative_abund() %>%
    reorder_matrix(axis = 'row',
                   order = colour_df$cpd_name)
  
  ###Plotting
  static_barplot(cpd_sample_matrix,
                 cols = colour_df$colours,
                 legend_df = readRDS('legend_df.rds'))
  
  pdf_barplot(cpd_sample_matrix,
              cols = colour_df$colours,
              legend_df = readRDS('legend_df.rds'),
              filename = my_prefix)
  
  library(highcharter)
  interactive_barplot(cpd_sample_matrix,
                      colour_df = colour_df,
                      filename = my_prefix)
  
  ###Output data matrix
  write.table(cpd_sample_matrix, file = paste0(my_prefix, '_output_matrix.tab'),
              sep = '\t',
              row.names = TRUE,
              col.names=NA)
  
  return(cpd_sample_matrix)
}


plot_ms <- function(spectrum1, spectrum2, title = ''){
  #useful function for plotting two spectra tail to tail for comparison
  spectrum1 <- spectrum1/sum(spectrum1)
  spectrum2 <- spectrum2/sum(spectrum2)
  plot(spectrum1, type = "h", main=title, cex.main=0.6,
       ylim=c(-max(spectrum2), max(spectrum1)))
  points(-spectrum2, type = "h", col="red")
  }
