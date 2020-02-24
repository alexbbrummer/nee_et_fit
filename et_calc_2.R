et.calc <- function(){
  ### This function is meant to perform a linear and non-linear fitting of h2o absorpotion data from LiCOR machines to provide estimates of evapotranspiraiton.  It is modeled off of a previous function designed to perform similar fitting of co2 data to estimate net ecosystem exchange.  The nee.calc code still exists as a future step will be to integrate these two fits to be performed simultaneously (the nee.fit function assumes that h20 remains constant).  All comments from previous version/s are maintained.
  
  ### Transpiration is only plants when stomata are open (daytime measurements).  Evapotranspiration is everything, plants and soils (daytime measurements), and evaporation is only soils, (nighttime or dark/respiration measurements).  Thus, we are measuring, directly, only evapotranspiration (ET) and evaporation (E).  But, we can extract Transpiration as ET = E + T (see Wang et al.).  Some report T/ET ratios.  Our estimates appear reasonable with others reported by Wang et al. and Burba et al.  Do see Wang et al.'s review paper for comments on reliability of nightime measurements as proxy for evaporation.  Specifically, he quotes studies that found measurable non-zero transpiration during plant respiration.
  
  ## This version is from 02.
  ## This function is meant to perform the linear and non-linear fitting to nee data.  It will produce a plot of the data, queery the user if they would like to modify the time interval over which the data is fit, then print the paste statement of data values.  Currently, the function is written so that it works for a given filename.  From there, we can extend to a folder/directory.
  
  ## These packages are necessary.
  library("proto")
  library("nls2")
  
  ####----Insert here the readline() command for querying user for the directory of LiCOR files.--------
  ## For querying user to define the working directory in which the LiCOR files are located.
  readline("Please set the working directory to the folder \n that contains the LiCOR files to be analyzed. \n Do so with the upcoming prompt. \n Note that you must choose a(any) file in the \n folder that you want to set as the working directory. \n Please press 'return' to continue.")
  setwd(dirname(file.choose()))
  
  ## Define directory as an object that contains the dir() item names as a vector.
  directory <- dir()
  
  ## For reading the .txt files, replace the ".txt" paste with ".txt"
  #   photo.names <- directory[grep(paste("[0-9]", ".txt", sep = ""), dir(), ignore.case = TRUE, value = FALSE)]
  #   ambient.names <- directory[grep(paste("[0-9]", "a", ".txt", sep = ""), dir(), ignore.case = TRUE, value = FALSE)]
  #   resp.names <- directory[grep("resp", dir(), ignore.case = TRUE, value = FALSE)]
  
  photo.names <- grep("[^resp].txt", grep("[^_a]\\.txt", dir(), value = TRUE), value = TRUE)
  ambient.names <- grep("a.txt", dir(), value = TRUE)
  resp.names <- grep("_[[:digit:]]resp.txt", dir(), value = TRUE)
  
  et.fit <- function(filename){
    ## For reading the .txt files, replace read.csv with read.table, and add a skip = 9 parameter to the function.
    input <- read.table(filename, header = FALSE, skip = 9)
    ## For reading the .txt files, these input subsetting commands may not be necessary.
    #     input <- input[-1,]
    #     input <- input[,-1]
    
    if(length(grep("resp", filename, ignore.case = TRUE, value = FALSE)) == 1){
      ambient <- read.table(paste(strsplit(filename, "resp.txt"), "a.txt", sep = ""), header = FALSE, skip = 9)
    } else{
      ambient <- read.table(paste(strsplit(filename, ".txt"), "a.txt", sep = ""), header = FALSE, skip = 9)
    }
    #     ambient <- ambient[,-1]
    #     ambient <- ambient[-1,]
    
    #  /// define constants - for Enquist Tent///
    vol = 2.197   # m^3, tent volume
    area = 1.69   # m^2, tent area
    R = 8.314472 	# J/mol K
    
    ## Define vectors to work with
    
    ## The data files are currently being read into R as factors, hence the awkward as.numeric(as.character()) class coercion.
    time <- as.numeric(as.character(input[,1])) #s
    co2 <- as.numeric(as.character(input[,8])) #umol/mol
    h2o <- as.numeric(as.character(input[,12])) #mmol/mol
    par <- as.numeric(as.character(input[,4])) # don't think there's any data here
    press <- as.numeric(as.character(input[,3])) #kPa
    temp <- as.numeric(as.character(input[,2])) #C
    
    #  /// average T and P per measurement ///
    tav <- mean(temp)
    pav <- mean(press)
    cav <- mean(co2)
    wav <- mean(h2o)
    # wav <- mean(h2o)  # notes from mtg with brian 1/21/20: Keep H20 here as is and leave NEE alone.  Separately, calculate H20 rates using similar approach as with C02.
    
    
    wprime <- h2o/(1-(h2o/10**3)) #dilution correction for gas mixture
    
    wamb <- mean(as.numeric(as.character(ambient[,12]))) #average h2o at ambient for leaky fit model.
    
    ## I think this is just the ideal gas law, but not clear to me why ideal gaw law is necessary for converting from whatever the old units were (find out and enter) to the new units.  The 18 represents the molar mass of H2O, 18 grams/mole.
    
    wamb <- wamb*R*(tav+273.15)/(18*pav)    # change to mmol/mol or ppm, or grams per mol?
    # camb
    
    #  /// Plotting the H20 vs. Time Curve //
    
    plot(wprime~(time), main = filename)
    
    ## Queery user for start time for fitting.  Default is set to 10 in the if() statement
    tstart <- readline("Enter preferred start time for fitting. \n Do not include units. \n Round to nearest integer second. \n Do not use 0. \n  If default of 10s is preferred, press 'return':")
    if(!grepl("^[0-9]+$", tstart)){
      tstart <- 10
    }
    tstart <- as.integer(tstart)
    
    ## Queery user for finish time for fitting.  Default is set to 80 in the if() statement
    tfinish <- readline("Enter preferred finish time for fitting. \n Do not include units. \n Round to nearest integer second. \n  If default of 80s is preferred, press 'return':")
    if(!grepl("^[0-9]+$", tfinish)){
      tfinish <- 80
    }
    tfinish <- as.integer(tfinish)
    
    ## Linear fitting code.
    
    linear.fit <- lm(wprime[tstart:tfinish]~(time[tstart:tfinish]))
    
    ## AIC value for linear model for later comparison with non-linear model.
    aic.lm <- AIC(linear.fit)
    
    # Calculate intercept
    inter<- as.numeric(linear.fit$coeff[1])
    
    # Calculate slope
    dwdt <- as.numeric(linear.fit$coeff[2])
    
    # Calculate r-squared (we're not reporting chi-squared significance from the non-linear fit, so this may not be so necessary)
    rsqd <- summary(linear.fit)$r.sq
    
    # Calculate flux from linear model
    flux_lm <- (vol*pav*(10**3)*dwdt) / (R*area*(tav + 273.15))	# in mmol/m2/s
    
    # Make plot of line.
    abline(inter,dwdt, col=6)
    
    ## comment/uncomment to use print statement for including non-linear fitting with aic scores
    ## Non-linear fitting code.
    # Set cnot to the actual first value of cprime used in the fitting time domain.
    wnot = wprime[tstart]

    # Define a temporary data frame from which the functional variables come from.
    df = data.frame(wprime, time)

    # Define a subset category from the tstart and tfinish variables.
    subsettime <- time > tstart & time < tfinish

    # Define boundaries of parameter grid.  May need to modify for evapotranspiration.
    strt <- data.frame(A = c(150, 850), B = c(0, 1000))

    # Use nls2() to scan through parameter grid, searching for "best" actual starting points.  control variable is set to prevent warnings from ending loop.
    optimize.start <- nls2(wprime ~ (wnot - A)*exp(-time/B) + A, data = df, start=strt, subset = subsettime, algorithm = "brute-force", control = nls.control(warnOnly = TRUE), trace = FALSE) #(A=375, B=40)

    # Run nls() with the optimized starting values from previous nls2().  Control variable is set to prevent warnings from ending loop.  However, they will still be printed at end of run.  When this happens, it is indicative of the fact that the function parameters (A and B) are large (non-physical) for the fitting, yet still produce a fit.  This is worth further investigation.  However, it appears that the nee value produced by the exponential model in such circumstances does not deviate from the linear model by much more than half a percent.  Add a "trace = TRUE" parameter setting to the nls() function to be able to watch the values of A and B change with each iteration.  A is Css and B is tau, from Saleska et al., and to translate to variables extracted from fit further down.
    uptake.fm <- nls(wprime ~ (wnot - A)*exp(-time/B) + A, data = df, start = coef(optimize.start), subset = subsettime, control = nls.control(warnOnly = TRUE), trace = FALSE)

    ##
    sigma <- summary(uptake.fm)$sigma

    ## AIC value for non-linear model for later comparison with linear model.
    aic.nlm <- AIC(uptake.fm)

    Wss = summary(uptake.fm)$param[1]
    tau = summary(uptake.fm)$param[2]
    # print(c(Wss, tau))
    flux_exp <- -((wamb-Wss)/(area*tau))*(vol*pav*(10**3-wav)/(R*(tav + 273.15))) #equation 4 in Saleska 1999

    curve((wnot - Wss)*exp(-(x-time[tstart])/tau) + Wss, col = 4, add = TRUE)	#equation 3 in Saleska 1999 to plot for visual inspection.##

    
    time <- "Evapotranspiration"
    
    if(length(grep("resp", filename, ignore.case = TRUE, value = FALSE)) == 1){
      time <- "Evaporation"
    }
    
    print(data.frame("tstart" = tstart, "tfinish" = tfinish, "time" = time, "flux_lm" = flux_lm, "flux_nlm" = flux_exp, "rsqd" = rsqd, "nlm_sigma" = sigma, "aic.lm" = aic.lm, "aic.nlm" = aic.nlm))
    # comment/uncomment to use print statement for including non-linear fitting with aic scores
    # print(data.frame("tstart" = tstart, "tfinish" = tfinish, "time" = time, "flux_lm" = flux_lm, "flux_exp" = flux_exp, "rsqd" = rsqd, "sigma" = sigma, "aic.lm" = aic.lm, "aic.nlm" = aic.nlm))
    
    replicate <- readline("Would you like to redo the fitting with \n a different time domain? (y/n)")
    if(replicate == "y"){
      et.fit(filename)
    } else {
      return(c(tstart,tfinish,time,wamb,tav,pav,cav,flux_lm,flux_exp,rsqd,sigma,aic.lm,aic.nlm))
      # comment/uncomment to use print statement for including non-linear fitting with aic scores
      # return(c(tstart,tfinish,time,wamb,tav,pav,cav,flux_lm,flux_exp,rsqd,sigma,aic.lm, aic.nlm))
    }
  }
  
  
  
  
  stats.df <- c()
  
  for (i in 1:length(photo.names)){
    stats.df <- rbind(stats.df, et.fit(photo.names[i]))
    
  }
  
  if (length(resp.names) > 1){
    for (i in 1:length(resp.names)){
      stats.df <- rbind(stats.df, et.fit(resp.names[i]))
    }
  }
  
  stats.df <- as.data.frame(stats.df)
  names.vec <- c("tstart", "tfinish", "time", "wamb", "tav", "pav", "cav", "flux_lm", "flux_nlm", "LM rsqd", "non-linear sigma", "aic_lm", "aic_nlm")
  # comment/uncomment to use print statement for including non-linear fitting with aic scores
  # names.vec <- c("tstart", "tfinish", "time", "wamb", "tav", "pav", "cav", "flux_lm", "flux_exp", "LM rsqd", "non-linear sigma", "aic_lm", "aic_nlm")
  for(i in 1:length(names.vec)){
    names(stats.df)[i] <- names.vec[i]
  }
  
  stats.df
  write.csv(stats.df, file = paste(paste(strsplit(getwd(), "/")[[1]][length(strsplit(getwd(), "/")[[1]])], "summary", sep = " "), ".csv", sep = ""))
  
}