library ("seewave");
library ("signal")
source ("data_generator.R")

########################################################################
findGmodes = function(r, um, dm){
########################################################################
  # data must not contain zeros, so 'r' must not contain NA values
  # finding g-modes

  # r  : output from specPdgrm
  # um : up   - to define the neighborhood
  # dm : down - to define the neighborhood
  
  x = r$x; y = r$y; z = r$z; # values from spectrogram
  
  ps = NULL;
  mp = which.max(z[dim(z)[1],]); # position maximum value

  for(i in (dim(z)[1]-1):1){
    
    #print(i);
    ps  = c(ps, mp);
    ind = (mp-dm):(mp+um);
    ind = ind[ind>0]; # preventing values below 0
    ind = ind[ind<dim(z)[2]]; # preventing from large values
    mp  = ind[which.max(z[i, ind])];
  }
  
  ps   = c(ps, mp);
  ps   = ps[length(ps):1]
  maxf = y[ps]; # FREQUENCIES FOR LARGEST POWER
  return(maxf);

}

########################################################################
movf = function(vec, n, f){
########################################################################

  if(n == 1){
    return(vec);
  }
  
  N   = length(vec);
  out = rep(NA, N);
  
  if( n %% 2 != 0){ # odd
    aux  = (n+1)/2;
    aux1 = aux - 1; 
    
    for(i in aux:(N-aux+1)){
      out[i] = f(vec[(i - aux1):(i + aux1)]);
    }
    
    index  = seq(1, n - 2, by = 2);
    index1 = seq(N - n + 1 + 2, N, by = 2);
    
    for(i in 1:length(index)){
      
      out[i] = f(vec[1:index[i]]); # calculating first values
      
      out[i + N - aux + 1] = f(vec[index1[i]:N]);
      
    }
    
  }else{
    stop("This version only works with odd n");
  }
  
  return(out);
  
}

    

########################################################################
ints = function(n, l, p = 00, eq = TRUE){
########################################################################
  # Function taken from psplinePsd function
  
  #' This function produces a matrix which each row contains the first and last indexes
  #'  to split a time series of length "n"
  #'
  #' n  = time series length
  #' l  = interval length
  #' p  = overlapping percentage
  #' eq = last interval has length "l"
  
  if( (n<=0) || (n %% 1 != 0)) stop("n must be a positive integer")
  if( (l<=0) || (l %% 1 != 0)) stop("l must be a positive integer")
  if( (p<0) || (p>=100)) stop("p must be a positive number between 0 and 100")
  if( l >= n ) stop("l must be lower than n")
  
  # This version yields eve interval lengths
  if (l %% 2 != 0) stop("this version of bsplinePsd must have l even")
  
  ovp  = round(l*p/100,0); # number of overlaped points
  a    = l - ovp + 1;
  col1 = c(1, seq(from = a, to = n - l + a -1, by = a-1));
  
  index = cbind(col1, col1 + l-1);
  # checking the lengths
  # index[,2] - index[,1] + 1 == l
  # diff(index[,1]) == a-1
  # diff(index[,2]) == a-1
  
  # fixing last interval
  index[dim(index)[1], 2] = n;
  
  colnames(index) = NULL;
  rownames(index) = NULL;
  
  #cat(paste("The number of data subsets is ",dim(index)[1], sep=""),"\n");
  
  if(eq == TRUE){
    
    lint = dim(index)[1]
    index[lint, 1] = n - l + 1;
    # index[,2] - index[,1] == l
    
    x = index[lint-1, 2] - index[lint, 1] + 1;
    x = round(100 * x / l, 2);
    
    if(x!=p){
      #cat(paste("Last interval overlaps ", x, "%", sep = ""), "\n");
    }
    
  }else{
    
    aux = index[dim(index)[1], 2] - index[dim(index)[1], 1] + 1;
    
    if(aux %% 2 != 0){
      
      # last interval is been made even
      index[dim(index)[1], 1] = index[dim(index)[1], 1] - 1;
      aux = index[dim(index)[1], 2] - index[dim(index)[1], 1] + 1;
    }
    
    if(aux != l){
      cat(paste("Last interval contains ", aux, " observations", sep = ""), "\n");
    }
  }
  
  return(index);
}

########################################################################
covpbb = function(data, mod, l=200, p=90, fs=16384, movGmode = 11, 
                  um = 3, dm = 3, movBand = 5, timeGmode = NULL, 
                  Obergaunlinguer_data, actPlot = FALSE){
########################################################################  
  # data     : dataset as matrix (with no zeros [specPdgrm])
  # mod      : model
  # l        : interval length in spectrogram
  # p        : overlaping percentage in spectrogram
  # movGmode : number of steps to smooth estimated g-modes 
  # um       : define upper neighborhood to find g-modes 
  # dm       : define lower neighborhood to find g-modes
  # movBand  : define the number of points to smooth the band
  # timeGMode: time interval to define g-modes
  # Oberganlinger_data: simulated R and M time evolution
  
  # Compute true ratios
  true_ratios = Obergaunlinguer_data$Mpns / (Obergaunlinguer_data$Rpns^(2))
  
  # spectrogram
  r = specPdgrm(data$V2, data$V1, l=l, p=p, fs=fs, actPlot=FALSE, logPow=TRUE,
                zoomFreq=c(0,1)); # generating the spectrogram
  
  n = length(data$V1);
  
  # starting & ending points of the intervals used in the spectrogram
  index = ints(n=n, l=l, p=p); # from psplinePsd
  
  # centred point for even "l"
  mindx = index + cbind(rep(l/2 -1, dim(index)[1]), rep(-(l/2-1), dim(index)[1]));
  
  # Ajusting data
  timedata = seq(from=data$V1[1], to=data$V1[n], length = n);
  
  # mean time (centred point) for our g-mode estimates
  timefreq = apply(mindx, 1, function(x) mean(timedata[c(x[1], x[2])]) );
  
  # g-modes
  if( !is.null(timeGmode)){
    #timeGmode = data0[c(1,length(data0[,1])), 1];
    out = apply(as.matrix(timefreq), 1, 
                function(x){
                  if((x >= timeGmode[1]) & (x <= timeGmode[2])){
                    return(TRUE);
                  }else{
                    return(FALSE);
                  }
                });
    
    #maxf     = maxf[out];
    timefreq = timefreq[out];
    
    r$x = r$x[out];
    r$z = r$z[out, ];
  }
  
  maxf = findGmodes(r, um=um, dm=dm);
  maxf = movf(maxf, movGmode, median); # smoothing g-mode estimates
  #points(timefreq,maxf, col = 'green')
  
  # prediction : pred$fit pred$lwr pred$upr
  new  = data.frame(f = maxf);
  pred = predict(mod, new, interval = "prediction"); # predictions
  
  #plot(new$f,pred[,1])
  #pred[pred[,2] < 0, 2] = 0; # discarding negative values in CI
  
  ### generating band function ### 
  # n = 1 is equal to system defined above
  
  # interpolating lower bound for predicted values
  fd = approxfun(x = timefreq, y = movf(pred[,2],n=movBand,median), method = "linear",
                 yleft = NA, yright = NA, rule = 1, f = 0, ties = "mean");
  # interpolating upper bound for predicted values
  fu = approxfun(x = timefreq, y = movf(pred[,3],n=movBand,median), method = "linear",
                 yleft = NA, yright = NA, rule = 1, f = 0, ties = "mean");
  
  # fd & fu use smooth confidence intervals by using "movf"
  aux = cbind(true_ratios,                    # true ratios
              fd(Obergaunlinguer_data$time),  # lower band (using predicted ratios)
              fu(Obergaunlinguer_data$time)); # upper band (using predicted ratios)
  
  # discarding the true values which are out of the range of the predicted values
  
  if(actPlot == TRUE){
    
    plot(new$f,pred[,1])
  
    points(timefreq,maxf, col = 'green')
  
    plot(Obergaunlinguer_data$time,aux[,1],xlab="time[s]",ylab="r",ylim=c(-0.000515,0.0037809),pch=1)
    points(Obergaunlinguer_data$time,aux[,2],col="red",pch=2)
    points(Obergaunlinguer_data$time,aux[,3],col="red",pch=3)
    leg <- c("true ratio", "pred lower","pred upper")
    col=c("black","red","red")
    legend(x=.25,y=0.0037,legend=leg,cex=.8,col=col,pch=c(1,2,3))
  
  # From plotOmeSim
    yaux = c(true_ratios, pred[,2:3]);
    plot(Obergaunlinguer_data$time, true_ratios, xlab = "Time",
       ylab = "Ratio", ylim = c(min(yaux), max(yaux)), type = "n");
    arrows(timefreq, pred[,2], timefreq, pred[,3], code=3, angle=90,
         length=0.05, col="gray",pch=3);
    points(Obergaunlinguer_data$time, true_ratios, col = "black",pch=1);
    points(timefreq, pred[,1], col = "red", cex = pred[,1]/max(pred[,1])+ 0.3,pch=2);
  
    leg <- c("true ratio", "pred ","pred uncertainty")
    col=c("black","red","gray")
    legend(x=.25,y=0.0037,legend=leg,cex=.8,col=col,pch=c(1,2,3))
  }
  
  # left side
  disc = which(is.na(aux[,2]));
  
  if(length(disc) != 0){
    aux = aux[-disc,];
  }
  
  # right side
  disc = which(is.na(aux[,3]));
  
  if(length(disc) != 0){
    aux = aux[-disc,];
  }
  
  # testing if the true ratios are inside the bands
  # testing if the true ratios are inside the bands
  prop = apply(aux, 1, 
               function(x){
                 if((x[1]>= x[2]) && (x[1] <= x[3])){
                   return(1);
                 }else{
                   return(0);
                 }
               });
  
  aux[aux[,2]<0,2] = 0; # it replaces negative values in lower limit
  
  l = aux[,3] - aux[,2];
  p = mean(prop);
  
  return(list(covpbb = p, medBandWidth = median(l)));
  
}


########################################################################
repcovpbb = function(wvf, duration, ampl, fcut,
                     mod, N, snr = NULL, movGmode = 11, um = 3, dm = 3,
                     movBand = 5, timeGmode = NULL, l=200, p=90, fs=16384,
                     Obergaunlinguer_data){
########################################################################  
  # data : matrix (time, signal)
  # mod  : lm object
  # N    : number of replications
  # ampl : distances (scalar or vector)

  if(length(ampl) == 1){
    
    out_cp = NULL;
    out_bl = NULL;
    
    for(i in 1:N){
      noisydata = data_generator(fs, duration, wvf, ampl, fcut, actPlot=FALSE);
      noisydata = data.frame("V1"=noisydata$t,"V2"=noisydata$y);
      
      aux = covpbb(noisydata, mod, l, p, fs, movGmode, 
                   um, dm, movBand, timeGmode, 
                   Obergaunlinguer_data, actPlot = FALSE);
      #print(aux$covpbb)        
      out_cp = c(out_cp, aux$covpbb);
      out_bl = c(out_bl, aux$medBandWidth);      
    }
    
    out = list(covpbb = out_cp, medBandWidth = out_bl);
    
  }else if(length(ampl) > 1){
    
    out_cp = list();
    out_bl = list();
    
    for(j in ampl){
    
      aux_cp = NULL; # to save coverage probabilities
      aux_bl = NULL; # to save band widths
      
      for(i in 1:N){
        
        noisydata = data_generator(fs, duration, wvf, ampl = j, fcut, actPlot=FALSE);
        noisydata = data.frame("V1"=noisydata$t,"V2"=noisydata$y);
        
        aux = covpbb(noisydata, mod, l, p, fs, movGmode, 
                     um, dm, movBand, timeGmode, 
                     Obergaunlinguer_data, actPlot = FALSE);
        #print(aux$covpbb)        
        aux_cp = c(aux_cp, aux$covpbb);
        aux_bl = c(aux_bl, aux$medBandWidth);
        
      }
      if(is.null(snr)){
        
        out_cp[[paste(j)]] = aux_cp;
        out_bl[[paste(j)]] = aux_bl;
        
      }else{
        
        out_cp[[paste(snr * j)]] = aux_cp;
        out_bl[[paste(snr * j)]] = aux_bl;  
        
      }
    }
    
    out = list(covpbb = out_cp, medBandWidth = out_bl);
    
  }
  
  return(out)
  
}

















##### Functions not yet checked by MAB




repcovpbbWN = function(data, mod, N, SDs, movGmode = 11, um, dm, movBand, 
                       timeGmode, l=200, p=90, fs=16384){
  
  # calculate coverage probabilities for the standarized data (+ white noise)
  
  newdata = data$V2;
  newdata = (newdata - mean(newdata)) / sd(newdata);
  cat("The data has been standarised", "\n")
  
  n       = length(newdata);
  
  if(length(SDs) == 1){
    
    out = NULL;
    
    for(i in 1:N){
      
      # Generate noise
      Y = rnorm(n, sd = SDs);
      
      # noisy data
      noisydata$V2 = newdata + Y; # signal + noise
      hpf = 50; # high pass filter
      noisydata$V2 = ffilter(noisydata$V2, f = fs, from = hpf);
      
      out = c(out, covpbb(noisydata, mod=mod, l=l, p=p, fs=fs, movGmode=movGmode, 
                          um=um, dm=dm, movBand=movBand, timeGmode=timeGmode));
      
    }
    
  }else if(length(SDs) > 1){
    
    out = list();
    
    for(s in SDs){
      
      cat(paste("sd =", s, "\n"));
      
      aux = NULL;
      
      for(i in 1:N){
        
        # Generate noise
        Y = rnorm(n, sd = s);
        
        # noisy data
        noisydata$V2 = newdata + Y; # signal + noise
        hpf = 50; # high pass filter
        noisydata$V2 = ffilter(noisydata$V2, f = fs, from = hpf);
        
        aux = c(aux, covpbb(noisydata, mod=mod, l=l, p=p, fs=fs, movGmode=movGmode, 
                            um=um, dm=dm, movBand=movBand, timeGmode=timeGmode));
        
      }
      
      out[[paste(s)]] = aux;
      
    }
  }
  
  return(out)
  
}

##########################
### BAYESIAN functions ###
##########################

covpbbBayes = function(data, postSamples, l=200, p=90, fs=16384, movGmode = 11, 
                       um = 3, dm = 3, movBand = 5, timeGmode = NULL,
                       Obergaunlinguer_data, actPlot = FALSE, fmean = 0){
  
  # The only difference with 'covpbb' function is how the predicted 
  #  values are calculated.
  
  # data     : dataset as matrix (with no zeros [specPdgrm])
  # mod      : model
  # l        : interval length in spectrogram
  # p        : overlaping percentage in spectrogram
  # movGmode : number of steps to smooth estimated g-modes 
  # um       : define upper neighborhood to find g-modes 
  # dm       : define lower neighborhood to find g-modes
  # movBand  : define the number of points to smooth the band
  # timeGMode: time interval to define g-modes
  # Oberganlinger_data: simulated R and M time evolution
  # fmean: mean of frequencies used to generate the model
  
  # Compute true ratios
  true_ratios = Obergaunlinguer_data$Mpns / (Obergaunlinguer_data$Rpns^(2));
  
  # spectrogram
  r = specPdgrm(data$V2, data$V1, l = l, p = p, fs = fs, actPlot = FALSE, logPow = T,
                zoomFreq=c(0,1)); # generating the spectrogram
  
  n     = length(data$V1);
  # starting & ending points of the intervals used in the spectrogram
  index = ints(n = n, l = l, p = p); # from psplinePsd 
  
  # centred point for even "l"
  mindx = index + cbind(rep(l/2 -1, dim(index)[1]), rep(-(l/2-1), dim(index)[1]));
  # Ajusting data / #it can be directly taken from the 'data'
  timedata = seq(from=data$V1[1], to=data$V1[n], length = n);
  
  # mean time (centred point) for our g-mode estimates
  timefreq = apply(mindx, 1, function(x) mean(timedata[c(x[1], x[2])]) );
  
  # g-modes
  
  if( !is.null(timeGmode)){
    
    # g-mode time range
    
    out = apply(as.matrix(timefreq), 1, 
                function(x){
                  if((x >= timeGmode[1]) & (x <= timeGmode[2])){
                    return(TRUE);
                  }else{
                    return(FALSE);
                  }
                });
    
    timefreq = timefreq[out];
    
    r$x = r$x[out];   # discarding values out of g-mode range
    r$z = r$z[out, ];
    
  }
  
  maxf = findGmodes(r, um=um, dm=dm);
  maxf = movf(maxf, movGmode, median); # smoothing g-mode estimates
  
  cmaxf = maxf - fmean; # centred g-modes (SEE MODEL USED IN JAGS)
  
  xs =  rbind(rep(1, length = length(cmaxf)), cmaxf, cmaxf^2, cmaxf^3); 
  
  # postSamples = (a0, a1, a2, a3, b1); # postSample structure
  postSamples = as.matrix(postSamples); # needed for matrix operations below
  
  # a0+ a1*(fnew[i]-mean(f))+a2*(fnew[i]-mean(f))^2+a3*(fnew[i]-mean(f))^3
  mus = postSamples[,1:4] %*% xs; # means
  
  # standard deviations
  sds = as.matrix(postSamples[,5], ncol=1) %*% t(as.matrix(maxf,ncol=1));
  sds = sqrt(sds);
  
  # simulation from Normal(mus, sds)
  n_maxf = length(maxf);
  X      = apply(cbind(mus, sds), 1, function(x)
    rnorm(n_maxf, mean = x[1:n_maxf], sd = x[-(1:n_maxf)]));
  X      = t(X);
  
  # median, quantiles 0.025 & 0.975
  pred = apply(X, 2, function(x) quantile(x ,probs = c(0.5,0.025,0.975)));
  pred = t(pred);
  
  #pred[pred[,2] < 0, 2] = 0; # discarding negative values in CI
  
  ### generating band function ### 
  
  # interpolating lower bound for predicted values
  fd = approxfun(x = timefreq, y = movf(pred[,2],n=movBand,median), method = "linear",
                 yleft = NA, yright = NA, rule = 1, f = 0, ties = "mean");
  # interpolating upper bound for predicted values
  fu = approxfun(x = timefreq, y = movf(pred[,3],n=movBand,median), method = "linear",
                 yleft = NA, yright = NA, rule = 1, f = 0, ties = "mean");
  
  # fd & fu use smooth confidence intervals by using "movf"
  
  aux = cbind(true_ratios,                    # true ratios
              fd(Obergaunlinguer_data$time),  # lower band (using predicted ratios)
              fu(Obergaunlinguer_data$time)); # upper band (using predicted ratios)
  
  # discarding the true values which are out of the range of the predicted values
  
  discL = which(is.na(aux[,2])); # left side
  discR = which(is.na(aux[,3])); # right side
  
  disc  = unique(c(discL, discR)); 
  
  if(length(disc) != 0){
    
    aux = aux[-disc, ];
    
    trueRatioTime = Obergaunlinguer_data$time[-disc];
    
  }
  
  out = NULL;
  
  if( !is.null(timeGmode)){
    
    # To test only true values inside g-mode time range
    
    out = apply(as.matrix(trueRatioTime), 1,
                function(x){
                  if((x >= timeGmode[1]) & (x <= timeGmode[2])){
                    return(TRUE);
                  }else{
                    return(FALSE);
                  }
                });
    
    aux = aux[out, ];
    
  }
  
  if(actPlot == TRUE){
    
    plot(Obergaunlinguer_data$time, true_ratios, 
         ylim=c(min(pred[,2]),max(pred[,3])), xlab = "Time", ylab = "Ratio",
         main = "Bayesian analysis");
    arrows(timefreq, pred[,2], timefreq, pred[,3], 
           code=3, angle=90, length=0.05, col="gray", pch = 3);
    
    points(Obergaunlinguer_data$time, true_ratios, col = "black",pch=1);
    points(timefreq, pred[,1], col = "red", cex = pred[,1]/max(pred[,1])+ 0.3,pch=2);
    
    leg <- c("true ratio", "pred ","pred uncertainty")
    col=c("black","red","gray")
    legend(x=.25,y=0.0037,legend=leg,cex=.8,col=col,pch=c(1,2,3))
    
  }
  
  # testing if the true ratios are inside the bands
  prop = apply(aux, 1, 
               function(x){
                 if((x[1]>= x[2]) && (x[1] <= x[3])){
                   return(1);
                 }else{
                   return(0);
                 }
               });
  
  #aux[aux[,2]<0,2] = 0; # it replaces negative values in lower limit
  
  l = aux[,3] - aux[,2];
  p = mean(prop);
  
  return(list(covpbb = p, medBandWidth = median(l)));
  
}


repcovpbbBayes = function(data, postSamples, N, dists, snr = NULL, movGmode = 11, um = 3, dm=3,
                          movBand = 5, timeGmode = NULL, l=200, p=90, fs=16384){
  
  # data : matrix (time, signal)
  # mod  : lm object
  # N    : number of replications
  # dists: distances (scalar or vector)
  
  Fs   = fs;   # 16384; # sampling frequency
  mean = 0.0;
  std  = 1.0;
  n    = dim(data)[1]; # 16384
  
  freq    = Fs*np.fft.fftfreq(n);  # two-sided frequency vector
  psd2    = aLIGO_PSD(freq,2);     # two-sided PSD
  newpsd2 = aLIGO_PSD_new(freq,2); # two-sided PSD              
  
  #ff     = freq[1:int(n/2)];    # one-sided frequency vector 
  #psd    = aLIGO_PSD(ff,1);     # one-sided PSD              
  #newpsd = aLIGO_PSD_new(ff,1); # one-sided PSD 
  
  noisydata = data;
  
  if(length(dists) == 1){
    
    out = NULL;
    
    newdata = data$V2 * dists;
    
    for(i in 1:N){
      
      # Generate noise
      X   = rnorm(n, mean = mean, sd = std); # Gaussian white noise
      XX  = sqrt(Fs) * fft(X) / (n*sqrt(n)); # FFT computing and normalization
      XXX = XX*sqrt(newpsd2);                # Coloring
      Y   = fft(XXX, inverse = TRUE);        # FFT inverse
      Y   = Re(Y)*sqrt(n);
      
      # noisy data
      noisydata$V2 = newdata + Y; # signal + noise
      hpf = 50; # high pass filter
      noisydata$V2 = ffilter(noisydata$V2, f = fs, from = hpf);
   
      out = c(out, covpbbBayes(noisydata, postSamples, l=l, p=p, fs=fs, movGmode = movGmode, 
                               um = um, dm = dm, movBand = movBand, timeGmode = timeGmode));
      
    }
    
  }else if(length(dists) > 1){
    
    out = list();
    
    for(j in dists){
      
      print(j);
      
      aux = NULL;
      
      newdata = data$V2 * j;
      
      for(i in 1:N){
        
        # Generate noise
        X   = rnorm(n, mean = mean, sd = std); # Gaussian white noise
        XX  = sqrt(Fs) * fft(X) / (n*sqrt(n)); # FFT computing and normalization
        XXX = XX*sqrt(newpsd2);                # Coloring
        Y   = fft(XXX, inverse = TRUE);        # FFT inverse
        Y   = Re(Y)*sqrt(n);
        
        # noisy data
        noisydata$V2 = newdata + Y; # signal + noise
        hpf = 50; # high pass filter
        noisydata$V2 = ffilter(noisydata$V2, f = fs, from = hpf);
        
        aux = c(aux, covpbbBayes(noisydata, postSamples, l=l, p=p, fs=fs, movGmode = movGmode, 
                                 um = um, dm = dm, movBand = movBand, timeGmode = timeGmode));
        
      }
      
      if(is.null(snr)){
        
        out[[paste(j)]] = aux;
        
      }else{
        
        out[[paste(snr * j)]] = aux;  
        
      }
    }
  }
  
  return(out)
  
}

plotOneSim = function(data, mod = NULL, postSamples=NULL, dist=1, l=200, p=90,
                      fs=16384, movGmode = 11, um = 3, dm = 3, movBand = 5,
                      timeGmode = NULL, filterHz = NULL){
 
  # data: original data. The functions injects coloured noise

  Fs   = fs; # 16384; # sampling frequency
  mean = 0.0;
  std  = 1.0;
  n    = dim(data)[1]; # 16384
 
  freq    = Fs*fftfreq(n);  # two-sided frequency vector
  psd2    = aLIGO_PSD(freq,2);     # two-sided PSD
  newpsd2 = aLIGO_PSD_new(freq,2); # two-sided PSD             
  
  ff     = freq[1:int(n/2)];    # one-sided frequency vector
  psd    = aLIGO_PSD(ff,1);     # one-sided PSD             
  newpsd = aLIGO_PSD_new(ff,1); # one-sided PSD
  
  noisydata = data;
  noisydata$V2 = noisydata$V2 * dist;

  out = NULL;

  # Generate noise
  X   = rnorm(n, mean = mean, sd = std); # Gaussian white noise
  XX  = sqrt(Fs) * fft(X) / (n*sqrt(n)); # FFT computing and normalization
  XXX = XX*sqrt(newpsd2);                # Coloring
  Y   = fft(XXX, inverse = TRUE);        # FFT inverse
  Y   = Re(Y)*sqrt(n);

  # noisy data
  noisydata$V2 = noisydata$V2 + Y; # signal + noise
 
  if(!is.null(filterHz)){
    
    noisydata$V2 = ffilter(noisydata$V2, f = fs, from = filterHz);
   
  }
 
  r = specPdgrm(noisydata$V2, l = l, p = p, fs = fs, actPlot = FALSE, logPow = T,
                zoomFreq=c(0,1)); # generating the spectrogram
 
  n     = length(data$V1);
  # starting & ending points of the intervals used in the spectrogram
  index = ints(n = n,l = l, p = p); # from psplinePsd
 
  # centred point for even "l"
  mindx = index + cbind(rep(l/2 -1, dim(index)[1]), rep(-(l/2-1), dim(index)[1]));
  # Ajusting data
  timedata = seq(from=data$V1[1], to=data$V1[n], length = n);
 
  # mean time (centred point) for our g-mode estimates
  timefreq = apply(mindx, 1, function(x) mean(timedata[c(x[1], x[2])]) );

  # g-modes
 
  if( !is.null(timeGmode) ){
   
    #timeGmode = data0[c(1,length(data0[,1])), 1];
   
    out = apply(as.matrix(timefreq), 1,
                function(x){
                  if((x >= timeGmode[1]) & (x <= timeGmode[2])){
                    return(TRUE);
                  }else{
                    return(FALSE);
                  }
                });
 
    #maxf     = maxf[out];
    timefreq = timefreq[out];
   
  }
 
  r$x = r$x[out];
  r$z = r$z[out, ];
 
  maxf = findGmodes(r, um=um, dm=dm);
  maxf = movf(maxf, movGmode, median); # smoothing g-mode estimates
  points(timefreq,maxf, col = 'green')
 
  if( class(mod) == "lm"){
   
    # prediction
    new   = data.frame(f = maxf);
    pred  = predict(mod, new, interval = "prediction"); # predictions
    title = "Frequentist analysis";
   
  }else if(class(postSamples) == "data.frame"){
   
    # prediction
    aux  = as.matrix(postSamples[,1:3]) %*% t(cbind(maxf, maxf^2, maxf^3));
   
    aux1 = apply(aux,2,median);
   
    sdb  = 1.96 * sqrt(maxf * median(postSamples[,4]));
   
    pred = cbind(apply(aux,2,median), aux1 - sdb, aux1 + sdb);
   
    title = "Bayesian analysis";
  }
 
  ### PLOT ###
 
  yaux = c(true_ratios, pred[,2:3]);
  plot(Obergaunlinguer_data$time, true_ratios, main = paste(title), xlab = "Time",
       ylab = "Ratio", ylim = c(min(yaux), max(yaux)), type = "n");
  arrows(timefreq, pred[,2], timefreq, pred[,3], code=3, angle=90,
         length=0.05, col="gray");
  points(Obergaunlinguer_data$time, true_ratios, col = "black");
  points(timefreq, pred[,1], col = "red", cex = pred[,1]/max(pred[,1])+ 0.3);
 
}
