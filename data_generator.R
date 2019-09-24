library ("signal")
library ("seewave")

########################################################################
data_generator = function (fs, duration, wvf=NULL, ampl=1, fcut=15, actPlot=TRUE){
########################################################################
  # fs: sampling frequency
  # duration: duration (in second) of the output time serie
  # wvf: dataframe (V1=time, V2=h(t)) that contains the signal waveform sampled at fs
  # ampl: multiplication factor to change the source distance
  # fcut: lower frequency cutoff of the high pass filter. If NULL no filter.
  
  # check that duration is an integer > 0

  if ((duration%%1)>0){
    print("duration must be a positive integer")
    return()
  }
 
  n=duration*fs 
  n_0=n
  wvf_opt=0

  # Prepration in case a signal is added  
  if (length(wvf)!=1) {
    wvf_opt=1
    wvf_size=length(wvf.df$V2)
    if (wvf_size>n){
      print("The waveform length is larger than the requested data vector size. Please check")
      return(NULL)
    }
    if (2*wvf_size>n){
      #print("vector's duration is doubled")
      n=2*duration*fs
    }
    if (wvf_size%%2>0){
      print("The waveform length should be an even number. Please check")
      return(NULL)
    }
  }
  
  # Noise generation
  freq2 = fs*fftfreq(n)   # two-sided frequency vector
  freq1=freq2[1:int(n/2)]        # one-sided frequency vector
  psd2=aLIGO_PSD_new(freq2,2)    # two-sided PSD          

  T = seq(1, n, by = 1)  
  X = rnorm(n, mean=0, sd=1);          # Gaussian white noise
  XX = sqrt(fs) * fft(X) / (n*sqrt(n)); # FFT computing and normalization
  XXX = XX*sqrt(psd2);                  # Coloring
  Y = fft(XXX, inverse = TRUE);         # FFT inverse
  Y = Re(Y)*sqrt(n);                    # noise in time domain

  if (actPlot==TRUE){  

  }

  # Signal addition (centered at the middle of the data vector to avoid filtering leakage
  # at the beggining and end).
  if (wvf_opt==1){
    ind1=(n-wvf_size)/2
    
    for (i in 1:wvf_size){
      Y[ind1+i]=Y[ind1+i]+ampl*wvf.df$V2[i]
    }
  }

  # filter noise + signal
  # filtfilt : zero phase filter (forward& backward)
  myfilter=butter(n=4,W=fcut/(fs/2),type="high")
  YY=filtfilt(filt=myfilter,x=Y)          

  # select the orginal data size
  Tf = wvf.df$V1
  Yf = seq(1, n_0, by = 1)
  YYf = seq(1, n_0, by = 1)
  
  for (i in 1:n_0){
    Yf[i]=Y[ind1+i]
    YYf[i]=YY[ind1+i]
  }
    
  if (actPlot==TRUE){  
    ts.plot(Y); # noise only
    points(T,YY,col="red",type="l",pch=1,panel.first = grid()); # (noise + signal) filtered
    T_wvf=seq(ind1,ind1+wvf_size-1,by=1)
    points(T_wvf,(wvf.df$V2)*ampl,col="green",type="l",pch=3);  # signal only
#    points(Tf,YYf,col="green",type="l",pch=3);  # signal only
    leg = c("noise+signal HP filtered", "signal only")
    col = c("red","green")
    legend(x=0,y=-5*10^-21,legend=leg,cex=.8,col=col,pch=c(1,3))
  }

  return(list(t=Tf,x=Yf,y=YYf))
  }

########################################################################
fftfreq = function(n, d = 1){
########################################################################
  # surogate for the numpy fft.fftfreq function that generates the two sided 
  # frequency vector. Defaults d=1 means sampling frequency is 1. 
  # https://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.fftfreq.html
  #
  # n: samples number
  # d: sample spacing (inverse of the sampling rate). Defaults to 1
  
  if(n%%2 == 0){# n is even
    out = c(seq(0, n/2-1, by = 1), seq(-n/2, -1, by=1)) / (d*n);
  }else{ # n is odd
    out = c(seq(0, (n-1)/2, by = 1), seq(-(n-1)/2, -1, by=1)) / (d*n);
  }
  
  return(out);
}

########################################################################
int = function(n){
########################################################################
  # https://stackoverflow.com/questions/31036098/what-is-the-difference-between-int-and-floor-in-python-3
  if(n < 0 ){
    return(-floor(abs(n)))
  }else{
    return(floor(n))
  }
}

########################################################################
aLIGO_PSD = function(f,type){
########################################################################
  # Original aLIGO PSD function used by Patricio
  cutoff = -109.35 + log(2e10);
  fn     = length(f);
  logpsd = rep(0, fn);
  if(f[1]==0){
    f[1]=f[2]
  }
  if(type == 1){
    for(i in 1:fn){
      x  = f[i]/215;
      x2 = x*x;
      logpsd[i] = log(1e-49) + log(x^(-4.14) -5/x2 + 111*(1-x2+0.5*x2*x2)/(1.+0.5*x2));
      
      if(logpsd[i]>cutoff){
        logpsd[i]=cutoff;
      }
    }
    output=exp(logpsd);
    
  }else{
    for(i in 1:(int(fn/2)+1)){
      x  = abs(f[i]/215);
      x2 = x*x;
      logpsd[i] = log(1e-49) + log(x^(-4.14) -5./x2 + 111.*(1-x2+0.5*x2*x2)/(1.+0.5*x2));
      
      if(logpsd[i]>cutoff){
        logpsd[i]=cutoff;
      }
      if(i>0){
        logpsd[fn-i]=logpsd[i];
      }
    }
    output=exp(logpsd)/2;            # Two sided PSD
  }
  return (output)
}

########################################################################
aLIGO_PSD_new = function(f,type){
########################################################################
  # aLIGO sensitivity curve: fit the data point from https://dcc.ligo.org/LIGO-T1800044/public
  # Type=1 --> one-sided PSD.
  # Type=2 --> two-sided PSD. 
  
  S1 = 5.0e-26;
  S2 = 1.0e-40;
  S3 = 1.4e-46;
  S4 = 2.7e-51;
  fcut = 10;
  cutoff=1e-42;
  fn     = length(f);
  output = rep(0, fn); #np.zeros(len(f))
  
  # to avoid issue with f=0
  if(f[1]==0){
    f[1]=f[2];
  }
  if(type == 1){
    for(i in 1:fn){
      x = abs(f[i]);
      output[i] =  S1/(x^20) + S2/(x^4.05) + S3/(x^.5) + S4*((x/fcut)^2);
      if(output[i]>cutoff){
        output[i]=cutoff
      }
    }
  }else{
    for(i in 1:(int(fn/2)+1)){ # range(int(len(f)/2)+1)
      x = abs(f[i]);
      output[i] = S1/(x^20) + S2/(x^4.05) + S3/(x^.5) + S4*((x/fcut)^2);
      if(output[i]>cutoff){
        output[i]=cutoff
      }
      if(i>0){
        output[fn-i]=output[i]
      }
    }  
    output=output/2;          # Two sided PSD
  }
  return(output);
}

    