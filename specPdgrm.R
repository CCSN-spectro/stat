source("functions.R")
#library(psplinePsd)
library(fields)

  specPdgrm = function(data, time = NULL, l, p, fs = 16384, actPlot = TRUE, logPow = TRUE, 
                     zoomFreq = c(0,1), main  = NULL){
  
  # fs: sampling frequency
  
  if((zoomFreq[1]<0) || (zoomFreq[2]>1) || (zoomFreq[1] > zoomFreq[2])){
    stop("zoomFreq must be a vector c(a,b) with values 0 <= a < b <= 1")
  }
  
  index = ints(n = length(data), l = l, p = p);
  
  R = NULL;
  
  for(i in 1:dim(index)[1]){
    
    auxdata = data[index[i,1]:index[i,2]];
    
    # Optimal rescaling to prevent numerical issues
    rescale = stats::sd(auxdata);
    auxdata = auxdata / rescale;  # Data now has standard deviation 1
    auxdata = auxdata - mean(auxdata)
    
    N       = length(auxdata)/2 + 1;
    pdgrm   = (abs(fft(auxdata))^2 / (2*pi*length(auxdata)))[1:N];
    aux     = pdgrm[-c(1,N)];
    pdgrm   = pdgrm * rescale^2;
    #aux     = log(pdgrm);
    R       = cbind(R, aux);
    
  }
  
  if(logPow == TRUE){R = log10(sqrt(R));}
  
  index    = round(dim(R)[1] * zoomFreq);
  index    = seq(index[1] + 1, index[2]);
  
  pc = 1-p/100;
  
  i <- 1:dim(R)[2];
  #i <- 1:(length(R)-1); # Rput form psplinePsd
  
  x <- round(l / fs * (i-1) * pc * 1000, 2) # number of intervals
  
  if(!is.null(time)){
    
    if(length(time)==2){
    
      x <- seq(from= time[1], to=time[2], length = length(x));
      
    }else{
      
      x = seq(from= time[1], to=time[length(time)], length = length(x));
      
    }
  }
  
  y <- seq(1,fs/2,length=dim(R)[1]);
  z <- t(R);
  #z <- z[, -c(1, dim(z)[2])]; #deleting first and last row in image
  #y <- y[-c(1, dim(z)[2])]; # adjustment
  
  y <- y[index]; # zoom
  z <- z[, index];
  rownames(z) =NULL;
  
  if(actPlot == TRUE){
    
    if( is.null(main) ){
      
      main = c("Spectrogram based on periodogram");
      
    }
    
    image.plot(x, y, z, xlab = "Time [ms]", ylab = "Frequency [Hz]", 
               main = paste(main)); # yaxt='n'
    
  }else{
    
    return(list(x = x, y = y, z = z));
    
  }
}
