source ("data_simulation.R")

dist=10    # kpc
fmin=10    # Hz

detectors=c("aLIGO", "ADV", "KAGRA", "CE1", "CE2", "ET_B", "ET_C", "ET_D")
signals=c("s11.2--LS220", "s15.0--LS220", "s15.0--SFHo", "s15.0--GShen", 
          "s20.0--LS220", "s20.0--SFHo", "s25.0--LS220", "s40.0--LS220")

#signals=c("s15.0--LS220", "s15.0--SFHo", "s15.0--GShen")
#det="aLIGO"
#signals="s11.2--LS220"

for (det in detectors){
print(det)
  for (sig in signals){
    snr=compute_SNR(sig, det, fmin, dist, pbOff=TRUE)
    print(sprintf('%s snr=%.3f', sig, snr))
  }
}

