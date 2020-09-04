source ("data_simulation.R")

dist=10    # kpc
fmin=10    # Hz

detectors=c("aLIGO", "ADV", "KAGRA", "CE1", "CE2", "ET_B", "ET_C", "ET_D")
signals=c("s20.0--SFHo", "s11.2--LS220", "s15.0--GShen", "s15.0--LS220", 
          "s15.0--SFHo", "s20.0--LS220", "s25.0--LS220", "s40.0--LS220")

for (det in detectors){
  print(det)
  for (sig in signals){
     compute_SNR(sig, det, fmin, dist)
  }
}

