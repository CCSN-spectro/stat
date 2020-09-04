fits_data = read.table("inputs/A-A_fits_data_g2.dat", sep = ",");# data to generate model
CoCo  = read.table("inputs/CoCo_fits_data_g2.dat", sep = ",");# data to generate model
colnames(fits_data) = c("r", "f");
colnames(CoCo) = c("r", "f");

mu3 = "~ f + I(f^2) + I(f^3) - 1";
s2  = "~ f + I(f^2) - 1";

Xm  = model.matrix(eval(parse(text=eval(parse(text=mu3)))), fits_data);
Xs  = model.matrix(eval(parse(text=eval(parse(text=s2)))), fits_data);

fit = lmvar(fits_data$r, X_mu = Xm, X_sigma = Xs, intercept_mu = TRUE);

pdf("model.pdf")
par(mar = c(4,4,1,1))
colCC = "green";
colAA = "red";
plot(c(fits_data$f, CoCo$f), c(fits_data$r, CoCo$r), col = "gray", type = "n",
     xlab = "Frequency", ylab = "Ratio");
points(CoCo$f, CoCo$r, col = colCC, pch = 1);
points(fits_data$f, fits_data$r, col = colAA, pch = 1);
set.seed(1)
i = sample(length(CoCo$r), size = length(CoCo$r)/2)
points(CoCo$f[i], CoCo$r[i], col = colCC, pch = 1);
pred = predict(fit, X_mu = Xm, X_sigma = Xs, interval = "prediction"); # predictions
lwr = cbind(fits_data$f, pred[,3]); # lower band
lwr = lwr[order(lwr[,1]), ];
lines(lwr[,1], lwr[,2], col = "black", lwd = 2, lty = 2);
upr = cbind(fits_data$f, pred[,4]); # upper band
upr = upr[order(upr[,1]), ];
lines(upr[,1], upr[,2], col = "black", lwd = 2, lty = 2);
med = cbind(fits_data$f, pred[,1]); # middle point
med = med[order(med[,1]), ];
lines(med[,1], med[,2], col = "black", lwd = 2, lty = 1);

legend("topleft", legend=c("CoCoNuT", "AenusAlcar"), col=c(colCC, colAA), 
       pch=1, cex=0.8);
dev.off()
