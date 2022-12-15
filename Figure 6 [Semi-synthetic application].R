##################
#### Figure 6 ####
##################

excesserrors <- cbind(c(200,400,600,800,1000,1200,2000,4000), 
                      c(18.09,25.61,31.43,32.79,28.00,22.22,14.71,17.86), 
                      c(20.00,29.33,34.92,37.04,32.56,27.03,15.38,10.00),
                      c(30,37,37,34,24,20,13,14), 
                      c(36.89,32.63,25.26,21.51,17.78,21.18,33.33,47.46))

colnames(excesserrors) <- c("Number of observations", "Excess RMSE", "Excess absolute bias",
                            "Excess coverage rate", "Excess length")

par(mfrow = c(1, 1))

plot(x = excesserrors[,1], y = excesserrors[,2], type = "o", pch = 1, col = "red", main = "",
     xlim = c(200,4000), ylim = c(7, 52), xlab = "Number of observations", ylab = "Excess (in %)")
lines(x = excesserrors[,1], y = excesserrors[,3], type = "o", pch = 2, col = "orange")
lines(x = excesserrors[,1], y = excesserrors[,4], type = "o", pch = 5, col = "violet")
lines(x = excesserrors[,1], y = excesserrors[,5], type = "o", pch = 22, col = "brown")
legend(x = 3000, y = 36, box.lty = 0, cex = 0.75, pch = c(1,2, 5, 22), legend = c("RMSE", "Absolute bias", "Coverage", "Length"),
       lwd = 1.5, col = c("red", "orange", "violet", "brown"), xpd = TRUE)
