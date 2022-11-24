
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

semisynthetic.rmse <- cbind(c(200,400,600,800,1000,1200,2000,4000),
                            c(0.124,0.106,0.099,0.086,0.082,0.078,0.064,0.051),
                            c(0.122,0.105,0.100,0.086,0.083,0.080,0.068,0.062),
                            c(0.112,0.091,0.091,0.078,0.071,0.069,0.056,0.048))

semisynthetic.bias <- cbind(c(200,400,600,800,1000,1200,2000,4000),
                            c(0.106,0.091,0.086,0.072,0.069,0.065,0.054,0.041),
                            c(0.103,0.089,0.084,0.070,0.068,0.066,0.055,0.050),
                            c(0.095,0.076,0.075,0.064,0.057,0.056,0.045,0.038))

semisynthetic.coverage <- cbind(c(200,400,600,800,1000,1200,2000,4000),
                            c(0.739,0.738,0.704,0.765,0.756,0.762,0.804,0.845),
                            c(0.843,0.870,0.856,0.901,0.904,0.905,0.934,0.940),
                            c(0.920,0.904,0.876,0.919,0.922,0.909,0.931,0.930))

semisynthetic.length <- cbind(c(200,400,600,800,1000,1200,2000,4000),
                            c(0.307,0.254,0.233,0.221,0.209,0.202,0.183,0.158),
                            c(0.380,0.341,0.327,0.317,0.309,0.301,0.286,0.264),
                            c(0.451,0.358,0.319,0.314,0.280,0.262,0.227,0.189))

par(mfrow = c(1, 1))

plot(x = semisynthetic.rmse[,1], y = semisynthetic.rmse[,2], type = "o", pch = 1, col = "red", main = "",
     xlim = c(100,4200), ylim = c(0.040, 0.130), xlab = "Dimension", ylab = "RMSE")
lines(x = semisynthetic.rmse[,1], y = semisynthetic.rmse[,3], type = "o", pch = 2, col = "purple")
lines(x = semisynthetic.rmse[,1], y = semisynthetic.rmse[,4], type = "o", pch = 22, col = "orange")
legend(x = 2900, y = 0.122, box.lty = 0, cex = 0.75, pch = c(1,2,22), legend = c("GRF", "CR.GRF","LLCF"),
       lwd = 1.5, col = c("red", "purple","orange"), xpd = TRUE)

plot(x = semisynthetic.bias[,1], y = semisynthetic.bias[,2], type = "o", pch = 1, col = "red", main = "",
     xlim = c(100,4200), ylim = c(0.030, 0.110), xlab = "Dimension", ylab = "Absolute bias")
lines(x = semisynthetic.bias[,1], y = semisynthetic.bias[,3], type = "o", pch = 2, col = "purple")
lines(x = semisynthetic.bias[,1], y = semisynthetic.bias[,4], type = "o", pch = 22, col = "orange")
legend(x = 2900, y = 0.105, box.lty = 0, cex = 0.75, pch = c(1,2,22), legend = c("GRF", "CR.GRF","LLCF"),
       lwd = 1.5, col = c("red", "purple","orange"), xpd = TRUE)

plot(x = semisynthetic.coverage[,1], y = semisynthetic.coverage[,2], type = "o", pch = 1, col = "red", main = "",
     xlim = c(100,4200), ylim = c(0.67, 0.965), xlab = "Dimension", ylab = "Coverage rate")
lines(x = semisynthetic.coverage[,1], y = semisynthetic.coverage[,3], type = "o", pch = 2, col = "purple")
lines(x = semisynthetic.coverage[,1], y = semisynthetic.coverage[,4], type = "o", pch = 22, col = "orange")
legend(x = 2900, y = 0.78, box.lty = 0, cex = 0.75, pch = c(1,2,22), legend = c("GRF", "CR.GRF","LLCF"),
       lwd = 1.5, col = c("red", "purple","orange"), xpd = TRUE)

plot(x = semisynthetic.length[,1], y = semisynthetic.length[,2], type = "o", pch = 1, col = "red", main = "",
     xlim = c(100,4200), ylim = c(0.14, 0.47), xlab = "Dimension", ylab = "Interval length")
lines(x = semisynthetic.length[,1], y = semisynthetic.length[,3], type = "o", pch = 2, col = "purple")
lines(x = semisynthetic.length[,1], y = semisynthetic.length[,4], type = "o", pch = 22, col = "orange")
legend(x = 2900, y = 0.45, box.lty = 0, cex = 0.75, pch = c(1,2,22), legend = c("GRF", "CR.GRF","LLCF"),
       lwd = 1.5, col = c("red", "purple","orange"), xpd = TRUE)




