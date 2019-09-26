puntosTodos = NULL
for (sig in c(seq(0.001, 0.015, 0.001), seq(0.02, 0.29, 0.02))) {
     # sig = 0.001
    # print(sig)
    analisisRtdo = analisis(lambda = lambda, v1 = v1, v2 = v2, mu1 = mu1, mu2 = mu2, K = K, beta = beta, sigma = sig, nsteps = nsteps)
    puntos = matrix(analisisRtdo$p0, nrow = 1)
    if (analisisRtdo$cantidad > 1) puntos = rbind(puntos, analisisRtdo$Pestrella)
    puntos = as.data.frame(puntos)
    colnames(puntos) = NULL
    puntos$sigma = sig
    puntosTodos = rbind(puntosTodos, puntos)
}
colnames(puntosTodos) = c("x", "y", "sigma")

library(ggplot2)
ggplot(puntosTodos, aes(x = sigma, y = y, group = sigma)) + 
    geom_line()
