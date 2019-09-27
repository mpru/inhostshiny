puntosTodos = NULL
# c(seq(0.001, 0.015, 0.001), seq(0.02, 0.29, 0.02))
for (sig in seq(0.001, 0.15, 0.001)) {
     # sig = 0.001
    # print(sig)
    analisisRtdo = analisis(lambda = lambda, v1 = v1, v2 = v2, mu1 = mu1, mu2 = mu2, K = K, beta = beta, sigma = sig, nsteps = nsteps)
    puntos = matrix(analisisRtdo$p0, nrow = 1)
    if (analisisRtdo$cantidad > 1) puntos = rbind(puntos, analisisRtdo$Pestrella)
    puntos = as.data.frame(puntos)
    colnames(puntos) = NULL
    puntos$sigma = sig
    if (nrow(puntos) == 1) {
        puntos$orden = 1
    } else if (nrow(puntos) == 2) {
        puntos$orden = c(1, 3)[order(puntos[, 2])]
    } else {
        puntos$orden = (1:3)[order(puntos[, 2])]
    }
    puntosTodos = rbind(puntosTodos, puntos)
    colnames(puntos) = NULL
}
colnames(puntosTodos) = c("x", "y", "sigma", "orden")

library(ggplot2)
ggplot(puntosTodos, aes(x = sigma, y = y, group = orden)) + 
    geom_line()


sigma0 = 0.25
sigmac = 0.5
sigma = 0.1

vertices = data.frame(
    xmin = c(0, sigma0, sigmac),
    xmax = c(sigma0, sigmac, 1),
    ymin = rep(-1, 3),
    ymax = rep(1, 3),
    color = factor(1:3)
)

ggplot(vertices, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = color)) +
    geom_rect() +
    scale_y_continuous(limits = c(-3, 3)) +
    geom_vline(xintercept = sigma, lty = 2) +
    scale_x_continuous(expression(sigma), expand = c(0, 0), limits = c(0, 1)) +
    scale_fill_manual("Estabilidad", values = c("darkgreen", "yellow", "red"), 
                       labels = c("Punto de equilibrio (PE) libre de infección estable. No hay PE con infección crónica.",
                                                 "PE libre de infección estable y dos PE con infección crónica, uno estable y otro inestable.",
                                                 "PE libre de infección inestable y un único punto de equilibrio con infección crónica"),
                      guide = guide_legend(nrow = 3)) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
