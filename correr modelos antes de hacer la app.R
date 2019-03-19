
####################################################################
#                        phaseR
###################################################################

library(phaseR)

# Elegir valores para los parámetros
params = c(lambda = 25, v1 = 0.038, v2 = 0.038, mu1 = 0.03, mu2 = 0.03, K = 1150, beta = 0.00103, sigma = 0.015)

# Especificacion del sistema
sistema <- function (t, y, parameters) {
    
    # Variables
    x = y[1]
    y = y[2]
    
    # Otros parametros
    lambda = parameters[1]
    v1 = parameters[2]
    v2 = parameters[3]
    mu1 = parameters[4]
    mu2 = parameters[5]
    K = parameters[6]
    beta = parameters[7]
    sigma = parameters[8]
    
    # Expresion de las derivadas
    dx = lambda + v1 * x * (1 - (x + y) / K) - mu1 * x - beta * x * y
    dy = sigma * beta * x * y + v2 * y * (1 - (x + y) / K) - mu2 * y
    
    # Devolver el sistema como una lista con un vector numerico
    return(list(c(dx, dy)))
}

# Crear el plano de fases
planoFases <- flowField(sistema,
                        xlim = c(0, 1000),
                        ylim = c(0, 1000),
                        parameters = params,
                        points = 50,
                        add = F)

# Isoclinas nulas
sistemaNullclines  <- nullclines(sistema, xlim = c(0, 1000), ylim = c(0, 1000), parameters = params)

# Buscar los puntos de equilibrio
findEquilibrium(sistema, parameters = params)
eq1 = c(145.4972, 164.6376)
eq2 = c(322.39858, 70.44923)
points(x = c(eq1[1], eq2[1]), y = c(eq1[2], eq2[2]), pch = 19)

# Clasificacion
# Determine the stability of the equilibrium points
stability(sistema, ystar = eq1, parameters = params)$classification
stability(sistema, ystar = eq2, parameters = params)$classification

# Puedo tratar de encontrarlos automaticamente?
# El tema es que hay que darle de donde arrancar... le voy a dar valores aleatorios
# hasta que los encuentre a todos, la cantidad que tiene que encontrar se la doy yo

nroEq <- 2
encontrados <- 0
nsteps = 800
puntos <- matrix(NA, nrow = nroEq, ncol = 2)
clases <- character(nroEq)
while(encontrados < nroEq) {
    y0 <-  runif(2, 1, nsteps)
    nuevo <- as.vector(findEquilibrium(sistema, parameters = params, y0 = y0)$ystar)
    nuevoClase <- findEquilibrium(sistema, parameters = params, y0 = y0)$classification
    if (encontrados == 0) {
        encontrados = encontrados + 1
        puntos[1, ] = nuevo
        clases[1] = nuevoClase
    } else {
        if (any(round(nuevo, 4) != round(puntos[encontrados, ], 4))) {
            encontrados = encontrados + 1
            puntos[encontrados, ] = nuevo
            clases[encontrados] = nuevoClase
        }
    }
}
rtdo <- as.data.frame(puntos)
colnames(rtdo) <- c("Coord.X.NoInf", "Coord.Y.Inf")
rtdo$Clasificacion <- clases
rtdo

####################################################################
#                        EpiModel
###################################################################

library(EpiModel)

# Definir el modelo
inHost <- function(t, t0, parms) {
    with(as.list(c(t0, parms)), {
        
        # Tamano poblacion
        num = x + y
        
        # Calculos dinamicos (opcionales)
        # prevalencia = y / num
        
        # Derivadas
        dx = lambda + v1 * x * (1 - (x + y) / K) - mu1 * x - beta * x * y
        dy = sigma * beta * x * y + v2 * y * (1 - (x + y) / K) - mu2 * y

        # Compartments and flows are part of the derivative vector
        # Other calculations to be output are outside the vector, but within the containing list
        list(c(dx, dy, 
               ax.flow = lambda + v1 * x * (1 - (x + y) / K),   # nuevas celulas no infectadas (arrival)
               ay.flow = v2 * y * (1 - (x + y) / K),            # nuevas celulas infectadas, por mitosis (arrival)
               si.flow = sigma * beta * x * y,                  # tasa de S a I, x es S, y es I
               ds.flow = mu1 * x,                               # removal rate para no infectadas
               di.flow = mu2 * x                                # removal rate para infectadas
               ),
             num = num,
             lambda = lambda,
             v1 = v1, 
             v2 = v2,
             K = K,
             beta = beta,
             sigma = sigma,
             # i.prev = i.num / num,
             prevalencia = y / num)
    })
}


# Parametrizacion
param <- param.dcm(lambda = 25, v1 = 0.038, v2 = 0.038, mu1 = 0.03, mu2 = 0.03, K = 1150, beta = 0.00103, sigma = 0.015)
init <- init.dcm(x = 1000, y = 50, ax.flow = 0, ay.flow = 0, si.flow = 0, ds.flow = 0, di.flow = 0)
control <- control.dcm(nsteps = 800, new.mod = inHost)

# Correr el modelo
mod <- dcm(param, init, control)
mod

# Graficos que salen automaticos
plot(mod, y = "x")
plot(mod, y = "y", add = T)
plot(mod, y = "ax.flow")
plot(mod, y = "ay.flow")
plot(mod, y = "si.flow")
plot(mod, y = "ds.flow")
plot(mod, y = "di.flow")
plot(mod, y = "num")
plot(mod, y = "lambda")
plot(mod, y = "beta")
plot(mod, y = "prevalencia")


# Mis propios graficos
# ----------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

res <- as.data.frame(mod$epi)
colnames(res) <- names(mod$epi)
res$Tiempo <- 1:nrow(res)
names(res)

res.x.y <- 
    res %>% 
    gather(key = Tipo, value = Numero, x, y) %>% 
    select(Tiempo, Tipo, Numero) %>% 
    mutate(Tipo = factor(Tipo, levels = c("x", "y"), labels = c("No infectadas", "Infectadas")))

res.x.y.prev <- 
    res %>% 
    mutate(prevalencia.x = 1 - prevalencia) %>% 
    gather(key = Tipo, value = Prevalencia, prevalencia, prevalencia.x) %>% 
    select(Tiempo, Tipo, Prevalencia) %>% 
    mutate(Tipo = factor(Tipo, levels = c("prevalencia.x", "prevalencia"), labels = c("No infectadas", "Infectadas")))

res.a.d <- 
    res %>% 
    gather(key = Tipo, value = Tasa, ax.flow, ay.flow, ds.flow, di.flow) %>% 
    mutate(Tipo = str_remove(Tipo, ".flow"), 
           Tipo = str_replace(Tipo, "s", "x"),
           Tipo = str_replace(Tipo, "i", "y")) %>% 
    separate(Tipo, into = c("Direccion", "Tipo"), sep = 1) %>% 
    mutate(Direccion = factor(Direccion, levels = c("a", "d"), labels = c("Crecimiento (arrival)", "Eliminación (departure)")),
           Tipo = factor(Tipo, levels = c("x", "y"), labels = c("Células no infectadas", "Células infectadas")))

# Num de celulas infectadas y no infectadas
ggplot(res.x.y, aes(x = Tiempo, y = Numero, group = Tipo, color = Tipo)) + 
    geom_line(lwd = 1.5) + 
    theme_bw() + 
    theme(legend.position = "bottom")

# Prevalencia
ggplot(res.x.y.prev, aes(x = Tiempo, y = Prevalencia, group = Tipo, color = Tipo)) + 
    geom_line(lwd = 1.5) + 
    theme_bw() + 
    theme(legend.position = "bottom")

# Incidencia
ggplot(res, aes(x = Tiempo, y = si.flow)) + 
    geom_line(lwd = 1.5) + 
    scale_y_continuous("Incidencia") +
    theme_bw()

# Tasa de aparición y eliminación de células
ggplot(res.a.d, aes(x = Tiempo, y = Tasa)) + 
    geom_line() + 
    facet_grid(Direccion ~ Tipo, scales = "free") + 
    theme_bw()

# Tamaño de la población
ggplot(res, aes(x = Tiempo, y  = num)) +
    geom_line() +
    scale_y_continuous("Tamaño de la población", limits = c(0, NA)) + 
    theme_bw()

# Resumen a un tiempo t
#----------------------

summary(mod, at = 150) # no anda con modelos nuevos definidos por uno

# resumen de estadsiticas en un tiempo dado
at = 152
summDig = 3
res2 = round(res[at, 1:8], summDig)

ver2 = paste0(
    "----------------------------------------------------
    Estadísticas del modelo para Tiempo = ", at, "\n",
    "----------------------------------------------------
    Población susceptible: ", res2$x,
    "\n    Población infectada: ", res2$y,
    "\n    Total: ", res2$num,
    "\n    Suc -> Inf: ", res2$si.flow,
    "\n    Nuevas células no infectadas: ", res2$ax.flow,
    "\n    Nuevas células infectadas: ", res2$ay.flow,
    "\n    Células no infectadas muertas: ", res2$ds.flow,
    "\n    Células no infectadas muertas: ", res2$di.flow
)

writeLines(ver2)

#----------------------------------
# ANALISIS CON LA TEORIA DEL LIBRO
#----------------------------------

library(Deriv)
library(ggplot2)   

lambda = 25
v1 = 0.038
v2 = 0.038
mu1 = 0.03
mu2 = 0.03
K = 1150
beta = 0.00103
sigma = 0.015

# Funciones que generan una funcion que solo depende de x, pero con los parametros dados
f1.generador <- function(lambda, v1, v2, mu1, mu2, K, beta, sigma) {
    function(x) lambda + (v1 - mu1) * x - v1 / K * x^2
}
f2.generador <- function(lambda, v1, v2, mu1, mu2, K, beta, sigma) {
    function(x) (1 / v2) * (v1 / K + beta) * ( (sigma * K * beta - v2) * x + K * (v2 - mu2) ) * x
}
f1menos2.gen <- function(lambda, v1, v2, mu1, mu2, K, beta, sigma) {
    function(x) {
       (lambda + (v1 - mu1) * x - v1 / K * x^2) - ((1 / v2) * (v1 / K + beta) * ( (sigma * K * beta - v2) * x + K * (v2 - mu2) ) * x)
    }
}
fy.generador <- function(lambda, v1, v2, mu1, mu2, K, beta, sigma) {
    function(x) (K / v2) * ( (sigma * beta - v2 / K) * x + (v2 - mu2) ) 
}
f1 = f1.generador(lambda, v1, v2, mu1, mu2, K, beta, sigma)
f2 = f2.generador(lambda, v1, v2, mu1, mu2, K, beta, sigma)
f1menos2 = f1menos2.gen(lambda, v1, v2, mu1, mu2, K, beta, sigma)
fy = fy.generador(lambda, v1, v2, mu1, mu2, K, beta, sigma)

# Punto de equilibrio libre de infeccion: raiz positiva de f1
#------------------------------------------------------------
# Siempre lo puedo sacar
# Necesito un valor inicial
inicial = valoresCercanos2(f1, 0, 5000, 1)[1, 1]
raiz = newtonRaphson(f1, inicial)
x0 = raiz$x[nrow(raiz)]
# El punto es:
c(x0, 0)

# Puntos de equilibrio con infeccion cronica, se dan donde se cruzan f1 y f2
#------------------------------------------------------------
# 1. Graficar funciones
ggplot(NULL, aes(x = c(0, nsteps))) + 
    stat_function(fun = f1, col = "blue") +
    stat_function(fun = f2, col = "red") +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_continuous("x") +
    theme_bw() + 
    ggtitle("Gráfico de las funciones f1 (azúl) y f2 (rojo)")

# Puede tener ptos con inf cronica si se verifica 5.17
(sigma * K * beta < v2) & (mu2 < v2)

# Si se cumple 5.17, seguro se cruzan f1 y f2 si se cumple 5.18:
(v1 + beta * K)^2 * (1 - mu2 / v2)^2 > (v1 - mu1)^2 + 4 * lambda * v1 / K

# Si se cumplen 5.17 y 5.18, se cruzan seguro, entonces procedo a buscar los puntos automaticamente
# Para buscarlas, busco las raices de la resta f1 - f2 = 0
# necesito valores cercanos
iniciales = valoresCercanos2(f1menos2, 0, 5000, 1)
raices = numeric(nrow(iniciales))
for (i in 1:length(raices)) {
    raiz = newtonRaphson(f1menos2, iniciales[i, 1])
    raices[i] = raiz$x[length(raices)]
}
raices
# Ahora busco los respectivas coordenadas de y
fy(raices)
cbind(raices, fy(raices)) # cada pto en una fila

# Teorema 5.2.1, si se verifican 5.17 y 5.18
sigma0 = (4 * lambda * beta * v2^2 - ((v1 - mu1) * v2 - (v1 + K * beta) * (v2 - mu2))^2) / (4 * lambda * beta * v2 * (v1 + K * beta))
sigmac = (v2 / (K * beta)) - (v2 - mu2) / (beta * x0)

if (0 < sigma & sigma < sigma0) {
    print("0 < sigma < sigma0 --> no hay puntos de equilibrio con infección crónica")
} else if (sigma == sigma0) {
    print("sigma = sigma0 --> hay un único punto de equilibrio con infección crónica")
} else if (sigma0 < sigma & sigma < sigmac) {
    print("sigma0 < sigma < sigmac --> hay dos puntos de equilibrio con infección crónica")
} else {
    print("sigma >= sigmac --> hay un único punto de equilibrio con infección crónica")
}

# Teorema 5.2.5, si se verifican 5.17 y 5.18, dinamica global
r0_sigma = (sigma * beta * x0 + v2 * (1 - x0 / K)) / mu2
r0_sigma0 = (sigma0 * beta * x0 + v2 * (1 - x0 / K)) / mu2

if (0 < r0_sigma & r0_sigma < r0_sigma0) {
    print("0 < R0(sigma) < RO(sigma0) --> El punto de equilibrio libre de infección es globalmente estable asintóticamente")
} else if (r0_sigma0 < r0_sigma & r0_sigma < 1) {
    print("R0(sigma0) < R0(sigma) < 1 --> el sistema tiene dos puntos de atracción, el punto de equilibrio libre de infección y uno de infección crónica, y otro punto infección crónica que es un punto de silla.")
} else {
    print("R0(sigma) > 1 --> el único punto de equilibrio con infección crónica es globalmente estable asintóticamente")
}

lambda = 25
v1 = 0.038
v2 = 0.038
mu1 = 0.03
mu2 = 0.03
K = 1150
beta = 0.00103
sigma = 0.015

