#=============================================================================
# ARCHIVO CON FUNCIONES AUXILIARES
#=============================================================================

#=============================================================================
# FUNCIÓN valoresCercanos2()
#-----------------------------------------------------------------------------
# Aplica el Teorema del Valor Medio o de Bolzano para encontrar intervalos que
# contengan raíces de la función deseada. Además, evalúa la posibilidad de existencia
# de raíces sin cambio de signo.
# 
# Argumentos:
#   - f: función de la cual se quiere hallar raíces
#   - from, to, by: inicio, fin y salto para crear los intervalos, es decir, para crear la secuencia de valores de x a evaluar
#   - epsilon: valor pequeño, la función debe arrojar un valor menor que él para considerar que puede haber raiz sin cambio de signo
# Salida:
#   - Mensaje en pantalla con los intervalos en los que hay cambio de signo
#   - Matriz de dos columnas con el límite inferior y superior de cada uno de los intervalos hallados
#=============================================================================

valoresCercanos2 <- function(f, from, to, by, epsilon = 0.0001) {
    x = seq(from, to, by)
    y = f(x)
    intervalos = NULL
    for (i in 2:length(x)) {
        if (y[i - 1] * y[i] < 0) {
            # cat("Hay un cambio de signo entre ", x[i - 1], " y ", x[i], "\n")
            intervalos = rbind(intervalos, x[(i - 1) : i])
        }
        if (i < length(x)) {
            if ((y[i] - y[i - 1]) * (y[i + 1] - y[i]) < 0 & abs(y[i]) < epsilon) {
                # cat("Puede haber una raiz entre ", x[i - 1], " y ", x[i + 1], "\n")
                intervalos = rbind(intervalos, x[c(i - 1, i + 1)])
            }
        }
    }
    if (is.null(intervalos)) cat("No se encontraron cambios de signos\n")
    return(intervalos)
}

#=============================================================================
# FUNCIÓN newtonRaphson()
#-----------------------------------------------------------------------------
# Método de Newton-Raphson para resolver ecuaciones no lineales
# 
# Argumentos:
#   - f: función expresada en la forma x = f(x)
#   - x0: valor inicial
#   - errorRel: valor lógico, por default (TRUE) calcula error relativo, caso contrario, absoluto
#   - iterMax: número máximo de iteraciones permitido
#   - errorCotaSup: para el criterio de divergencia
#   - errorCotaInf: para el criterio de convergencia
# Salida:
#   - Mensaje en pantalla con el resultado
#   - Matriz con todos los pasos iterativos (nro de orden, xi, error)
#=============================================================================

newtonRaphson <- function(f, x0, errorRel = TRUE, iterMax = 1000, errorCotaSup = 1000, errorCotaInf = 1e-6) {
    
    # Derivada
    d = Deriv(f)
    
    # Valores iniciales
    iter <- 0
    error <- 1
    rtdo <- c(iter, x0, NA)
    
    # Iterar
    while (error >= errorCotaInf & error <= errorCotaSup & iter < iterMax) {
        xn <- x0 - f(x0) / d(x0)
        error <- ifelse(errorRel, abs((x0 - xn) / x0), abs(x0 - xn))
        x0 <- xn
        iter <- iter + 1
        rtdo <- rbind(rtdo, c(iter, x0, error))
    }
    
    # Mostrar resultado
    if (error < errorCotaInf) {
        # cat("Raiz aproximada:", xn, "\n")
        # cat("Iteraciones:", iter, "\n")
        # cat("Error:", error, "\n")
        rtdo <- as.data.frame(rtdo)
        colnames(rtdo) <- c("Iteracion", "x", "Error")
        rownames(rtdo) <- NULL
        return(rtdo)
    } else {
        cat("El método no converge")
    }
}

#=============================================================================
# FUNCIÓN analisis()
#-----------------------------------------------------------------------------
# Implementación del análisis de los puntos de equilibrio y la dinámica global
# del sistema como está propuesto en el libro de Li, sección 5.2.
#=============================================================================

analisis <- function(lambda, v1, v2, mu1, mu2, K, beta, sigma, nsteps) {
    
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
    # fy devuelve la coordenada en Y para los ptos de equilibrio con infeccion cronica (los x se obtienen aparte, son donde se cruzan f1 y f2)
    
    # Punto de equilibrio libre de infeccion: raiz positiva de f1
    #------------------------------------------------------------
    # Siempre lo puedo sacar
    # Necesito un valor inicial
    inicial = valoresCercanos2(f1, 0, 5000, 1)[1, 1]
    raiz = newtonRaphson(f1, inicial)
    x0 = raiz$x[nrow(raiz)]
    # El punto es:
    p0 = c(x0, 0)
    lineas = "----------------------------------------------------------\n"
    lineas = paste0(lineas, "Puntos de equilibrio (x, y): \n\n")
    lineas = paste0(lineas, "- Punto de equilibrio libre de infección: ", toString(round(p0, 2)), "\n") 
                    
    
    # Puntos de equilibrio con infeccion cronica, se dan donde se cruzan f1 y f2
    #------------------------------------------------------------
    # 1. Graficar funciones
    d =
        tibble(x = seq(0, nsteps, 1), f1 = f1(x), f2 = f2(x)) %>%
        gather(key = funcion, value = valor, f1, f2)
    g = ggplot(d, aes(x = x, y = valor, group = funcion, color = funcion)) +
        geom_line() +
        scale_y_continuous("f", limits = c(0, NA)) +
        scale_x_continuous("x") +
        theme_bw() +
        ggtitle("Existencia de puntos de equilibrios con infección crónica",
                subtitle = "La o las intersecciones entre ambas curvas indican la coordenada en x de los puntos de equilibrio con infección crónica.")

    # Puede tener ptos con inf cronica si se verifica 5.17
    cond5.17 = (sigma * K * beta < v2) & (mu2 < v2)
    # Si se cumple 5.17, seguro se cruzan f1 y f2 si se cumple 5.18:
    cond5.18 = (v1 + beta * K)^2 * (1 - mu2 / v2)^2 > (v1 - mu1)^2 + 4 * lambda * v1 / K
    # Se verifican?
    verifican = cond5.17 & cond5.18
    # cantidad = 0
    
    # Si se cumplen 5.17 y 5.18, se cruzan seguro, entonces procedo a buscar los puntos automaticamente
    if (verifican & sigma != 0.01) {
        
        lineas = paste0(lineas, "- Existen otros puntos de equilibrio con infección crónica:\n\n")
        
        # Para buscar las coordenadas de x, busco las raices de la resta f1 - f2 = 0
        # necesito valores cercanos
        iniciales = valoresCercanos2(f1menos2, 0, 5000, 1)
        raices = numeric(nrow(iniciales))
        for (i in 1:length(raices)) {
            raiz = newtonRaphson(f1menos2, iniciales[i, 1])
            raices[i] = raiz$x[length(raices)]
        }
        # Ahora busco los respectivas coordenadas de y con fy(raices)
        # y directamente combino
        Pestrella = cbind(raices, fy(raices)) # cada pto en una fila
        positivos = apply(Pestrella, 1, function(x) all(x >= 0))
        Pestrella = Pestrella[positivos, , drop = F]
        cantidad = 1 + nrow(Pestrella)
        for (i in 1:nrow(Pestrella)) {
            lineas = paste0(lineas, "       + Punto de equilibrio con infección crónica = ", toString(round(Pestrella[i, ], 2)), "\n")
        }
        
    } 
    # } else {
    #     lineas = paste0(lineas, "- No existen otros puntos de equilibrio.\n\n")
    #     cantidad = 1 # el de libre infeccion solo
    # }
    lineas = paste0(lineas, "----------------------------------------------------------\n")
    # Ademas, si se verifian 5.17 y 5.18, podemos ver los otros rtdos
    # Teorema 5.2.1
    sigma0 = (4 * lambda * beta * v2^2 - ((v1 - mu1) * v2 - (v1 + K * beta) * (v2 - mu2))^2) / (4 * lambda * beta * v2 * (v1 + K * beta))
    sigmac = (v2 / (K * beta)) - (v2 - mu2) / (beta * x0)
    lineas = paste0(lineas, "Fracción de células infectadas que sobreviven:\n\n")
    lineas = paste0(lineas, "sigma = ", round(sigma, 4), "\n")
    lineas = paste0(lineas, "Límite inferior sigma_0 = ", round(sigma0, 4), "\n")
    lineas = paste0(lineas, "Límite superior sigma_c = ", round(sigmac, 4), "\n")
    if (0 < sigma & sigma < sigma0) {
        # cantidad = 1
        lineas = paste0(lineas, "0 < sigma < sigma0 --> no hay puntos de equilibrio con infección crónica\n\n")
    } else if (sigma == sigma0) {
        # cantidad = 2
        lineas = paste0(lineas, "sigma = sigma0 --> hay un único punto de equilibrio con infección crónica\n\n")
    } else if (sigma0 < sigma & sigma < sigmac) {
        # cantidad = 3
        lineas = paste0(lineas, "sigma0 < sigma < sigmac --> hay dos puntos de equilibrio con infección crónica\n\n")
    } else {
        # cantidad = 2
        lineas = paste0(lineas, "sigma >= sigmac --> hay un único punto de equilibrio con infección crónica\n\n")
    }
    lineas = paste0(lineas, "----------------------------------------------------------\n")
    
    # Teorema 5.2.5, dinamica global
    lineas = paste0(lineas, "Número básico de reproducción R0\n\n")
    r0_sigma = (sigma * beta * x0 + v2 * (1 - x0 / K)) / mu2
    r0_sigma0 = (sigma0 * beta * x0 + v2 * (1 - x0 / K)) / mu2
    lineas = paste0(lineas, "R0 = ", round(r0_sigma, 4), "\n")
    lineas = paste0(lineas, "R0_sigma0 = ", round(r0_sigma0, 4), "\n")
    if (0 < r0_sigma & r0_sigma < r0_sigma0) {
        lineas = paste0(lineas, "0 < R0(sigma) < RO(sigma0) --> El punto de equilibrio libre de infección es globalmente estable asintóticamente\n")
    } else if (r0_sigma0 < r0_sigma & r0_sigma < 1) {
        lineas = paste0(lineas, "R0(sigma0) < R0 < 1 --> El sistema tiene dos puntos de atracción, un punto de equilibrio libre de infección y otro de infección crónica, más un punto de equilibrio inestable (punto de silla).\n")
    } else {
        lineas = paste0(lineas, "R0(sigma) > 1 --> el único punto de equilibrio con infección crónica es globalmente estable asintóticamente\n")
    }
    lineas = paste0(lineas, "----------------------------------------------------------\n")
    

    return(list(g = g, lineas = lineas, cantidad = cantidad))
    
}

