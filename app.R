#--------------------------------------
# Paquetes
#--------------------------------------

library(shiny)
library(EpiModel)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(phaseR)
library(Deriv)

#--------------------------------------
# Cargar funciones auxiliares
#--------------------------------------

source("funcionesAuxiliares.R")

#--------------------------------------
# Definir el modelo con EpiModel
#--------------------------------------

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

#--------------------------------------
# Definir el modelo con EpiModel
#--------------------------------------

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

#--------------------------------------
# UI (user interface)
#--------------------------------------

ui <- fluidPage(
    
    titlePanel("Modelo In-host: Modelo Compartimental para Infecciones Virales en Poblaciones de Células T CD4+"),
    withMathJax(),
    sidebarLayout(
        sidebarPanel(
            
            h3("Marcos Prunello", style = "margin-top: 0px"),
            h3("Asignatura: Modelos Matemáticos Continuos", style = "margin-top: 0px"),
            helpText("Hacer click en Correr Modelo después de cambiar",
                     "los parámetros."),
            actionButton(inputId = "runMod", "Correr modelo"),
            
            h4("Conditiones iniciales", style = "margin-top: 25px"),
            numericInput(inputId = "x",
                         label = "Población susceptible",
                         value = 1000,
                         min = 0),
            numericInput(inputId = "y",
                         label = "Población infectada",
                         value = 50,
                         min = 0),
            h4("Tiempo", style = "margin-top: 25px"),
            numericInput(inputId = "nsteps",
                         label = "Seguimiento (días)",
                         value = 800,
                         min = 0),
            h4("Parámetros", style = "margin-top: 25px"),
            numericInput(inputId = "lambda",
                         label = HTML("Tasa constante de generación de nuevas células T CD4+, no infectadas (&lambda;)"),
                         min = 0,
                         value = 25),
            numericInput(inputId = "v1",
                         label = "Tasa constante de crecimiento de las células T CD4+ no infectadas (v1)",
                         min = 0,
                         value = 0.038),
            numericInput(inputId = "v2",
                         label = "Tasa constante de crecimiento de las células T CD4+ infectadas (v2)",
                         min = 0,
                         value = 0.038),
            numericInput(inputId = "mu1",
                         label = HTML("Tasa de eliminación de las células T CD4+ no infectadas (&mu;1)"),
                         min = 0,
                         value = 0.03),
            numericInput(inputId = "mu2",
                         label = HTML("Tasa de eliminación de las células T CD4+ infectadas (&mu;2)"),
                         min = 0,
                         value = 0.03),
            numericInput(inputId = "K",
                         label = "Nivel de células T CD4+ en el que se detiene la división mitótica (K)",
                         min = 0,
                         value = 1150),
            numericInput(inputId = "beta",
                         label = HTML("Coeficiente de transmisión (&beta;)"),
                         min = 0,
                         value = 0.00103),
            numericInput(inputId = "sigma",
                         label = HTML("Fracción de células infectadas por contacto directo que resisten respuestas inmunes (&sigma;)"),
                         min = 0,
                         max = 1,
                         value = 0.015)
        ), # End sidebarPanel
        
        mainPanel(
            tabsetPanel(
                tabPanel("Sobre el modelo",
                         includeHTML("SobreElModelo.html")
                ), # Fin panel Sobre el modelo
                
                tabPanel("Gráficos",
                         h4("Resultados del modelo (uso del paquete EpiModel)"),
                         plotOutput(outputId = "MainPlot"),
                         br(),
                         wellPanel(
                             h4("Opciones gráficas"),
                             fluidRow(
                                 column(5,
                                        selectInput(inputId = "compsel",
                                                    label = strong("Seleccionar gráfico"),
                                                    choices = c("Población de células infectada y no infectada",
                                                                "Prevalencia",
                                                                "Incidencia",
                                                                "Tasa de aparición y eliminación de células",
                                                                "Tamaño de la población"))))
    
                         ) # End wellPanel
                ), # End tabPanel
                
                tabPanel("Resumen",
                         h4("Resumen del modelo en un tiempo específico (uso del paquete EpiModel)"),
                         helpText("Seleccionar el tiempo de interés y la cantidad
                                    de decimales para el redondeo de resultados"),
                         fluidRow(
                             column(5,
                                    numericInput(inputId = "summTs",
                                                 label = strong("Tiempo"),
                                                 value = 1,
                                                 min = 1)),
                             column(5,
                                    numericInput(inputId = "summDig",
                                                 label = strong("Cantidad de decimales"),
                                                 value = 3,
                                                 min = 0,
                                                 max = 8))),
                         fluidRow(verbatimTextOutput(outputId = "outSummary1")),
                         fluidRow(verbatimTextOutput(outputId = "outSummary2"))
                        ), # end tabPanel Summary
                
                tabPanel("Datos",
                         h4("Datos del modelo (uso del paquete EpiModel)"),
                         dataTableOutput("outData"),
                         fluidRow(
                             column(4,
                                    numericInput(inputId = "tabdig",
                                                 label = "Cantidad de decimales",
                                                 min = 0,
                                                 value = 2))),
                         downloadButton(outputId = "dlData",
                                        label = "Descargar datos")
                ), # end tabPanel Data

                tabPanel("Análisis",
                         h4("Resultados análiticos según las relaciones entre los parámetros descriptas en el libro de Li"),
                         verbatimTextOutput(outputId = "analisis"),
                         plotOutput(outputId = "f1f2", width = "80%"),
                         h4("Resultados obtenidos a través de la automatización del análisis del modelo con el paquete phaseR."),
                         helpText("Advertencia: se ha automatizado la búsqueda de puntos de
                                  equilibrio con arranques aleatorios, si no se encuentran en 
                                  los primeros 500 intentos abandona la búsqueda. Hacer click
                                  nuevamente en Correr Modelo para darle otra oportunidad. 
                                  Esto debe ser mejorado."),
                         # helpText("Advertencia: se ha automatizado la búsqueda de puntos de
                         #          equilibrio con arranques aleatorios, pero falta implementar
                         #          un control para informar sólo puntos que sean de interés en
                         #          el problema, es decir, dentro del rango de tiempo estudiado.
                         #          Por esta razón podrían aparecer puntos con coordenadas negativas,
                         #          por ejemplo."),
                         numericInput(inputId = "cifrasEq",
                                      label = "Cantidad de decimales",
                                      min = 0,
                                      max = 8,
                                      value = 2),
                         dataTableOutput("outDataEq"),
                         plotOutput(outputId = "planoFases", width = "80%")
                ) # Fin panel Plano de fases
            ) # end tabsetPanel
    ) # end mainPanel
) # end sidebarLayout
)

#--------------------------------------
# SERVER
#--------------------------------------



# Define server logic required to draw a histogram
server <- function(input, output) {
    
    ## Main reactive functions
    
    # Parametrización
    param <- reactive({
        param.dcm(lambda = input$lambda, 
                  v1 = input$v1, 
                  v2 = input$v2, 
                  mu1 = input$mu1, 
                  mu2 = input$mu2, 
                  K = input$K, 
                  beta = input$beta, 
                  sigma = input$sigma)
        })
    init <- reactive({
        init.dcm(x = input$x, y = input$y, ax.flow = 0, ay.flow = 0, si.flow = 0, ds.flow = 0, di.flow = 0)
    })
    control <- reactive({
        control.dcm(nsteps = input$nsteps, new.mod = inHost)
    })
    
    # Correr el modelo
    mod <- reactive({
        input$runMod
        isolate(dcm(param(), init(), control()))
    })
    
    # Preparar datasets para graficar
    res <- reactive({
        lista <- mod()$epi
        res <- as.data.frame(lista)
        colnames(res) <- names(lista)
        res$Tiempo <- 1:nrow(res)
        return(res)
    })
    res.x.y <- reactive({
        res() %>% 
            gather(key = Tipo, value = Numero, x, y) %>% 
            select(Tiempo, Tipo, Numero) %>% 
            mutate(Tipo = factor(Tipo, levels = c("x", "y"), labels = c("No infectadas", "Infectadas")))
    })
    res.x.y.prev <- reactive({
        res() %>% 
            mutate(prevalencia.x = 1 - prevalencia) %>% 
            gather(key = Tipo, value = Prevalencia, prevalencia, prevalencia.x) %>% 
            select(Tiempo, Tipo, Prevalencia) %>% 
            mutate(Tipo = factor(Tipo, levels = c("prevalencia.x", "prevalencia"), labels = c("No infectadas", "Infectadas")))
    })
    res.a.d <- reactive({
        res() %>% 
            gather(key = Tipo, value = Tasa, ax.flow, ay.flow, ds.flow, di.flow) %>% 
            mutate(Tipo = str_remove(Tipo, ".flow"), 
                   Tipo = str_replace(Tipo, "s", "x"),
                   Tipo = str_replace(Tipo, "i", "y")) %>% 
            separate(Tipo, into = c("Direccion", "Tipo"), sep = 1) %>% 
            mutate(Direccion = factor(Direccion, levels = c("a", "d"), labels = c("Crecimiento (arrival)", "Eliminación (departure)")),
                   Tipo = factor(Tipo, levels = c("x", "y"), labels = c("Células no infectadas", "Células infectadas")))
    })
    
    ## Pestaña del grafico principal
    output$MainPlot <- renderPlot({

        if (input$compsel == "Población de células infectada y no infectada") {
            g = ggplot(res.x.y(), aes(x = Tiempo, y = Numero, group = Tipo, color = Tipo)) + 
                geom_line(lwd = 1.5) + 
                theme_bw() + 
                theme(legend.position = "bottom")
        }
        if (input$compsel == "Prevalencia") {
            g = ggplot(res.x.y.prev(), aes(x = Tiempo, y = Prevalencia, group = Tipo, color = Tipo)) + 
                geom_line(lwd = 1.5) + 
                theme_bw() + 
                theme(legend.position = "bottom")
        }
        if (input$compsel == "Incidencia") {
            g = ggplot(res(), aes(x = Tiempo, y = si.flow)) + 
                geom_line(lwd = 1.5) + 
                scale_y_continuous("Incidencia") +
                theme_bw()
        }
        if (input$compsel == "Tasa de aparición y eliminación de células") {
            g = ggplot(res.a.d(), aes(x = Tiempo, y = Tasa)) + 
                geom_line() + 
                facet_grid(Direccion ~ Tipo, scales = "free") + 
                theme_bw()
        }
        if (input$compsel == "Tamaño de la población") {
            g = ggplot(res(), aes(x = Tiempo, y  = num)) +
                geom_line() +
                scale_y_continuous("Tamaño de la población", limits = c(0, NA)) + 
                theme_bw()
        }
        g 
    }, height = 400, width = 600)
    
    ## Pestaña con el resumen del modelo
    
    # Impresion del modelo
    output$outSummary1 <- renderPrint({
        print(mod())
    })
    
    # Estadisticas eligiendo un tiempo
    output$outSummary2 <- renderPrint({
        res2 = round(res()[input$summTs, 1:8], input$summDig)
        lineas = paste0(
            "----------------------------------------------------
    Estadísticas del modelo para Tiempo = ", input$summTs, "\n",
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
        writeLines(lineas)
    })
    
    ## Pestaña de los datos
    output$outData <- renderDataTable({
        round(as.data.frame(mod()), input$tabdig)
    }, options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 10))
    
    output$dlData <- downloadHandler(
        filename = "ModelData.csv",
        content = function(file) {
            write.csv(as.data.frame(mod()), file)
        }
    )
    
    #-----------------------------------
    # Pestaña con el analisis
    #-----------------------------------
    
    # Preparar parametros
    
    # Crear una lista con todos los parametros para poder usarlos cuando necesite
    params <- reactive({
            input$runMod
            isolate(
            c(lambda = input$lambda, v1 = input$v1, v2 = input$v2, mu1 = input$mu1, mu2 = input$mu2,
                    K = input$K, beta = input$beta, sigma = input$sigma, nsteps = input$nsteps, x = input$x, y = input$y)
            )
    })
    
    # Analisis de equilibrio, estabilidad y dinamica global segun el libro de Li
    #----------------------------------------------------------------------------
    analisisRtdo <- reactive({
        input$runMod
        isolate(
            analisis(lambda = input$lambda, v1 = input$v1, v2 = input$v2, mu1 = input$mu1, mu2 = input$mu2, K = input$K, beta = input$beta, sigma = input$sigma, nsteps = input$nsteps) 
        )
    })
    
    output$analisis <- renderPrint({
        writeLines(analisisRtdo()$lineas)
    })
    
    output$f1f2 <- renderPlot({
        analisisRtdo()$g
    }, height = 400, width = 600)
    
    # Analisis de equilibrio, estabilidad con el paquete phaseR
    #----------------------------------------------------------------------------
    
    # Encontrar los puntos de equilibrio y su clasificacion
    puntosEq <- reactive({
        input$runMod
        isolate({
            nroEq <- analisisRtdo()$cantidad
            encontrados <- 0
            puntos <- matrix(NA, nrow = nroEq, ncol = 2)
            clases <- character(nroEq)
            ids <- character(nroEq)
            it <- 0
            while(encontrados < nroEq & it < 500) {
                it = it + 1
                y0 <-  runif(2, 0, input$nsteps)
                eq <- findEquilibrium(sistema, parameters = params(), y0 = y0)
                nuevo <- as.vector(eq$ystar)
                nuevoClase <- eq$classification
                id <- toString(round(nuevo, 2))
                if (all(nuevo >= 0)) {
                    if (encontrados == 0) {
                        encontrados = encontrados + 1
                        puntos[1, ] = nuevo
                        clases[1] = nuevoClase
                        ids[1] = id
                    } else {
                        if ((!id %in% ids)) {
                            encontrados = encontrados + 1
                            puntos[encontrados, ] = nuevo
                            clases[encontrados] = nuevoClase
                            ids[encontrados] = id
                        }
                    }    
                }
            }
            rtdo <- as.data.frame(round(puntos, input$cifrasEq))
            colnames(rtdo) <- c("Coord.X.NoInf", "Coord.Y.Inf")
            rtdo$Clasificacion <- clases
            cbind(Punto = 1:nrow(rtdo), rtdo)
        })
    })
    output$outDataEq <- renderDataTable({
        puntosEq()
    })
    
    # definir el campo vectorial para graficarlo despues
    flowField.aux <- reactive({
        flowField(sistema,
                  xlim = c(0, input$nsteps),
                  ylim = c(0, input$nsteps),
                  parameters = params(),
                  points = 40,
                  # asp = 1,
                  add = F
                  )
    })

    # Grafico
    output$planoFases <- renderPlot({
        # Voy armando el grafico, primero el campo vectorial, luego se le agregan las nullclines y por ultimo los puntos
        flowField.aux()
        nullclines(sistema, xlim = c(0, input$nsteps), ylim = c(0, input$nsteps), parameters = params())
        points(puntosEq()$Coord.X.NoInf, puntosEq()$Coord.Y.Inf, pch = 19)
    }, height = 400, width = 600)
}

#--------------------------------------
# RUN - EJECUTAR
#--------------------------------------

shinyApp(ui = ui, server = server)
