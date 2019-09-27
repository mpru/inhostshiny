---
title       : "In-Host Shiny: análisis de modelos compartimentales in-host para infecciones virales en poblaciones celulares"
subtitle: "Marcos Prunello"
author      :
    - "Universidad Nacional de Rosario"
job: "@mrqtsp | github.com/mpru | marcosprunello@gmail.com"
logo        : LOGO-UNR-NEGRO.png
framework   : io2012        # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      # 
widgets     : [mathjax, quiz, bootstrap, interactive] # {mathjax, quiz, bootstrap}
ext_widgets : {rCharts: [libraries/nvd3, libraries/leaflet, libraries/dygraphs]}
mode        : selfcontained # {standalone, draft}
knit        : slidify::knit2slides
assets      : {assets: ../../assets}
--- .class #id

<style type="text/css">
body {background:grey transparent;
}
</style>






## Introducción

- Una **enfermedad** es *infecciosa* si el agente causante, por ejemplo, un virus o bacteria, puede ser transmitido de un individuo a otro a través de algún medio de contagio (contacto directo, vía respiratorio, madre a recién nacido, etc.).
- El objetivo de la **modelización matemática** de una enfermedad infecciosa es describir el *proceso de transmisión* de la enfermedad, que se puede resumir como:
    
    1. Se introducen individuos infectados en una población susceptible --> transmisión. 
    2. Un individuo puede ser asintomático durante una etapa temprana. 
    3. Individuos infectados pueden o no recuperarse y ganar cierto grado de inmunidad.
    4. Según los mecanismos en juego: *brote*, *epidemia*, *pandemia*.

- La **modelización matemática** de procesos infecciosos ha sido capaz de proveer claridad sobre su transmisión y expansión, ayudando a estimar severidad y diseñar prevención.


--- .class #id

## Introducción

- Este mismo enfoque puede aplicarse para describir la **diseminación de una infección viral en una población de células**, resultando en modelos conocidos como **in-host**.
- En este trabajo se presenta una aplicación Shiny que permite estudiar el proceso infeccioso por el virus HTLV-I en células T CD4+, responsable de varias enfermedades como linfoma/leucemia y mielopatía.

--- .class #id

## Introducción

- Organización:
    
    1. Descripción de los modelos determinísticos compartimentales SIR.
    2. Adaptación de los modelos para el estudio de una población de células.
    3. Presentación de la aplicación.

--- .segue bg:grey

# 1. Modelos SIR

--- .class #id

<!--
## Modelos determinísticos compartimentales SIR

- Se basan en la solución de sistemas de ecuaciones diferenciales que representan sistemas epidémicos a través del tiempo.
- Son **determinísticos** porque las soluciones son funciones fijas de los parámetros y las condiciones iniciales, sin estudiar variabilidad existente en el proceso o población.
- Son **compartimentales** porque particionan a la población en grupos mutuamente exclusivos que representan distintos estados de la enfermedad. Por ejemplo:

    - **S**: individuos susceptibles
    - **I**: individuos infectados
    - **R**: individuos recuperados
-->

## Modelos determinísticos compartimentales SIR

<img src="assets/img/diagrama1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" width="70%" style="display: block; margin: auto;" />

- El objetivo de la modelización es llevar registro de las cantidades de individuos en cada compartimento en cualquier momento $t$: $S(t)$, $I(t)$ y $R(t)$. 
- Se establece un pequeño intervalo de tiempo $[t, t+\Delta t]$ para evaluar los cambios en dichas cantidades: $\Delta S(t)$, $\Delta I(t)$  y $\Delta R(t)$.
- $\Delta S(t) / \Delta t = S'(t)$ con $\Delta t \rightarrow 0$: podemos plantear un sistema de ecuaciones diferenciales en términos de incidencia, tasa de recuperación, tasa de remoción, etc.

<!--
Estas variaciones pueden calcularse como:
- $\Delta S(t)$ = nuevos susceptibles + ingresos desde R - nuevas infecciones - egresos desde S
- $\Delta I(t)$ = nuevos infectados (desde S) - egresos hacia R - egresos desde I
- $\Delta R(t)$ = ingresos desde I - egresos hacia S - egresos desde R
-->

--- .class #id

## Modelos determinísticos compartimentales SIR

<img src="assets/img/diagrama2.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" width="120%" style="display: block; margin: auto;" />

<!--
proportional rate assumption that the birth or death rate is proportional to the population size. 
Here b is the natural birth rate constant, and d1, d2, and d3 are
death rate constants for compartments S, I, and R, respectively.
Rate constant d2 may include both natural and disease-caused
death.

Para describir los movimientos poblacionales, se deben establecer supuestos sobre las tasas de transferencias entre compartimentos. Por ejemplo:

- Tasa de recuperación (o tasa de transferencia de $I$ a $R$): $\gamma I$
- Tasa de transferencia de $R$ a $S$: $\delta R$, pensando en que los individuos recuperados tienen cierto período de inmunidad.
- Incidencia (número de nuevas infecciones por unidad de tiempo): $\beta I S$, llamada **incidencia bilineal**. Se basa en el supuesto de que la mezcla o intercambio entre individuos es homogéneo y por lo tanto aplica la *Ley de la acción de masas*: el número de contacto depende del número de individuos en cada compartimento. $\beta$ recibe el nombre de **coeficiente de transmisión**.
- Tasas de nacimiento, muerte, inmigración y emigración para modelar el crecimiento de la población.

Para algunas enfermedades en las cuales se da un periodo de latencia entre la infección y la manifestación de los síntomas, se puede pensar en una modificación del modelo como lo muestra la siguiente figura. El compartimento de los infectados se separa en dos, el **compartimento latente** $E$ y el infeccioso $I$, y se puede también asumir que la transferencia de $E$ a $I$ es proprocional al tamaño de $E$, $\epsilon E$:
    
Cuando la transmisión ocurre por contacto directo entre individuos infectados y susceptibles se dice que la transmisión en **horizontal**. Sin embargo, pueden existir otros medios, como la **transmisión vertical**, cuando los patógenos pasan a un recién nacido directamente de un progenitor infectado. El siguiente diagrama muestra un modelo para una enfermedad con transmisión tanto horizontal como vertical:
-->

<!--
--- .class #id

## Modelos determinísticos compartimentales SIR

- **Número básico de reproducción** $\mathbf{R}_0$: promedio de infecciones secundarias causadas por un único individuo infectado en una población enteramente susceptible durante un período infeccioso medio. 
- $\mathbf{R}_0 < 1 \implies$ no habrá epidemia; $\mathbf{R}_0 > 1 \implies$ habrá epidemia.
- $\mathbf{R}_0$ puede deducirse a partir del análisis de la estabilidad de los puntos de equilibrio libres de enfermedad que presente el sistema.
-->

--- .segue bg:grey

# 2. Adaptación para poblaciones celulares

--- .class #id

## Modelo in-host para infecciones por HTLV-I en células T CD4+

- El mismo tipo de enfoque puede aplicarse para describir la diseminación de una infección viral en una población de células: modelos **in-host**
- Aplicación: análisis de la infección de células T CD4+ por el virus HTLV-I (responsable de linfoma/leucemia, mielopatía, etc.) 
- Es un **retrovirus**, fuera de la célula no genera infección, requiere contacto célula a célula, se transmite verticalmente.

<!-- emplea una proteina viral para sintetizar segmentos de ADN que se integran al ADN de la célula hospedadora -->

--- .class #id

## Modelo in-host para infecciones por HTLV-I en células T CD4+

<img src="assets/img/diagrama3.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" width="70%" style="display: block; margin: auto;" />

<!--
Parametros:
El cuerpo humano genera células T CD4+ a una tasa de λ=25 células/mm3.
La tasa natural de muerte es μ1≈μ2≈0.03 días−1.
En ausencia de infección, el número de células se espera que sea constante con valor normal alrededor de 1000/mm3 y una capacidad de carga constante igual a K=1150/mm3. Esto implica fijar x0≈1000/mm3, es decir, f1(1000)≈0, de lo cual se deduce una proliferación constante de v2≈v1≈0.038 día−1.
β≈1.03×10−3mm3/célula/día.
-->

--- .segue bg:grey

# 3. In-Host Shiny App

--- .class #id

## In-Host Shiny App

- Permite estudiar el proceso infeccioso por el virus HTLV-I en células T CD4+.
- Interés: llega la infección en algún momento a algún nivel de estabilidad? Se extinguirá o se expanderá? 
- Para responder esto, la app estudia y clasifica los puntos de equilibrio $(x, y)$, que son estados en los que sistema no cambia.
- Un equilibrio es (asintóticamente) **estable** si el sistema siempre vuelve a ese punto luego de pequeñas perturbaciones o inestable si no vuelve.

--- .segue bg:grey

# Muchas gracias

### Marcos Prunello.

### @mrqtsp | github.com/mpru

### mpru.shinyapps.io/inhostshiny
