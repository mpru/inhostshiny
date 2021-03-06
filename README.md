# inhostshiny
Shiny App for the implementation in R of the **in-host model**, an epidemiologic compartimental model for the study of the infection of CD4+ T cells by Human T-cell Lymphotropic Virus Type 1. The model is described by Michael Li in section 5.2 of the book **An Introduction to Mathematical Modeling of Infectious Diseases** (Springer, 2018). The evolution of the infectious process, as well as the stability of the equilibrium of the associated systems of differential equations, can be studied in the app by choosing different values for the parameters of the model.


The app is deployed at https://mpru.shinyapps.io/inhostshiny/ or can be run locally doing:


library(shiny)

runGitHub("inhostshiny", "mpru")

In this case, you should have installed the following R packages: EpiModel, phaseR, ggplot2, tidyr, dplyr, stringr and Deriv.

