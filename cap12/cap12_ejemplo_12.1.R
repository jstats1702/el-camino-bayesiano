# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Cargar la librería necesaria
suppressMessages(suppressWarnings(library(gtools)))

# Parámetros
k <- 10000  # Número de puntos

alpha_values <- c(0.1, 1, 10, 100)  # Diferentes valores de alpha

G0 <- function(x) pnorm(x)  # Medida base: función de distribución Normal estándar

set.seed(123)

# Simulación para diferentes valores de alpha
for (alpha in alpha_values) {
     pdf(
          file      = paste0("dp_simulacion_particion_alpha_", alpha, ".pdf"),
          width     = 5,
          height    = 5,
          pointsize = 15
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     # Configurar el gráfico
     plot(
          NA,
          NA,
          xlim = c(-3, 3),
          ylim = c(0, 1),
          xlab = "x",
          ylab = "G(x)",
          main = ""
     )
     
     # Generar múltiples realizaciones
     for (l in 1:10) {
          # Simular y ordenar los valores de x
          x <- sort(
               runif(
                    n   = k,
                    min = -3,
                    max = 3
               )
          )
          
          # Calcular los parámetros de concentración
          a <- numeric(k + 1)
          
          a[1]     <- alpha * G0(x[1])
          a[k + 1] <- alpha * (1 - G0(x[k]))
          
          for (j in 2:k) {
               a[j] <- alpha * (G0(x[j]) - G0(x[j - 1]))
          }
          
          # Simular a partir de la distribución Dirichlet
          u <- gtools::rdirichlet(
               n     = 1,
               alpha = a
          )
          
          # Graficar la suma acumulada de los pesos simulados
          lines(
               x,
               cumsum(u)[-length(u)],
               type = "l",
               col  = which(alpha_values == alpha)
          )
     }
     
     # Agregar la curva de la medida base
     curve(
          G0,
          from = -3,
          to   = 3,
          n    = 1000,
          lwd  = 2,
          add  = TRUE
     )
     
     dev.off()
}

# Fin --------------------------------------------------------------------------