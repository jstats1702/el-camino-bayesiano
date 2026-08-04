# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Cargar la librería necesaria
suppressMessages(suppressWarnings(library(gtools)))

# Parámetros
alpha <- 10

H_values <- c(5, 10, 50)

G0 <- function(y) punif(y, min = 0, max = 1)

set.seed(123)

# Simular una aproximación Dirichlet finita para cada valor de H
for (H in H_values) {
     # Simular los átomos
     vartheta <- runif(
          n   = H,
          min = 0,
          max = 1
     )
     
     # Simular los pesos
     omega <- as.numeric(
          gtools::rdirichlet(
               n     = 1,
               alpha = rep(alpha / H, H)
          )
     )
     
     # Ordenar los átomos y los pesos correspondientes
     ord          <- order(vartheta)
     vartheta_ord <- vartheta[ord]
     omega_ord    <- omega[ord]
     
     # Construir la función de distribución de G_H
     G_H <- stepfun(
          x     = vartheta_ord,
          y     = c(0, cumsum(omega_ord)),
          right = FALSE
     )
     
     # Crear un archivo independiente para cada valor de H
     pdf(
          file = paste0(
               "dp_aproximacion_dirichlet_finita_H_",
               H,
               ".pdf"
          ),
          width     = 5,
          height    = 5,
          pointsize = 18
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     # Graficar la realización de G_H
     plot(
          G_H,
          xlim      = c(0, 1),
          ylim      = c(0, 1),
          do.points = FALSE,
          verticals = TRUE,
          xlab      = "y",
          ylab      = expression(G[H](y)),
          main      = "",
          lwd       = 1,
          col       = 4
     )
     
     # Agregar la función de distribución base
     curve(
          G0,
          from = 0,
          to   = 1,
          n    = 1000,
          lwd  = 2,
          add  = TRUE
     )
     
     dev.off()
}

# Fin --------------------------------------------------------------------------