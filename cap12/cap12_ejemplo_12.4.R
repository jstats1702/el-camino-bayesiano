# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

set.seed(123)

n      <- 500
alphas <- c(0.5, 1, 5, 10)

rG0 <- function(n) {
     rnorm(
          n,
          mean = 0,
          sd   = 1
     )
}

# Función para simular una sucesión de Pólya -----------------------------------

simular_polya <- function(n, alpha, rG0) {
     theta             <- numeric(n)
     valores_distintos <- numeric(n)
     frecuencias       <- integer(n)
     K_n               <- integer(n)
     
     theta[1]             <- rG0(1)
     valores_distintos[1] <- theta[1]
     frecuencias[1]       <- 1L
     numero_grupos        <- 1L
     K_n[1]               <- 1L
     
     if (n >= 2) {
          for (i in 2:n) {
               prob_nuevo <- alpha / (alpha + i - 1)
               
               if (runif(1) < prob_nuevo) {
                    numero_grupos <- numero_grupos + 1L
                    theta[i] <- rG0(1)
                    valores_distintos[numero_grupos] <- theta[i]
                    frecuencias[numero_grupos] <- 1L
               } else {
                    grupo <- sample.int(
                         n    = numero_grupos,
                         size = 1,
                         prob = frecuencias[seq_len(numero_grupos)]
                    )
                    
                    theta[i] <- valores_distintos[grupo]
                    frecuencias[grupo] <- frecuencias[grupo] + 1L
               }
               
               K_n[i] <- numero_grupos
          }
     }
     
     list(
          theta             = theta,
          valores_distintos = valores_distintos[seq_len(numero_grupos)],
          frecuencias       = frecuencias[seq_len(numero_grupos)],
          K_n               = K_n
     )
}

# Simulación para diferentes valores de alpha ----------------------------------

simulaciones <- lapply(
     alphas,
     function(alpha) {
          simular_polya(
               n     = n,
               alpha = alpha,
               rG0   = rG0
          )
     }
)

names(simulaciones) <- paste0("alpha_", alphas)

# Valores distintos y frecuencias para alpha = 1 -------------------------------

simulacion_alpha_1 <- simulaciones[["alpha_1"]]

tabla_grupos <- data.frame(
     valor_distinto = simulacion_alpha_1$valores_distintos,
     frecuencia     = simulacion_alpha_1$frecuencias
)

tabla_grupos <- tabla_grupos[
     order(tabla_grupos$frecuencia, decreasing = TRUE),
]

round(tabla_grupos, 5)

# Evolución del número de valores distintos ------------------------------------

K_matriz <- do.call(
     cbind,
     lapply(
          simulaciones,
          function(resultado) {
               resultado$K_n
          }
     )
)

pdf(
     file      = "polya_evolucion_valores_distintos.pdf",
     width     = 7,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

matplot(
     x    = seq_len(n),
     y    = K_matriz,
     type = "l",
     lty  = 1,
     lwd  = 2,
     col  = seq_along(alphas),
     xlab = "Número de observaciones",
     ylab = "Número de valores distintos"
)

legend(
     "topleft",
     legend = paste0("alpha = ", alphas),
     lty    = 1,
     lwd    = 2,
     col    = seq_along(alphas),
     bty    = "n"
)

dev.off()

# Tamaños de los grupos --------------------------------------------------------

max_tamano <- max(
     vapply(
          simulaciones,
          function(x) {
               max(x$frecuencias)
          },
          numeric(1)
     )
)

for (k in seq_along(alphas)) {
     tamanos <- sort(
          simulaciones[[k]]$frecuencias,
          decreasing = TRUE
     )
     
     pdf(
          file      = paste0(
               "polya_tamanos_grupos_alpha_",
               alphas[k],
               ".pdf"
          ),
          width     = 5,
          height    = 5,
          pointsize = 15
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     plot(
          x    = seq_along(tamanos),
          y    = tamanos,
          type = "h",
          lwd  = 6,
          ylim = c(0, max_tamano),
          xlab = "Grupo ordenado",
          ylab = "Tamaño",
          main = ""
     )
     
     box()
     
     dev.off()
}

# Estimación de la probabilidad de empate --------------------------------------

B_empates <- 10000

probabilidad_estimada <- numeric(length(alphas))

for (k in seq_along(alphas)) {
     indicadores_empate <- logical(B_empates)
     
     for (b in seq_len(B_empates)) {
          muestra <- simular_polya(
               n     = 2,
               alpha = alphas[k],
               rG0   = rG0
          )$theta
          
          indicadores_empate[b] <- muestra[1] == muestra[2]
     }
     
     probabilidad_estimada[k] <- mean(indicadores_empate)
}

probabilidad_teorica <- 1 / (alphas + 1)

tabla_empates <- data.frame(
     alpha                  = alphas,
     probabilidad_estimada  = probabilidad_estimada,
     probabilidad_teorica   = probabilidad_teorica
)

tabla_empates

# Fin --------------------------------------------------------------------------