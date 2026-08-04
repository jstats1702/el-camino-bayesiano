# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

set.seed(123)

n      <- 500
alphas <- c(0.5, 1, 5, 10)

# Función para simular el CRP --------------------------------------------------

simular_crp <- function(n, alpha) {
     xi      <- integer(n)
     K_i     <- integer(n)
     tamanos <- integer(n)
     
     xi[1]      <- 1L
     K_i[1]     <- 1L
     tamanos[1] <- 1L
     
     numero_mesas <- 1L
     
     if (n >= 2) {
          for (i in 2:n) {
               prob_nueva <- alpha / (alpha + i - 1)
               
               if (runif(1) < prob_nueva) {
                    numero_mesas <- numero_mesas + 1L
                    
                    xi[i]                  <- numero_mesas
                    tamanos[numero_mesas] <- 1L
               } else {
                    mesa <- sample.int(
                         n    = numero_mesas,
                         size = 1,
                         prob = tamanos[seq_len(numero_mesas)]
                    )
                    
                    xi[i]         <- mesa
                    tamanos[mesa] <- tamanos[mesa] + 1L
               }
               
               K_i[i] <- numero_mesas
          }
     }
     
     list(
          xi      = xi,
          K_i     = K_i,
          tamanos = tamanos[seq_len(numero_mesas)]
     )
}

# Simulación para diferentes valores de alpha ----------------------------------

simulaciones_crp <- lapply(
     alphas,
     function(alpha) {
          simular_crp(
               n     = n,
               alpha = alpha
          )
     }
)

names(simulaciones_crp) <- paste0(
     "alpha_",
     alphas
)

# Tabla de asignaciones para alpha = 1 -----------------------------------------

resultado_alpha_1 <- simulaciones_crp[["alpha_1"]]

tabla_asignaciones <- data.frame(
     cliente        = seq_len(15),
     mesa           = resultado_alpha_1$xi[seq_len(15)],
     mesas_ocupadas = resultado_alpha_1$K_i[seq_len(15)]
)

tabla_asignaciones

# Esperanza teórica de K_i -----------------------------------------------------

esperanzas_K <- lapply(
     alphas,
     function(alpha) {
          cumsum(
               alpha / (
                    alpha + 0:(n - 1)
               )
          )
     }
)

K_matriz <- do.call(
     cbind,
     lapply(
          simulaciones_crp,
          function(resultado) {
               resultado$K_i
          }
     )
)

esperanzas_matriz <- do.call(
     cbind,
     esperanzas_K
)

# Evolución del número de mesas ------------------------------------------------

cols <- c(
     "#1B9E77",
     "#D95F02",
     "#7570B3",
     "#E7298A"
)

pdf(
     file      = "crp_evolucion_mesas.pdf",
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
     col  = cols,
     xlab = "Número de clientes",
     ylab = "Número de mesas ocupadas"
)

matlines(
     x   = seq_len(n),
     y   = esperanzas_matriz,
     lty = 2,
     lwd = 2,
     col = cols
)

leyenda_alpha <- legend(
     "topleft",
     legend = paste0("alpha = ", alphas),
     lty    = 1,
     lwd    = 2,
     col    = cols,
     ncol   = 1,
     bty    = "n"
)

legend(
     x      = leyenda_alpha$rect$left + leyenda_alpha$rect$w,
     y      = leyenda_alpha$rect$top,
     legend = c(
          "Simulado",
          "Esperanza teórica"
     ),
     lty   = c(1, 2),
     lwd   = 2,
     col   = 1,
     ncol  = 1,
     xjust = 0,
     yjust = 1,
     bty   = "n"
)

box()

dev.off()

# Tabla del número final de mesas ==============================================

resumen_mesas <- data.frame(
     alpha = alphas,
     K_n_simulado = vapply(
          simulaciones_crp,
          function(resultado) {
               tail(resultado$K_i, 1)
          },
          numeric(1)
     ),
     esperanza_K_n = vapply(
          esperanzas_K,
          function(esperanza) {
               tail(esperanza, 1)
          },
          numeric(1)
     ),
     mesa_mayor = vapply(
          simulaciones_crp,
          function(resultado) {
               max(resultado$tamanos)
          },
          numeric(1)
     )
)

resumen_mesas

# Tamaños finales de las mesas -------------------------------------------------

max_tamano <- max(
     vapply(
          simulaciones_crp,
          function(resultado) {
               max(resultado$tamanos)
          },
          numeric(1)
     )
)

max_numero_mesas <- max(
     vapply(
          simulaciones_crp,
          function(resultado) {
               length(resultado$tamanos)
          },
          numeric(1)
     )
)

for (k in seq_along(alphas)) {
     tamanos <- sort(
          simulaciones_crp[[k]]$tamanos,
          decreasing = TRUE
     )
     
     pdf(
          file = paste0(
               "crp_tamanos_mesas_alpha_",
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
          lwd  = 3,
          col  = cols[k],
          xlim = c(1, max_numero_mesas),
          ylim = c(0, max_tamano),
          xlab = "Mesa ordenada",
          ylab = "Número de clientes",
          main = ""
     )
     
     dev.off()
}

# Fin --------------------------------------------------------------------------