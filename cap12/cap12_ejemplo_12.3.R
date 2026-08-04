# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Cargar la librería necesaria
suppressMessages(suppressWarnings(library(stats)))

# Parámetros
alpha <- 10
eps   <- 1e-6

# Seleccionar el nivel de truncación
H <- ceiling(
     log(eps) /
          log(alpha / (alpha + 1))
)

# Distribución base
G0 <- function(x) {
     pnorm(x)
}

set.seed(123)

# Simular los átomos
vartheta <- rnorm(
     n    = H,
     mean = 0,
     sd   = 1
)

# Simular las variables del proceso de partición secuencial
V <- rbeta(
     n      = H,
     shape1 = 1,
     shape2 = alpha
)

# Calcular los pesos antes de reasignar la masa residual
omega <- V * c(
     1,
     cumprod(1 - V[-H])
)

# Calcular la masa residual después de H particiones
R_H <- prod(1 - V)

# Asignar la masa residual al último peso
omega[H] <- omega[H] + R_H

# Evaluar la función de distribución de la realización
x_vals <- seq(
     from       = -3,
     to         = 3,
     length.out = 10000
)

G <- sapply(
     X   = x_vals,
     FUN = function(x) {
          sum(omega[vartheta <= x])
     }
)

# Visualización de los pesos
pdf(
     file      = "dp_particion_secuencial_pesos.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Panel izquierdo: átomos y pesos
plot(
     x    = vartheta,
     y    = omega,
     type = "h",
     xlim = c(-3, 3),
     xlab = expression(vartheta),
     ylab = expression(omega),
     main = "",
     col  = 4,
     lend = 1
)

dev.off()

# Visualización de la realización
pdf(
     file      = "dp_particion_secuencial_realizacion.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

# Configurar la ventana gráfica
par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Panel derecho: función de distribución de la realización
plot(
     x    = x_vals,
     y    = G,
     type = "l",
     xlim = c(-3, 3),
     ylim = c(0, 1),
     xlab = "x",
     ylab = "G(x)",
     main = "",
     col  = 4
)

# Agregar la función de distribución base
curve(
     expr = G0,
     from = -3,
     to   = 3,
     n    = 1000,
     lwd  = 2,
     add  = TRUE
)

dev.off()

# Fin --------------------------------------------------------------------------