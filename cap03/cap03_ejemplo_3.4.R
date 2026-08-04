# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Parámetros -------------------------------------------------------------------

a <- 1
b <- 1
n <- 80
s <- 69

# Parámetros de la distribución posterior
ap <- a + s
bp <- b + n - s

# MAP posterior
theta_map <- (ap - 1) / (ap + bp - 2)

# Varianza aproximada por BCLT
hess <- -(ap - 1) / theta_map^2 - (bp - 1) / (1 - theta_map)^2

round(hess, 4)

var_map <- -1 / hess

round(var_map, 4)

# Puntos para graficar ---------------------------------------------------------

theta <- seq(0.7, 1, length.out = 10000)

posterior_exact  <- dbeta(theta, shape1 = ap, shape2 = bp)
posterior_approx <- dnorm(theta, mean = theta_map, sd = sqrt(var_map))

# Visualización ----------------------------------------------------------------

pdf(file = "victimas_laplace.pdf", width = 5, height = 5, pointsize = 15)

par(mar = c(3.25, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

plot(
     theta,
     posterior_exact,
     type = "l",
     col  = 4,
     lwd  = 2,
     ylab = "Densidad",
     xlab = expression(theta),
     main = ""
)

lines(theta, posterior_approx, col = 2, lwd = 2, lty = 2)

legend(
     "topleft",
     legend = c("Exacta", "Aproximada"),
     col    = c(4, 2),
     lwd    = 2,
     lty    = c(1, 2),
     bty    = "n"
)

dev.off()

# Fin --------------------------------------------------------------------------