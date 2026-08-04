# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Parámetros de la simulación
set.seed(123)

a <- 3
b <- 2
B <- 1000

# Muestra Monte Carlo
theta_mc_3 <- rgamma(B, shape = a, rate = b)

# Resúmenes --------------------------------------------------------------------

# Media
media_theta <- mean(theta_mc_3)

round(media_theta, 3)

# Desviación estándar
desv_theta <- sd(theta_mc_3)

round(desv_theta, 3)

# Error estándar de Monte Carlo
error_estandar <- desv_theta / sqrt(length(theta_mc_3))

round(error_estandar, 3)

# Coeficiente de variación de Monte Carlo
cv_mc <- error_estandar / abs(media_theta)

round(cv_mc, 3)

# Tamaño de muestra necesario para un margen de error de 0.01
epsilon     <- 0.01
B_necesario <- (1.96 * desv_theta / epsilon)^2

round(B_necesario)

# Fin --------------------------------------------------------------------------