# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Datos ------------------------------------------------------------------------

# Datos: length of stay
LoS <- c(1, 2, 1, 1, 4, 1, 2, 2, 0, 3, 6, 2, 1, 3)

(n <- length(LoS))
(y <- sum(LoS))

# Distribución previa ----------------------------------------------------------

a <- 1
b <- 1 / 3

# Media previa
round(a / b, 3)

# Coeficiente de variación previo
round(1 / sqrt(a), 3)

# Distribución posterior -------------------------------------------------------

# Parámetros posteriores
(ap <- a + y)
(bp <- b + n)

# Inferencia -------------------------------------------------------------------

# Media posterior
round(ap / bp, 3)

# Coeficiente de variación posterior
round(1 / sqrt(ap), 3)

# Intervalo de credibilidad al 95%
round(qgamma(p = c(0.025, 0.975), shape = ap, rate = bp), 3)

# Prueba unilateral: H0: theta >= 3 vs H1: theta < 3 --------------------------

theta0 <- 3

# Probabilidades posteriores de las hipótesis
post_H0 <- pgamma(q = theta0, shape = ap, rate = bp, lower.tail = FALSE)
post_H1 <- pgamma(q = theta0, shape = ap, rate = bp, lower.tail = TRUE)

# P(H0 | y)
round(post_H0, 3)

# P(H1 | y)
round(post_H1, 3)

# Probabilidades previas de las hipótesis
prior_H0 <- pgamma(q = theta0, shape = a, rate = b, lower.tail = FALSE)
prior_H1 <- pgamma(q = theta0, shape = a, rate = b, lower.tail = TRUE)

# P(H0)
round(prior_H0, 3)

# P(H1)
round(prior_H1, 3)

# Factor de Bayes a favor de H1 frente a H0
B10 <- (post_H1 / post_H0) / (prior_H1 / prior_H0)

round(B10, 3)

# Factor de Bayes a favor de H0 frente a H1
B01 <- 1 / B10

round(B01, 3)

# Fin --------------------------------------------------------------------------