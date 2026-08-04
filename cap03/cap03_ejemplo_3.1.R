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

# Prueba puntual: H0: theta = 3 vs H1: theta != 3 ------------------------------

theta0 <- 3

# Probabilidades previas de las hipótesis
p_H0 <- 0.5
p_H1 <- 0.5

# Verosimilitud marginal bajo H0
py_H0 <- dpois(x = y, lambda = n * theta0)

round(py_H0, 4)

# Verosimilitud marginal bajo H1
py_H1 <- gamma(a + y) / (gamma(a) * gamma(y + 1)) *
     (b / (b + n))^a *
     (n / (b + n))^y

round(py_H1, 4)

# Factor de Bayes a favor de H1 frente a H0
B10 <- py_H1 / py_H0

round(B10, 3)

# Factor de Bayes a favor de H0 frente a H1
B01 <- py_H0 / py_H1

round(B01, 3)

# Probabilidades posteriores de las hipótesis
post_H1 <- (py_H1 * p_H1) / (py_H0 * p_H0 + py_H1 * p_H1)
post_H0 <- (py_H0 * p_H0) / (py_H0 * p_H0 + py_H1 * p_H1)

round(post_H1, 3)
round(post_H0, 3)

# Fin --------------------------------------------------------------------------