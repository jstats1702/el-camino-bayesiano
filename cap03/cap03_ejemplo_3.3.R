# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Datos ------------------------------------------------------------------------

n      <- 1000
y      <- 532
theta0 <- 0.5

# Enfoque frecuentista ---------------------------------------------------------

# Estadístico de prueba usando aproximación Normal
z <- (y - n * theta0) / sqrt(n * theta0 * (1 - theta0))

round(z, 3)

# Valor p bilateral aproximado
p_val <- 2 * pnorm(q = abs(z), lower.tail = FALSE)

round(p_val, 3)

# Enfoque Bayesiano ------------------------------------------------------------

# Probabilidades previas de las hipótesis
p_H0 <- 0.5
p_H1 <- 0.5

# Bajo H0: theta = 0.5
py_H0 <- dbinom(x = y, size = n, prob = theta0)

round(py_H0, 4)

# Bajo H1: theta | H1 ~ Beta(a, b)
a <- 1
b <- 1

# Verosimilitud marginal bajo H1 usando la distribución Beta-Binomial
py_H1 <- exp(
     lchoose(n, y) +
          lbeta(a + y, b + n - y) -
          lbeta(a, b)
)

round(py_H1, 4)

# Forma equivalente para a = b = 1
round(1 / (n + 1), 4)

# Factor de Bayes a favor de H0 frente a H1
B01 <- py_H0 / py_H1

round(B01, 3)

# Factor de Bayes a favor de H1 frente a H0
B10 <- 1 / B01

round(B10, 3)

# Probabilidades posteriores de las hipótesis ---------------------------------

post_H0 <- (py_H0 * p_H0) / (py_H0 * p_H0 + py_H1 * p_H1)
post_H1 <- (py_H1 * p_H1) / (py_H0 * p_H0 + py_H1 * p_H1)

round(post_H0, 3)
round(post_H1, 3)

# Fin --------------------------------------------------------------------------