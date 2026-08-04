# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Hiperparámetros de la distribución previa
a <- 1
b <- 1

# Tamaño de muestra y número de infectados observados
n <- 20
y <- 1

# Parámetros de la distribución posterior
a_post <- a + y
b_post <- b + n - y

# Media posterior
post_mean <- a_post / (a_post + b_post)

round(post_mean, 3)

# Visualización
pdf(file = "prevalencia_posterior.pdf", width = 5, height = 5, pointsize = 15)

par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

curve(
     expr = dbeta(x, shape1 = a_post, shape2 = b_post),
     from = 0,
     to   = 1,
     n    = 1000,
     col  = "red",
     lwd  = 2,
     xlab = expression(theta),
     ylab = expression(paste("p(", theta, " | ", y, ")")),
     main = ""
)

curve(
     expr = dbeta(x, shape1 = a, shape2 = b),
     from = 0,
     to   = 1,
     n    = 1000,
     col  = "royalblue",
     lwd  = 2,
     add  = TRUE
)

legend(
     "topright",
     legend = c("Posterior", "Previa"),
     col    = c("red", "royalblue"),
     lty    = 1,
     lwd    = 2,
     bty    = "n"
)

dev.off()

# Fin --------------------------------------------------------------------------