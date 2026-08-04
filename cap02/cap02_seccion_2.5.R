# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Tratamiento de datos ---------------------------------------------------------

# Datos
df <- read.delim("victimas.txt", stringsAsFactors = FALSE)

# Dimensión de la base de datos
dim(df)

# Variables disponibles
names(df)

# Frecuencias de sexo
table(df$sexo, useNA = "ifany")

# Proporción sin información en sexo
round(mean(df$sexo == "Sin Informacion", na.rm = TRUE), 4)

# Codificación de sexo
df <- df[df$sexo != "Sin Informacion", ]

df$sexo <- ifelse(df$sexo == "Mujer", 1, ifelse(df$sexo == "Hombre", 0, NA))

df$sexo <- as.numeric(df$sexo)

# Sexo en el año 2016
y <- df[df$agno == 2016, "sexo"]

# Frecuencias de sexo en el año 2016
table(y)

# Tamaño de muestra
(n <- length(y))

# Estadístico suficiente
(s <- sum(y))

# Modelamiento -----------------------------------------------------------------

# Hiperparámetros
a <- 1
b <- 1

# Parámetros de la posterior
(ap <- a + s)
(bp <- b + n - s)

# Visualización
pdf(file = "victimas_posterior.pdf", width = 5, height = 5, pointsize = 15)

par(mar = c(3.25, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

# Gráfico de la distribución Beta posterior
curve(
     expr = dbeta(x, shape1 = ap, shape2 = bp),
     from = 0,
     to   = 1,
     n    = 10000,
     col  = "red",
     lwd  = 2,
     xlab = expression(theta),
     ylab = expression(paste("p(", theta, " | ", y, ")")),
     main = ""
)

# Distribución previa Beta(1,1)
curve(
     expr = dbeta(x, shape1 = a, shape2 = b),
     from = 0,
     to   = 1,
     n    = 10000,
     col  = "royalblue",
     lwd  = 2,
     add  = TRUE
)

# Leyenda
legend(
     "topleft",
     legend = c("Previa", "Posterior"),
     col    = c("royalblue", "red"),
     lty    = 1,
     lwd    = 2,
     bty    = "n"
)

dev.off()

# Inferencia -------------------------------------------------------------------

# Cálculo de estadísticas posteriores
media    <- ap / (ap + bp)
mediana  <- qbeta(0.5, shape1 = ap, shape2 = bp)
moda     <- (ap - 1) / (ap + bp - 2)
varianza <- (ap * bp) / ((ap + bp)^2 * (ap + bp + 1))
cv       <- sqrt(varianza) / media
ic_perc  <- qbeta(p = c(0.025, 0.975), shape1 = ap, shape2 = bp)

# Crear tabla de resultados
out <- data.frame(
     Media       = media,
     Mediana     = mediana,
     Moda        = moda,
     CV          = cv,
     `Q2.5%`     = ic_perc[1],
     `Q97.5%`    = ic_perc[2],
     check.names = FALSE
)

round(out, 4)

# Visualización
pdf(file = "victimas_intervalos.pdf", width = 5, height = 5, pointsize = 15)

par(mar = c(3.25, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

curve(
     expr = dbeta(x, shape1 = ap, shape2 = bp),
     from = 0.7,
     to   = 1,
     n    = 10000,
     col  = "red",
     lwd  = 2,
     xlab = expression(theta),
     ylab = expression(paste("p(", theta, " | ", y, ")")),
     main = ""
)

# Intervalo de credibilidad basado en percentiles del 95%
ic_perc <- qbeta(
     p      = c(0.025, 0.975),
     shape1 = ap,
     shape2 = bp
)

abline(
     v   = ic_perc,
     col = "darkgreen",
     lty = 4,
     lwd = 2
)

# Media posterior
media <- ap / (ap + bp)

abline(
     v   = media,
     col = "black",
     lty = 3,
     lwd = 2
)

dev.off()

# Fin --------------------------------------------------------------------------