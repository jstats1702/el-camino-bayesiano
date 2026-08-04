# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Datos ------------------------------------------------------------------------

df <- read.csv("CNPV2018.txt", stringsAsFactors = FALSE)

# Dimensiones
dim(df)

# Recodificación del nivel educativo ------------------------------------------

# 0: Sin educación superior
# 1: Con educación superior
df$educ_sup <- ifelse(
     df$P_NIVEL_ANOSR %in% c(8, 9),
     1,
     ifelse(is.na(df$P_NIVEL_ANOSR) | df$P_NIVEL_ANOSR == 99, NA, 0)
)

# Frecuencias: indicador de educación superior
table(df$educ_sup, useNA = "ifany")

round(100 * prop.table(table(df$educ_sup, useNA = "ifany")), 1)

# Hijos(as) nacidos vivos ------------------------------------------------------

# PA1_THNV: Total de hijos(as) nacidos vivos
table(df$PA1_THNV, useNA = "ifany")

# Recodificación: NA se interpreta como 0 hijos por salto del cuestionario
df$hijos <- as.numeric(replace(df$PA1_THNV, is.na(df$PA1_THNV), 0))

# Frecuencias: número de hijos
table(df$hijos, useNA = "ifany")

# 0 hijos
sum(df$hijos == 0)
round(100 * sum(df$hijos == 0) / length(df$hijos), 3)

# 1 o más hijos
sum((df$hijos != 0) & (df$hijos != 99))
round(100 * sum((df$hijos != 0) & (df$hijos != 99)) / length(df$hijos), 3)

# No informa hijos
sum(df$hijos == 99)
round(100 * sum(df$hijos == 99) / length(df$hijos), 3)

# Remover datos faltantes o no informados -------------------------------------

df <- subset(df, !is.na(educ_sup) & hijos != 99)

# Filtro de selección ----------------------------------------------------------

filtro <- with(
     df,
     (P_PARENTESCOR == 1) &
          (P_SEXO == 2) &
          (P_EDADR == 9) &
          (PA1_GRP_ETNIC == 6) &
          (PA_LUG_NAC %in% c(2, 3)) &
          (PA_VIVIA_5ANOS %in% c(2, 3)) &
          (PA_HNV %in% c(1, 2)) &
          (P_ALFABETA == 1)
)

# Número de hijos según nivel educativo ---------------------------------------

y1 <- df$hijos[filtro & df$educ_sup == 0]  # Sin educación superior
y2 <- df$hijos[filtro & df$educ_sup == 1]  # Con educación superior

# Tamaños de muestra
(n1 <- length(y1))
(n2 <- length(y2))

# Estadísticos suficientes
(s1 <- sum(y1))
(s2 <- sum(y2))

# Distribución de frecuencias --------------------------------------------------

pdf(file = "hijos_barras.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

y_vals <- 0:6

freq_y1 <- table(factor(y1, levels = y_vals)) / n1
freq_y2 <- table(factor(y2, levels = y_vals)) / n2

freq_mat <- rbind(
     "Sin" = freq_y1,
     "Con" = freq_y2
)

barplot(
     freq_mat,
     beside      = TRUE,
     ylim        = c(0, 0.4),
     names.arg   = y_vals,
     ylab        = "F. relativa",
     xlab        = "No. de hijos",
     legend.text = rownames(freq_mat),
     args.legend = list(
          bty    = "n",
          x      = "topright",
          border = c(2, 4)
     ),
     col         = c(2, 4),
     border      = NA
)

dev.off()

# Previa Gamma(2, 1) -----------------------------------------------------------

a <- 2
b <- 1

# Media previa de theta
round(a / b, 3)

# Coeficiente de variación previo de theta
round(1 / sqrt(a), 3)

# Parámetros posteriores -------------------------------------------------------

# Sin educación superior
(ap1 <- a + s1)
(bp1 <- b + n1)

# Con educación superior
(ap2 <- a + s2)
(bp2 <- b + n2)

# Visualización ----------------------------------------------------------------

pdf(file = "hijos_posterior.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

# Distribución previa y distribuciones posteriores
theta <- seq(0, 5, length.out = 1000)

plot(
     NA,
     NA,
     xlim = c(0, 4),
     ylim = c(0, 5.5),
     xlab = expression(theta),
     ylab = expression(paste("p(", theta, " | ", y, ")")),
     main = ""
)

lines(theta, dgamma(theta, shape = ap1, rate = bp1), col = 2, lwd = 2)
lines(theta, dgamma(theta, shape = ap2, rate = bp2), col = 4, lwd = 2)
lines(theta, dgamma(theta, shape = a, rate = b), col = 1, lwd = 1)

abline(h = 0, col = 1)

legend(
     "topright",
     legend = c("Sin", "Con", "Previa"),
     bty    = "n",
     lwd    = c(2, 2, 1),
     col    = c(2, 4, 1)
)

dev.off()

# Inferencia -------------------------------------------------------------------

# Media posterior de theta
theta_hat_1 <- ap1 / bp1
theta_hat_2 <- ap2 / bp2

# Coeficiente de variación posterior de theta
cv_1 <- 1 / sqrt(ap1)
cv_2 <- 1 / sqrt(ap2)

# Intervalo de credibilidad al 95% para theta
ic95_1 <- qgamma(p = c(0.025, 0.975), shape = ap1, rate = bp1)
ic95_2 <- qgamma(p = c(0.025, 0.975), shape = ap2, rate = bp2)

# Probabilidad posterior de theta > 2
pr_theta_1 <- pgamma(q = 2, shape = ap1, rate = bp1, lower.tail = FALSE)
pr_theta_2 <- pgamma(q = 2, shape = ap2, rate = bp2, lower.tail = FALSE)

# Tabla de resultados ----------------------------------------------------------

tab <- rbind(
     "Sin superior" = c(theta_hat_1, cv_1, ic95_1, pr_theta_1),
     "Con superior" = c(theta_hat_2, cv_2, ic95_2, pr_theta_2)
)

colnames(tab) <- c("Media", "CV", "Q2.5%", "Q97.5%", "Pr. > 2")

round(tab, 3)

# Prueba de hipótesis ----------------------------------------------------------

# H0: theta1 = theta2, con theta ~ Gamma(a1, b1)
# H1: theta1 ~ Gamma(a1, b1), theta2 ~ Gamma(a2, b2)

# Hiperparámetros bajo H1: theta1 ~ Gamma(a1, b1), theta2 ~ Gamma(a2, b2)
a1 <- a
b1 <- b
a2 <- a
b2 <- b

# Factor de Bayes B01 en escala logarítmica
log_B01 <- (
     a1 * log(b1) - lgamma(a1) +
          lgamma(a1 + s1 + s2) -
          (a1 + s1 + s2) * log(b1 + n1 + n2)
) - (
     a1 * log(b1) - lgamma(a1) +
          lgamma(a1 + s1) -
          (a1 + s1) * log(b1 + n1) +
          a2 * log(b2) - lgamma(a2) +
          lgamma(a2 + s2) -
          (a2 + s2) * log(b2 + n2)
)

# Factores de Bayes
B01 <- exp(log_B01)
B10 <- exp(-log_B01)

round(B01, 6)
round(B10, 3)

# Probabilidades posteriores, asumiendo P(H0) = P(H1) = 0.5
post_H1 <- B10 / (1 + B10)
post_H0 <- 1 / (1 + B10)

round(post_H1, 6)
round(post_H0, 6)

# Valor p frecuentista ---------------------------------------------------------

# Valor p = Pr(observar una diferencia tan extrema o más extrema bajo H0)
yb1 <- mean(y1)
yb2 <- mean(y2)

sd1 <- sd(y1)
sd2 <- sd(y2)

z_welch <- (yb1 - yb2) / sqrt(sd1^2 / n1 + sd2^2 / n2)
p_welch <- 2 * pnorm(q = abs(z_welch), lower.tail = FALSE)

round(z_welch, 3)
p_welch

# Fin --------------------------------------------------------------------------