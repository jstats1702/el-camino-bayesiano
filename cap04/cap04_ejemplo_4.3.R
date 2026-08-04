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
     df$P_NIVEL_ANOSR %in% c(8, 9), 1,
     ifelse(is.na(df$P_NIVEL_ANOSR) | df$P_NIVEL_ANOSR == 99, NA, 0)
)

# Frecuencias del indicador de educación superior
table(df$educ_sup, useNA = "ifany")
round(100 * prop.table(table(df$educ_sup, useNA = "ifany")), 1)

# Hijos(as) nacidos vivos ------------------------------------------------------

# PA1_THNV: total de hijos(as) nacidos vivos
table(df$PA1_THNV, useNA = "ifany")

# Recodificación: NA se interpreta como 0 hijos por salto del cuestionario
df$hijos <- as.numeric(replace(df$PA1_THNV, is.na(df$PA1_THNV), 0))

# Frecuencias del número de hijos
table(df$hijos, useNA = "ifany")

# 0 hijos
sum(df$hijos == 0)
round(100 * sum(df$hijos == 0) / length(df$hijos), 3)

# 1 o más hijos
sum(df$hijos != 0 & df$hijos != 99)
round(100 * sum(df$hijos != 0 & df$hijos != 99) / length(df$hijos), 3)

# No informa hijos
sum(df$hijos == 99)
round(100 * sum(df$hijos == 99) / length(df$hijos), 3)

# Remover datos faltantes o no informados -------------------------------------

df <- subset(df, !is.na(educ_sup) & hijos != 99)

# Filtro de selección ----------------------------------------------------------

filtro <- with(
     df,
     P_PARENTESCOR == 1 &
          P_SEXO == 2 &
          P_EDADR == 9 &
          PA1_GRP_ETNIC == 6 &
          PA_LUG_NAC %in% c(2, 3) &
          PA_VIVIA_5ANOS %in% c(2, 3) &
          PA_HNV %in% c(1, 2) &
          P_ALFABETA == 1
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

# Simulación de Monte Carlo ----------------------------------------------------

# Número de muestras
B <- 10000

# Muestras posteriores de theta_1 y theta_2
set.seed(123)
th1_mc <- rgamma(B, shape = ap1, rate = bp1)
th2_mc <- rgamma(B, shape = ap2, rate = bp2)

# Muestra posterior de Delta = theta_1 - theta_2
Delta <- th1_mc - th2_mc

# Resúmenes posteriores de Delta ----------------------------------------------

# Media posterior de Delta
round(mean(Delta), 3)

# Coeficiente de variación posterior de Delta
round(sd(Delta) / abs(mean(Delta)), 3)

# Intervalo de credibilidad al 95% para Delta
round(quantile(Delta, probs = c(0.025, 0.975)), 3)

# Probabilidad posterior de que Delta > 0
round(mean(Delta > 0), 3)

# Gráficos ---------------------------------------------------------------------

pdf(file = "hijos_diferencia.pdf", width = 5, height = 5, pointsize = 15)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

hist(
     Delta,
     breaks = 30,
     freq   = FALSE,
     col    = "gold2",
     border = "gold2",
     xlab   = expression(Delta),
     ylab   = expression(paste("p(", Delta, " | ", y, ")")),
     main   = ""
)

lines(density(Delta), col = "black", lwd = 2)

dev.off()

# Simulación predictiva posterior ---------------------------------------------

# Muestras de la distribución predictiva posterior
set.seed(123)
y1_mc <- rpois(B, lambda = th1_mc)
y2_mc <- rpois(B, lambda = th2_mc)

# Muestra predictiva posterior de d = y*_1 - y*_2
d_mc <- y1_mc - y2_mc

# Media predictiva posterior de d
round(mean(d_mc), 3)

# Intervalo predictivo posterior al 95% para d
round(quantile(d_mc, probs = c(0.025, 0.975)), 3)

# Probabilidad predictiva posterior de que d > 0
round(mean(d_mc > 0), 3)

# Gráfico de las distribuciones predictivas posteriores ------------------------

pdf(file = "hijos_predictiva.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(2, 0.75, 0))

y_vals <- 0:6

pred_y1 <- dnbinom(y_vals, size = ap1, mu = ap1 / bp1)
pred_y2 <- dnbinom(y_vals, size = ap2, mu = ap2 / bp2)

pred_mat <- rbind(
     "Sin" = pred_y1,
     "Con" = pred_y2
)

barplot(
     pred_mat,
     beside      = TRUE,
     ylim        = c(0, 0.4),
     names.arg   = y_vals,
     xlab        = expression(y^"*"),
     ylab        = expression(p(y^"*" ~ "|" ~ bold(y))),
     main        = "",
     legend.text = rownames(pred_mat),
     args.legend = list(
          bty    = "n",
          x      = "topright",
          border = c(2, 4)
     ),
     col         = c(2, 4),
     border      = NA
)

dev.off()

# Gráfico de la diferencia predictiva posterior --------------------------------

pdf(file = "hijos_predictiva_diferencia.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(2, 0.75, 0))

d_vals <- min(d_mc):max(d_mc)

d_freq <- table(factor(d_mc, levels = d_vals)) / B

barplot(
     d_freq,
     ylim      = c(0, max(d_freq) * 1.1),
     names.arg = d_vals,
     xlab      = expression(d == y[1]^"*" - y[2]^"*"),
     ylab      = expression(p(d ~ "|" ~ bold(y))),
     main      = "",
     col       = "darkgrey",
     border    = NA
)

dev.off()

# Fin --------------------------------------------------------------------------