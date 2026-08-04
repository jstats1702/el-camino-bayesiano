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

# Estadísticos observados ------------------------------------------------------

t_obs_1 <- c(media = mean(y1), de = sd(y1))
t_obs_2 <- c(media = mean(y2), de = sd(y2))

round(t_obs_1, 3)
round(t_obs_2, 3)

# Simulación de Monte Carlo ----------------------------------------------------

# Número de muestras
B <- 10000

# Muestras posteriores de theta_1 y theta_2
set.seed(123)
th1_mc <- rgamma(B, shape = ap1, rate = bp1)
th2_mc <- rgamma(B, shape = ap2, rate = bp2)

# Inicializar matrices para almacenar estadísticos de prueba
set.seed(123)
t_mc_1 <- matrix(NA, nrow = B, ncol = 2)
t_mc_2 <- matrix(NA, nrow = B, ncol = 2)

colnames(t_mc_1) <- c("media", "de")
colnames(t_mc_2) <- c("media", "de")

# Distribución predictiva posterior
set.seed(123)

for (i in seq_len(B)) {
     # Datos replicados
     y1_rep <- rpois(n = n1, lambda = th1_mc[i])
     y2_rep <- rpois(n = n2, lambda = th2_mc[i])
     
     # Estadísticos de prueba
     t_mc_1[i, ] <- c(mean(y1_rep), sd(y1_rep))
     t_mc_2[i, ] <- c(mean(y2_rep), sd(y2_rep))
}

# Valores p predictivos posteriores -------------------------------------------

ppp_1 <- c(
     media = mean(t_mc_1[, "media"] <= t_obs_1["media"]),
     de    = mean(t_mc_1[, "de"] <= t_obs_1["de"])
)

ppp_2 <- c(
     media = mean(t_mc_2[, "media"] <= t_obs_2["media"]),
     de    = mean(t_mc_2[, "de"] <= t_obs_2["de"])
)

round(ppp_1, 3)
round(ppp_2, 3)

# Visualización ----------------------------------------------------------------

# Colores
col1 <- 2
col2 <- 4

# Límites de los ejes
xlim_media <- range(
     t_mc_1[, "media"],
     t_mc_2[, "media"],
     t_obs_1["media"],
     t_obs_2["media"]
)

xlim_de <- range(
     t_mc_1[, "de"],
     t_mc_2[, "de"],
     t_obs_1["de"],
     t_obs_2["de"]
)

ylim_media <- c(0, 4)
ylim_de    <- c(0, 6.5)

# Histograma de la media, grupo sin educación superior -------------------------

pdf(file = "hijos_chequeo_media_sin.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

hist(
     x      = t_mc_1[, "media"],
     freq   = FALSE,
     col    = col1,
     border = col1,
     xlim   = xlim_media,
     ylim   = ylim_media,
     xlab   = "t",
     ylab   = expression(p(t ~ "|" ~ bold(y))),
     main   = ""
)

abline(v = t_obs_1["media"], col = 1, lwd = 2)

legend(
     "topleft",
     legend = c("Sin", "Con", "Obs"),
     bty    = "n",
     fill   = c(2, 4, 1),
     border = c(2, 4, 1),
     col    = c(2, 4, 1)
)

dev.off()

# Histograma de la media, grupo con educación superior -------------------------

pdf(file = "hijos_chequeo_media_con.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

hist(
     x      = t_mc_2[, "media"],
     freq   = FALSE,
     col    = col2,
     border = col2,
     xlim   = xlim_media,
     ylim   = ylim_media,
     xlab   = "t",
     ylab   = expression(p(t ~ "|" ~ bold(y))),
     main   = ""
)

abline(v = t_obs_2["media"], col = 1, lwd = 2)

dev.off()

# Histograma de la desviación estándar, grupo sin educación superior -----------

pdf(file = "hijos_chequeo_de_sin.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

hist(
     x      = t_mc_1[, "de"],
     freq   = FALSE,
     col    = col1,
     border = col1,
     xlim   = xlim_de,
     ylim   = ylim_de,
     xlab   = "t",
     ylab   = expression(p(t ~ "|" ~ bold(y))),
     main   = ""
)

abline(v = t_obs_1["de"], col = 1, lwd = 2)

legend(
     "topleft",
     legend = c("Sin", "Con", "Obs"),
     bty    = "n",
     fill   = c(2, 4, 1),
     border = c(2, 4, 1),
     col    = c(2, 4, 1)
)

dev.off()

# Histograma de la desviación estándar, grupo con educación superior -----------

pdf(file = "hijos_chequeo_de_con.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

hist(
     x      = t_mc_2[, "de"],
     freq   = FALSE,
     col    = col2,
     border = col2,
     xlim   = xlim_de,
     ylim   = ylim_de,
     xlab   = "t",
     ylab   = expression(p(t ~ "|" ~ bold(y))),
     main   = ""
)

abline(v = t_obs_2["de"], col = 1, lwd = 2)

dev.off()

# Fin --------------------------------------------------------------------------