# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Paquetes
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(sf)))
suppressMessages(suppressWarnings(library(stringi)))

# Cargar datos
dat <- read.csv(
     file             = file.path("Examen_Saber_11_20251.txt"),
     sep              = ";",
     stringsAsFactors = FALSE
)

dim(dat)

# Revisión inicial por departamento
table(dat$estu_depto_reside, useNA = "ifany")

round(100 * table(dat$estu_depto_reside, useNA = "ifany") / nrow(dat), 1)

# Limpiar variable de departamento
dat$estu_depto_reside <- trimws(dat$estu_depto_reside)

# Remover registros sin departamento y registros del exterior
dat <- dat[
     !is.na(dat$estu_depto_reside) &
          dat$estu_depto_reside != "" &
          !dat$estu_depto_reside %in% c("EXTRANJERO", "EXTRANGERO"),
     ,
     drop = FALSE
]

# Verificar departamentos restantes
table(dat$estu_depto_reside)
unique(dat$estu_depto_reside)

# Muestreo proporcional por departamento --------------------------------------

prop_muestra  <- 0.05
departamentos <- sort(unique(dat$estu_depto_reside))

set.seed(123)

idx_muestra <- integer(0)

for (depto in departamentos) {
     idx_depto     <- which(dat$estu_depto_reside == depto)
     n_depto_total <- length(idx_depto)
     
     if (n_depto_total < 100) {
          n_depto <- min(10, n_depto_total)
     } else {
          n_depto <- max(1, floor(n_depto_total * prop_muestra))
     }
     
     idx_muestra <- c(
          idx_muestra,
          sample(idx_depto, size = n_depto, replace = FALSE)
     )
}

# Base final
dat <- dat[idx_muestra, , drop = FALSE]

dim(dat)
table(dat$estu_depto_reside)

# Estadísticos suficientes -----------------------------------------------------

# m : número de grupos (departamentos)
deptos <- sort(unique(dat$estu_cod_reside_depto))
m      <- length(deptos)

# n : número de individuos (estudiantes)
n <- nrow(dat)

# Tratamiento de datos
# y  : puntaje de los estudiantes (c)
# Y  : puntaje de los estudiantes por departamento (list)
# g  : identificador secuencial de los departamentos (c)
# nj : número de estudiantes por departamento (c)
# yb : promedios por departamento (c)
# s2 : varianzas por departamento (c)

y <- dat$punt_matematicas
Y <- vector(mode = "list", length = m)
g <- rep(NA_integer_, n)

for (j in seq_len(m)) {
     idx    <- dat$estu_cod_reside_depto == deptos[j]
     g[idx] <- j
     Y[[j]] <- y[idx]
}

# Tabla
estadisticos <- dat %>%
     group_by(estu_cod_reside_depto) %>%
     summarise(
          codigo  = first(estu_cod_reside_depto),
          nombre  = first(estu_depto_reside),
          nj      = dplyr::n(),
          yb      = mean(punt_matematicas, na.rm = TRUE),
          med     = median(punt_matematicas, na.rm = TRUE),
          s2      = var(punt_matematicas, na.rm = TRUE),
          s       = sd(punt_matematicas, na.rm = TRUE),
          min     = min(punt_matematicas, na.rm = TRUE),
          max     = max(punt_matematicas, na.rm = TRUE),
          .groups = "drop"
     ) %>%
     arrange(codigo) %>%
     select(codigo, nombre, nj, yb, med, s2, s, min, max)

as.data.frame(estadisticos[, c(2, 3, 4, 5, 7, 8, 9)])

# Vectores resumen
nj <- estadisticos$nj
yb <- estadisticos$yb
s2 <- estadisticos$s2

# Promedio global
round(mean(y), 3)

# Mapa con puntajes promedio ---------------------------------------------------

# Leer GeoJSON de departamentos de Colombia
url_geojson <- paste0(
     "https://raw.githubusercontent.com/caticoa3/colombia_mapa/master/",
     "co_2018_MGN_DPTO_POLITICO.geojson"
)

map_col <- sf::st_read(url_geojson, quiet = TRUE)

# Preparar nombres del mapa
map_col <- map_col %>%
     mutate(
          nombre_key = DPTO_CNMBR %>%
               stringi::stri_trans_general("Latin-ASCII") %>%
               toupper()
     )

# Preparar datos departamentales
estadisticos_mapa <- estadisticos %>%
     mutate(
          nombre_key = nombre %>%
               stringi::stri_trans_general("Latin-ASCII") %>%
               toupper(),
          nombre_key = case_when(
               nombre_key == "BOGOTA"          ~ "BOGOTA, D.C.",
               nombre_key == "VALLE"           ~ "VALLE DEL CAUCA",
               nombre_key == "NORTE SANTANDER" ~ "NORTE DE SANTANDER",
               TRUE                            ~ nombre_key
          )
     )

# Unir mapa con datos
map_col_yb <- map_col %>%
     left_join(
          estadisticos_mapa,
          by = "nombre_key"
     )

# Revisar departamentos sin emparejar
map_col_yb %>%
     filter(is.na(yb)) %>%
     select(DPTO_CNMBR, nombre_key)

# Excluir San Andrés y Providencia para evitar espacio vacío a la izquierda
map_col_yb_cont <- map_col_yb %>%
     filter(
          !grepl(
               pattern     = "SAN ANDRES|PROVIDENCIA|SANTA CATALINA",
               x           = nombre_key,
               ignore.case = TRUE
          )
     )

# Límites del mapa
bbox <- sf::st_bbox(map_col_yb_cont)

pad_x <- 0.35
pad_y <- 0.35

# Mapa
p_mapa <- ggplot(map_col_yb_cont) +
     geom_sf(
          aes(fill = yb),
          color     = "gray40",
          linewidth = 0.15
     ) +
     scale_fill_gradientn(
          colors   = c("#f3e8ff", "#d8b4fe", "#c084fc", "#a855f7", "#7e22ce"),
          na.value = "gray90",
          name     = expression(bar(y)[j])
     ) +
     coord_sf(
          xlim   = c(bbox["xmin"] - pad_x, bbox["xmax"] + pad_x),
          ylim   = c(bbox["ymin"] - pad_y, bbox["ymax"] + pad_y),
          expand = FALSE
     ) +
     labs(
          x = "Longitud",
          y = "Latitud"
     ) +
     theme_bw(base_size = 13) +
     theme(
          legend.position = "right",
          panel.border = element_rect(
               color     = "black",
               fill      = NA,
               linewidth = 0.6
          ),
          panel.grid.major = element_line(
               color     = "gray85",
               linewidth = 0.25
          ),
          panel.grid.minor = element_blank(),
          axis.title       = element_text(size = 12),
          axis.text        = element_text(size = 10),
          plot.title       = element_text(hjust = 0),
          plot.subtitle    = element_text(hjust = 0),
          plot.caption     = element_text(hjust = 1),
          plot.margin      = margin(t = 6, r = 6, b = 6, l = 6)
     )

p_mapa

# Guardar el mapa
ggsave(
     filename = "matematicas_modelo_jerarquico_normal_mapa_promedio.pdf",
     plot     = p_mapa,
     width    = 5.8,
     height   = 7.2,
     units    = "in"
)

dev.off()

# Ranking frecuentista ---------------------------------------------------------

pdf(
     file      = "matematicas_modelo_jerarquico_normal_ranking_promedio.pdf",
     width     = 5.8,
     height    = 7.2,
     pointsize = 13
)

# Ranking basado en el promedio muestral
par(
     mfrow = c(1, 1),
     mar   = c(4, 9.5, 1.5, 1),
     mgp   = c(2.5, 0.75, 0)
)

ord <- order(yb)

# Tonos violeta
col_puntos <- adjustcolor("#7e22ce", alpha.f = 0.5)
col_media  <- "#7e22ce"
col_linea  <- adjustcolor("#a855f7", alpha.f = 0.45)
col_ref    <- adjustcolor("gray", alpha.f = 0.85)

# Inicializar gráfico vacío
plot(
     x    = c(0, 100),
     y    = c(0.5, m + 0.5),
     type = "n",
     xlab = "Puntaje",
     ylab = "",
     main = "",
     yaxt = "n",
     yaxs = "i"
)

# Líneas horizontales de referencia
abline(
     h   = 1:m,
     col = "lightgray",
     lwd = 1
)

# Línea vertical en el promedio nacional
abline(
     v   = 50,
     col = col_ref,
     lwd = 3
)

# Puntos individuales por departamento según el orden del promedio
for (l in 1:m) {
     j <- ord[l]
     
     points(
          x   = Y[[j]],
          y   = rep(l, nj[j]),
          pch = 16,
          cex = 0.4,
          col = col_puntos
     )
}

# Conectar promedios muestrales
lines(
     x    = yb[ord],
     y    = 1:m,
     type = "l",
     col  = col_linea,
     lwd  = 1.5
)

# Añadir puntos del promedio muestral
points(
     x   = yb[ord],
     y   = 1:m,
     pch = 16,
     cex = 1.1,
     col = col_media
)

# Etiquetas con nombres de los departamentos
axis(
     side   = 2,
     at     = 1:m,
     labels = estadisticos$nombre[ord],
     las    = 2
)

dev.off()

# Histograma de los promedios --------------------------------------------------

pdf(
     file      = "matematicas_modelo_jerarquico_normal_histograma_promedios.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Histograma del promedio por grupo
hist(
     x      = yb,
     freq   = FALSE,
     breaks = 15,
     col    = adjustcolor("#d8b4fe", alpha.f = 0.65),
     border = "white",
     xlab   = "Promedio",
     ylab   = "Densidad",
     main   = ""
)

abline(
     v   = mean(y),
     col = "gray",
     lwd = 3
)

dev.off()

# Dispersograma del promedio frente al tamaño ---------------------------------

pdf(
     file      = "matematicas_modelo_jerarquico_normal_promedio_vs_tamano.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Diagrama de dispersión: tamaño del grupo vs. promedio
plot(
     x    = nj,
     y    = yb,
     xlab = "Tamaño del grupo",
     ylab = "Promedio",
     pch  = 16,
     cex  = 1.2,
     col  = adjustcolor("#7e22ce", alpha.f = 0.60)
)

abline(
     h   = mean(y, na.rm = TRUE),
     col = "gray",
     lwd = 3
)

dev.off()

# Distribución previa ----------------------------------------------------------

# Hiperparámetros
mu0 <- 50
g20 <- 10^2

eta0 <- 1
t20  <- 10^2

nu0 <- 1
s20 <- 10^2

# Muestreador de Gibbs ---------------------------------------------------------

mcmc <- function(B, y, nj, yb, s2, mu0, g20, eta0, t20, nu0, s20) {
     # Frecuencia para imprimir progreso
     ncat <- max(1, floor(B / 10))
     
     # Número total de observaciones y grupos
     n <- sum(nj)
     m <- length(nj)
     
     # Almacenamiento de la cadena
     THETA <- matrix(NA_real_, nrow = B, ncol = m + 4)
     
     # Valores iniciales
     theta <- yb
     sig2  <- mean(s2, na.rm = TRUE)
     mu    <- mean(theta, na.rm = TRUE)
     tau2  <- var(theta, na.rm = TRUE)
     
     # Cadena de MCMC
     for (b in seq_len(B)) {
          # Actualizar theta_j
          v_theta <- 1 / (1 / tau2 + nj / sig2)
          m_theta <- v_theta * (mu / tau2 + nj * yb / sig2)
          theta   <- rnorm(n = m, mean = m_theta, sd = sqrt(v_theta))
          
          # Actualizar sigma^2
          ss_sig2 <- sum((nj - 1) * s2 + nj * (yb - theta)^2)
          a_sig2  <- 0.5 * (nu0 + n)
          b_sig2  <- 0.5 * (nu0 * s20 + ss_sig2)
          sig2    <- 1 / rgamma(n = 1, shape = a_sig2, rate = b_sig2)
          
          # Actualizar mu
          v_mu <- 1 / (1 / g20 + m / tau2)
          m_mu <- v_mu * (mu0 / g20 + sum(theta) / tau2)
          mu   <- rnorm(n = 1, mean = m_mu, sd = sqrt(v_mu))
          
          # Actualizar tau^2
          ss_tau2 <- sum((theta - mu)^2)
          a_tau2  <- 0.5 * (eta0 + m)
          b_tau2  <- 0.5 * (eta0 * t20 + ss_tau2)
          tau2    <- 1 / rgamma(n = 1, shape = a_tau2, rate = b_tau2)
          
          # Log-verosimilitud usando estadísticos suficientes por grupo
          ll <- sum(
               -0.5 * nj * log(2 * pi * sig2) -
                    0.5 * ((nj - 1) * s2 + nj * (yb - theta)^2) / sig2
          )
          
          # Almacenar iteración
          THETA[b, ] <- c(theta, sig2, mu, tau2, ll)
          
          # Imprimir progreso
          if (b %% ncat == 0) {
               cat(100 * round(b / B, 1), "% completado ... \n", sep = "")
          }
     }
     
     # Salida
     colnames(THETA) <- c(paste0("theta", seq_len(m)), "sig2", "mu", "tau2", "ll")
     THETA <- as.data.frame(THETA)
     
     return(list(THETA = THETA))
}

# Ajuste del modelo ------------------------------------------------------------

set.seed(123)
chain <- mcmc(B = 10000, y, nj, yb, s2, mu0, g20, eta0, t20, nu0, s20)

# Función para graficar cadenas ------------------------------------------------

plot_cadena <- function(x, ylab, col, cex, file) {
     pdf(
          file      = file,
          width     = 6,
          height    = 4,
          pointsize = 15
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3.1, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     plot(
          x    = seq_along(x),
          y    = x,
          type = "p",
          pch  = 16,
          col  = adjustcolor(col, alpha.f = 0.1),
          cex  = cex,
          xlab = "Iteración",
          ylab = ylab,
          main = ""
     )
     
     dev.off()
}

# Convergencia -----------------------------------------------------------------

# Cadenas para la log-verosimilitud
plot_cadena(
     x    = chain$THETA$ll,
     ylab = "Log-verosimilitud",
     col  = 1,
     cex  = 0.5,
     file = "matematicas_modelo_jerarquico_normal_cadena_logverosimilitud.pdf"
)

# Tamaños efectivos de muestra
neff <- coda::effectiveSize(
     coda::as.mcmc(chain$THETA)
)

# Error estándar de Monte Carlo
EEMC <- apply(
     X      = chain$THETA,
     MARGIN = 2,
     FUN    = sd,
     na.rm  = TRUE
) / sqrt(neff)

# Coeficiente de variación de Monte Carlo
media_post <- colMeans(
     chain$THETA,
     na.rm = TRUE
)

CVMC <- EEMC / abs(media_post)

round(summary(neff), digits = 1)
round(summary(EEMC), digits = 4)
round(summary(CVMC), digits = 4)

# Graficar posterior -----------------------------------------------------------

plot_posterior <- function(x, xlab, file, legend = TRUE, legend_pos = "topleft") {
     pdf(
          file      = file,
          width     = 5,
          height    = 5,
          pointsize = 15
     )
     
     on.exit(dev.off())
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     hist(
          x      = x,
          freq   = FALSE,
          breaks = 30,
          col    = adjustcolor("gray90", alpha.f = 0.65),
          border = "white",
          xlab   = xlab,
          ylab   = "Densidad",
          main   = ""
     )
     
     abline(
          v   = mean(x),
          col = 2,
          lwd = 2,
          lty = 2
     )
     
     abline(
          v   = quantile(x, probs = c(0.025, 0.975)),
          col = 4,
          lwd = 2,
          lty = 2
     )
     
     if (legend) {
          legend(
               x      = legend_pos,
               legend = c("Media", "IC 95%"),
               col    = c(2, 4),
               fill   = c(2, 4),
               border = c(2, 4),
               bty    = "n"
          )
     }
}

# Inferencia -------------------------------------------------------------------

# Cadenas de eta, mu, sigma y tau
PAR <- cbind(
     eta   = chain$THETA$sig2 / (chain$THETA$sig2 + chain$THETA$tau2),
     mu    = chain$THETA$mu,
     sigma = sqrt(chain$THETA$sig2),
     tau   = sqrt(chain$THETA$tau2)
)

# Posterior mu
plot_posterior(
     x          = PAR[, "mu"],
     xlab       = expression(mu),
     file       = "matematicas_modelo_jerarquico_normal_posterior_mu.pdf",
     legend     = TRUE,
     legend_pos = "topleft"
)

# Posterior tau
plot_posterior(
     x          = PAR[, "tau"],
     xlab       = expression(tau),
     file       = "matematicas_modelo_jerarquico_normal_posterior_tau.pdf",
     legend     = FALSE,
     legend_pos = "topleft"
)

# Posterior sigma
plot_posterior(
     x          = PAR[, "sigma"],
     xlab       = expression(sigma),
     file       = "matematicas_modelo_jerarquico_normal_posterior_sigma.pdf",
     legend     = FALSE,
     legend_pos = "topleft"
)

# Posterior eta
plot_posterior(
     x          = PAR[, "eta"],
     xlab       = expression(eta),
     file       = "matematicas_modelo_jerarquico_normal_posterior_eta.pdf",
     legend     = FALSE,
     legend_pos = "topleft"
)

# Resumen posterior: media, CV (%) y cuantiles
tab <- cbind(
     `Media`  = colMeans(PAR),
     `CV`     = abs(apply(PAR, 2, sd) / colMeans(PAR)),
     `Q2.5%`  = apply(PAR, 2, quantile, probs = 0.025),
     `Q97.5%` = apply(PAR, 2, quantile, probs = 0.975)
)

round(tab, 3)

# Ranking Bayesiano ------------------------------------------------------------

pdf(
     file      = "matematicas_modelo_jerarquico_normal_ranking_bayesiano.pdf",
     width     = 5.8,
     height    = 7.2,
     pointsize = 13
)

# Ranking Bayesiano
par(
     mfrow = c(1, 1),
     mar   = c(4, 9.5, 1.5, 1),
     mgp   = c(2.5, 0.75, 0)
)

# Resúmenes posteriores
THETA_theta <- as.matrix(chain$THETA[, seq_len(m)])

ids2 <- estadisticos$nombre
that <- colMeans(THETA_theta)

ic1 <- apply(
     X      = THETA_theta,
     MARGIN = 2,
     FUN    = quantile,
     probs  = c(0.025, 0.975)
)

# Ordenar por media posterior
ord <- order(that)

ids2 <- ids2[ord]
that <- that[ord]
ic1  <- ic1[, ord]

# Colores según posición del intervalo frente a 50
colo <- rep(2, m)
colo[ic1[1, ] > 50] <- 1
colo[ic1[2, ] < 50] <- 3

colo <- c("royalblue", "black", "red")[colo]

# Color de referencia
col_ref <- adjustcolor("gray", alpha.f = 0.85)

# Inicializar gráfico vacío
plot(
     x    = c(0, 100),
     y    = c(0.5, m + 0.5),
     type = "n",
     xlab = "Puntaje",
     ylab = "",
     main = "",
     yaxt = "n",
     yaxs = "i"
)

# Líneas horizontales de referencia
abline(
     h   = 1:m,
     col = "lightgray",
     lwd = 1
)

# Línea vertical en el promedio de referencia
abline(
     v   = 50,
     col = col_ref,
     lwd = 3
)

# Intervalos de credibilidad y medias posteriores
for (j in 1:m) {
     segments(
          x0  = ic1[1, j],
          y0  = j,
          x1  = ic1[2, j],
          y1  = j,
          col = colo[j],
          lwd = 1.5
     )
     
     points(
          x   = that[j],
          y   = j,
          pch = 16,
          cex = 1.1,
          col = colo[j]
     )
}

# Etiquetas con nombres de los departamentos
axis(
     side   = 2,
     at     = 1:m,
     labels = ids2,
     las    = 2
)

dev.off()

# CV de los theta --------------------------------------------------------------

# Coeficiente de variación posterior de los theta_j
that <- apply(
     X      = chain$THETA[, seq_len(m)],
     MARGIN = 2,
     FUN    = mean,
     na.rm  = TRUE
)

shat <- apply(
     X      = chain$THETA[, seq_len(m)],
     MARGIN = 2,
     FUN    = sd,
     na.rm  = TRUE
)

cv_b <- abs(shat / that)

round(summary(cv_b), 3)

# Mapa con medias posteriores de los theta -------------------------------------

# Leer GeoJSON de departamentos de Colombia
url_geojson <- paste0(
     "https://raw.githubusercontent.com/caticoa3/colombia_mapa/master/",
     "co_2018_MGN_DPTO_POLITICO.geojson"
)

map_col <- sf::st_read(url_geojson, quiet = TRUE)

# Media posterior de cada theta_j
theta_cols <- paste0("theta", seq_len(m))

theta_post <- colMeans(
     as.matrix(chain$THETA[, theta_cols]),
     na.rm = TRUE
)

# Preparar nombres del mapa
map_col <- map_col %>%
     mutate(
          nombre_key = DPTO_CNMBR %>%
               stringi::stri_trans_general("Latin-ASCII") %>%
               toupper()
     )

# Preparar datos departamentales
estadisticos_mapa <- estadisticos %>%
     mutate(
          theta_post = as.numeric(theta_post),
          nombre_key = nombre %>%
               stringi::stri_trans_general("Latin-ASCII") %>%
               toupper(),
          nombre_key = case_when(
               nombre_key == "BOGOTA"          ~ "BOGOTA, D.C.",
               nombre_key == "VALLE"           ~ "VALLE DEL CAUCA",
               nombre_key == "NORTE SANTANDER" ~ "NORTE DE SANTANDER",
               TRUE                            ~ nombre_key
          )
     )

# Unir mapa con datos
map_col_theta <- map_col %>%
     left_join(
          estadisticos_mapa,
          by = "nombre_key"
     )

# Revisar departamentos sin emparejar
map_col_theta %>%
     filter(is.na(theta_post)) %>%
     select(DPTO_CNMBR, nombre_key)

# Excluir San Andrés y Providencia para evitar espacio vacío a la izquierda
map_col_theta_cont <- map_col_theta %>%
     filter(
          !grepl(
               pattern     = "SAN ANDRES|PROVIDENCIA|SANTA CATALINA",
               x           = nombre_key,
               ignore.case = TRUE
          )
     )

# Límites del mapa
bbox <- sf::st_bbox(map_col_theta_cont)

pad_x <- 0.35
pad_y <- 0.35

# Mapa
p_mapa <- ggplot(map_col_theta_cont) +
     geom_sf(
          aes(fill = theta_post),
          color     = "gray40",
          linewidth = 0.15
     ) +
     scale_fill_gradientn(
          colors   = c("#f3e8ff", "#d8b4fe", "#c084fc", "#a855f7", "#7e22ce"),
          na.value = "gray90",
          name     = expression(hat(theta)[j])
     ) +
     coord_sf(
          xlim   = c(bbox["xmin"] - pad_x, bbox["xmax"] + pad_x),
          ylim   = c(bbox["ymin"] - pad_y, bbox["ymax"] + pad_y),
          expand = FALSE
     ) +
     labs(
          x = "Longitud",
          y = "Latitud"
     ) +
     theme_bw(base_size = 13) +
     theme(
          legend.position = "right",
          panel.border = element_rect(
               color     = "black",
               fill      = NA,
               linewidth = 0.6
          ),
          panel.grid.major = element_line(
               color     = "gray85",
               linewidth = 0.25
          ),
          panel.grid.minor = element_blank(),
          axis.title       = element_text(size = 12),
          axis.text        = element_text(size = 10),
          plot.title       = element_text(hjust = 0),
          plot.subtitle    = element_text(hjust = 0),
          plot.caption     = element_text(hjust = 1),
          plot.margin      = margin(t = 6, r = 6, b = 6, l = 6)
     )

p_mapa

# Guardar el mapa
ggsave(
     filename = "matematicas_modelo_jerarquico_normal_mapa_media_posterior_theta.pdf",
     plot     = p_mapa,
     width    = 5.8,
     height   = 7.2,
     units    = "in"
)

dev.off()

# Contracción ------------------------------------------------------------------

# Estimación posterior vs. promedio muestral
pdf(
     file      = "matematicas_modelo_jerarquico_normal_theta_hat_vs_promedio.pdf",
     width     = 5,
     height    = 5,
     pointsize = 20
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3.3, 1.4, 1.4),
     mgp   = c(2, 0.75, 0)
)

# Medias posteriores de theta_j
theta_hat <- colMeans(
     as.matrix(chain$THETA[, seq_len(m)]),
     na.rm = TRUE
)

plot(
     x    = yb,
     y    = theta_hat,
     xlim = range(yb, theta_hat, na.rm = TRUE),
     ylim = range(yb, theta_hat, na.rm = TRUE),
     xlab = expression(bar(italic(y))[j]),
     ylab = expression(hat(theta)[j]),
     main = "",
     pch  = 16,
     cex  = 1.2,
     col  = adjustcolor("#7e22ce", alpha.f = 0.60)
)

abline(
     a   = 0,
     b   = 1,
     col = "gray",
     lwd = 3
)

dev.off()

# Diferencia entre promedio muestral y media posterior vs. tamaño de muestra
pdf(
     file      = "matematicas_modelo_jerarquico_normal_theta_hat_vs_tamano.pdf",
     width     = 5,
     height    = 5,
     pointsize = 20
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3.3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

d_theta <- yb - theta_hat

plot(
     x    = nj,
     y    = d_theta,
     xlim = range(nj, na.rm = TRUE),
     ylim = c(-1, 1) * max(abs(d_theta), na.rm = TRUE),
     xlab = "Tamaño del grupo",
     ylab = expression(bar(italic(y))[j] - hat(theta)[j]),
     main = "",
     pch  = 16,
     cex  = 1.2,
     col  = adjustcolor("#7e22ce", alpha.f = 0.60)
)

abline(
     h   = 0,
     col = "gray",
     lwd = 3
)

dev.off()

# Contracción: comparación visual entre theta_hat y promedio muestral
pdf(
     file      = "matematicas_modelo_jerarquico_normal_contraccion_theta_promedio.pdf",
     width     = 5,
     height    = 5,
     pointsize = 20
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Medias posteriores de theta_j
theta_hat <- colMeans(
     as.matrix(chain$THETA[, seq_len(m)]),
     na.rm = TRUE
)

# Límites del eje x
x_lim <- range(
     c(yb, theta_hat),
     na.rm = TRUE
)

x_pad <- 0.05 * diff(x_lim)

plot(
     x        = NA,
     y        = NA,
     xlim     = c(x_lim[1] - x_pad, x_lim[2] + x_pad),
     ylim     = c(1, 4),
     xlab     = "Puntaje",
     ylab     = "",
     main     = "",
     yaxt     = "n",
     cex.axis = 0.8
)

axis(
     side   = 2,
     at     = c(2, 3),
     labels = c(expression(hat(theta)[j]), expression(bar(y)[j])),
     las    = 1
)

abline(
     h   = c(2, 3),
     col = c(4, 2),
     lwd = 2
)

for (j in seq_len(m)) {
     segments(
          x0  = theta_hat[j],
          y0  = 2,
          x1  = yb[j],
          y1  = 3,
          col = adjustcolor("black", alpha.f = 0.45),
          lwd = 1
     )
}

dev.off()

# Fin --------------------------------------------------------------------------