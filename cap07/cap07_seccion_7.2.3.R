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

# Distribución previa modelo homocedástico -------------------------------------

# Hiperparámetros
mu0 <- 50
g20 <- 10^2

eta0 <- 1
t20  <- 10^2

nu0 <- 1
s20 <- 10^2

# Muestreador de Gibbs modelo homocedástico ------------------------------------

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

# Ajuste del modelo homocedástico ----------------------------------------------

set.seed(123)
chain_1 <- mcmc(B = 10000, y, nj, yb, s2, mu0, g20, eta0, t20, nu0, s20)

# Distribución previa ----------------------------------------------------------

# Hiperparámetros
mu0 <- 50
g20 <- 10^2

eta0 <- 1
t20  <- 10^2

lam0 <- 1

al0 <- 1
be0 <- 1 / 10^2

nus0 <- 1:50  # Rango para p(nu | rest)

# Muestreador de Gibbs modelo heterocedástico ----------------------------------

mcmc <- function(B, y, nj, yb, s2, mu0, g20, eta0, t20, lam0, al0, be0, nus0) {
     # Frecuencia para imprimir progreso
     ncat <- max(1, floor(B / 10))
     
     # Número total de observaciones y grupos
     n <- sum(nj)
     m <- length(nj)
     
     # Corregir varianzas no definidas o nulas
     s2 <- ifelse(is.na(s2), 0, s2)
     s2 <- pmax(s2, .Machine$double.eps)
     
     # Almacenamiento de la cadena
     THETA <- matrix(NA_real_, nrow = B, ncol = 2 * m + 5)
     
     # Valores iniciales
     theta <- yb
     sig2  <- s2
     mu    <- mean(theta, na.rm = TRUE)
     tau2  <- max(var(theta, na.rm = TRUE), .Machine$double.eps)
     nu    <- min(nus0)
     ups2  <- 100
     
     # Cadena MCMC
     for (b in seq_len(B)) {
          # Actualizar theta_j
          v_theta <- 1 / (1 / tau2 + nj / sig2)
          m_theta <- v_theta * (mu / tau2 + nj * yb / sig2)
          theta   <- rnorm(n = m, mean = m_theta, sd = sqrt(v_theta))
          
          # Actualizar sigma_j^2
          ss_sig2 <- (nj - 1) * s2 + nj * (yb - theta)^2
          a_sig2  <- 0.5 * (nu + nj)
          b_sig2  <- 0.5 * (nu * ups2 + ss_sig2)
          sig2    <- 1 / rgamma(n = m, shape = a_sig2, rate = b_sig2)
          
          # Actualizar mu
          v_mu <- 1 / (1 / g20 + m / tau2)
          m_mu <- v_mu * (mu0 / g20 + sum(theta) / tau2)
          mu   <- rnorm(n = 1, mean = m_mu, sd = sqrt(v_mu))
          
          # Actualizar tau^2
          ss_tau2 <- sum((theta - mu)^2)
          a_tau2  <- 0.5 * (eta0 + m)
          b_tau2  <- 0.5 * (eta0 * t20 + ss_tau2)
          tau2    <- 1 / rgamma(n = 1, shape = a_tau2, rate = b_tau2)
          
          # Actualizar nu
          lpnu <- 0.5 * m * nus0 * log(0.5 * nus0 * ups2) -
               m * lgamma(0.5 * nus0) -
               0.5 * nus0 * sum(log(sig2)) -
               nus0 * (lam0 + 0.5 * ups2 * sum(1 / sig2))
          
          pnu <- exp(lpnu - max(lpnu))
          pnu <- pnu / sum(pnu)
          
          nu <- sample(
               x    = nus0,
               size = 1,
               prob = pnu
          )
          
          # Actualizar ups2
          a_ups2 <- al0 + 0.5 * m * nu
          b_ups2 <- be0 + 0.5 * nu * sum(1 / sig2)
          ups2   <- rgamma(n = 1, shape = a_ups2, rate = b_ups2)
          
          # Log-verosimilitud usando estadísticos suficientes por grupo
          ll <- sum(
               -0.5 * nj * log(2 * pi * sig2) -
                    0.5 * ((nj - 1) * s2 + nj * (yb - theta)^2) / sig2
          )
          
          # Almacenar iteración
          THETA[b, ] <- c(theta, sig2, mu, tau2, nu, ups2, ll)
          
          # Imprimir progreso
          if (b %% ncat == 0) {
               cat(100 * round(b / B, 1), "% completado ... \n", sep = "")
          }
     }
     
     # Salida
     colnames(THETA) <- c(
          paste0("theta", seq_len(m)),
          paste0("sig2", seq_len(m)),
          "mu", "tau2", "nu", "ups2", "ll"
     )
     
     THETA <- as.data.frame(THETA)
     
     return(list(THETA = THETA))
}

# Ajuste del modelo heterocedástico --------------------------------------------

set.seed(123)
chain_2 <- mcmc(B = 10000, y, nj, yb, s2, mu0, g20, eta0, t20, lam0, al0, be0, nus0)

# Función para graficar cadenas ------------------------------------------------

plot_cadena <- function(x, ylim, ylab, col, cex, file) {
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
          ylim = ylim,
          ylab = ylab,
          main = ""
     )
     
     dev.off()
}

# Convergencia -----------------------------------------------------------------

# Cadenas para la log-verosimilitud
plot_cadena(
     x    = chain_1$THETA$ll,
     ylim = range(chain_1$THETA$ll, chain_2$THETA$ll),
     ylab = "Log-verosimilitud",
     col  = 2,
     cex  = 0.5,
     file = "matematicas_modelo_jerarquico_normal_homocedastico_cadena_logverosimilitud.pdf"
)

# Cadenas para la log-verosimilitud
plot_cadena(
     x    = chain_2$THETA$ll,
     ylim = range(chain_1$THETA$ll, chain_2$THETA$ll),
     ylab = "Log-verosimilitud",
     col  = 4,
     cex  = 0.5,
     file = "matematicas_modelo_jerarquico_normal_heterocedastico_cadena_logverosimilitud.pdf"
)

# Tamaños efectivos de muestra
neff <- coda::effectiveSize(
     coda::as.mcmc(chain_2$THETA)
)

# Error estándar de Monte Carlo
EEMC <- apply(
     X      = chain_2$THETA,
     MARGIN = 2,
     FUN    = sd,
     na.rm  = TRUE
) / sqrt(neff)

# Coeficiente de variación de Monte Carlo
media_post <- colMeans(
     chain_2$THETA,
     na.rm = TRUE
)

CVMC <- EEMC / abs(media_post)

round(summary(neff), digits = 1)
round(summary(EEMC), digits = 4)
round(summary(CVMC), digits = 4)

# Comparación de la log-verosimilitud de los modelos
round(mean(chain_1$THETA$ll), 3)
round(mean(chain_2$THETA$ll), 3)

# Graficar posteriores de forma comparativa ------------------------------------

col <- RColorBrewer::brewer.pal(9, "Set1")[1:9]

# Posterior de mu
pdf(
     file      = "matematicas_modelo_jerarquico_normal_heterocedastico_posterior_mu.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

d_mu_1 <- density(chain_1$THETA$mu)
d_mu_2 <- density(chain_2$THETA$mu)

plot(
     d_mu_1,
     col  = col[1],
     lwd  = 3,
     xlab = expression(mu),
     ylab = "Densidad",
     main = "",
     xlim = range(d_mu_1$x, d_mu_2$x),
     ylim = range(0, d_mu_1$y, d_mu_2$y)
)

lines(d_mu_2, col = col[2], lwd = 2)

legend(
     "topright",
     legend = c("Hom.", "Het."),
     col    = col[1:2],
     fill   = col[1:2],
     border = col[1:2],
     bty    = "n"
)

dev.off()

# Posterior de tau
pdf(
     file      = "matematicas_modelo_jerarquico_normal_heterocedastico_posterior_tau.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

d_tau_1 <- density(sqrt(chain_1$THETA$tau2))
d_tau_2 <- density(sqrt(chain_2$THETA$tau2))

plot(
     d_tau_1,
     col  = col[1],
     lwd  = 3,
     xlab = expression(tau),
     ylab = "Densidad",
     main = "",
     xlim = range(d_tau_1$x, d_tau_2$x),
     ylim = range(0, d_tau_1$y, d_tau_2$y)
)

lines(d_tau_2, col = col[2], lwd = 2)

legend(
     "topright",
     legend = c("Hom.", "Het."),
     col    = col[1:2],
     fill   = col[1:2],
     border = col[1:2],
     bty    = "n"
)

dev.off()

# Posterior de nu
pdf(
     file      = "matematicas_modelo_jerarquico_normal_heterocedastico_posterior_nu.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

nu_vals   <- 1:30
nu_counts <- table(factor(chain_2$THETA$nu, levels = nus0))
nu_freq   <- as.numeric(nu_counts[as.character(nu_vals)]) / nrow(chain_2$THETA)
nu_freq[is.na(nu_freq)] <- 0

plot(
     nu_vals,
     nu_freq,
     type = "h",
     lwd  = 3,
     col  = col[2],
     xlab = expression(nu),
     ylab = "Frecuencia relativa",
     main = "",
     ylim = range(0, nu_freq)
)

abline(h = 0, col = "lightgray")

dev.off()

# Posterior de sigma
pdf(
     file      = "matematicas_modelo_jerarquico_normal_heterocedastico_posterior_sigma.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

d_sigma <- density(sqrt(chain_2$THETA$ups2))

plot(
     d_sigma,
     col  = col[2],
     lwd  = 3,
     xlab = expression(sigma),
     ylab = "Densidad",
     main = ""
)

dev.off()

# Resumen posterior de sigma^2 global
ups2_post <- sqrt(chain_2$THETA$ups2)

tab <- c(
     Media    = mean(ups2_post, na.rm = TRUE),
     CV       = sd(ups2_post, na.rm = TRUE) / abs(mean(ups2_post, na.rm = TRUE)),
     `Q2.5%`  = quantile(ups2_post, probs = 0.025, na.rm = TRUE),
     `Q97.5%` = quantile(ups2_post, probs = 0.975, na.rm = TRUE)
)

round(tab, 3)

# Máximo a posteriori de nu
table(chain_2$THETA$nu)

# Medias posteriores de theta_j
theta_hat_1 <- colMeans(chain_1$THETA[, seq_len(m)], na.rm = TRUE)
theta_hat_2 <- colMeans(chain_2$THETA[, seq_len(m)], na.rm = TRUE)

# Correlación entre medias posteriores
cor_theta <- cor(theta_hat_1, theta_hat_2, use = "complete.obs")

round(cor_theta, 3)

# Intervalos de credibilidad del 95%
ic_1 <- apply(
     chain_1$THETA[, seq_len(m)],
     MARGIN = 2,
     FUN    = quantile,
     probs  = c(0.025, 0.975),
     na.rm  = TRUE
)

ic_2 <- apply(
     chain_2$THETA[, seq_len(m)],
     MARGIN = 2,
     FUN    = quantile,
     probs  = c(0.025, 0.975),
     na.rm  = TRUE
)

# Longitudes de los intervalos de credibilidad
l_1 <- ic_1[2, ] - ic_1[1, ]
l_2 <- ic_2[2, ] - ic_2[1, ]

# Correlación entre longitudes de intervalos
cor_l <- cor(l_1, l_2, use = "complete.obs")

round(cor_l, 3)

# Mapa con medias posteriores de theta_j bajo el modelo heterocedástico --------

# Leer GeoJSON de departamentos de Colombia
url_geojson <- paste0(
     "https://raw.githubusercontent.com/caticoa3/colombia_mapa/master/",
     "co_2018_MGN_DPTO_POLITICO.geojson"
)

map_col <- sf::st_read(url_geojson, quiet = TRUE)

# Media posterior de cada theta_j
theta_cols <- paste0("theta", seq_len(m))

theta_post <- colMeans(
     as.matrix(chain_2$THETA[, theta_cols, drop = FALSE]),
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
     filename = "matematicas_modelo_jerarquico_normal_heterocedastico_mapa_media_posterior_theta.pdf",
     plot     = p_mapa,
     width    = 5.8,
     height   = 7.2,
     units    = "in"
)

dev.off()

# Mapa con medias posteriores de sigma_j ---------------------------------------

# Leer GeoJSON de departamentos de Colombia
url_geojson <- paste0(
     "https://raw.githubusercontent.com/caticoa3/colombia_mapa/master/",
     "co_2018_MGN_DPTO_POLITICO.geojson"
)

map_col <- sf::st_read(url_geojson, quiet = TRUE)

# Media posterior de cada sigma_j
sig2_cols <- paste0("sig2", seq_len(m))

sigma_post <- colMeans(
     sqrt(as.matrix(chain_2$THETA[, sig2_cols, drop = FALSE])),
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
          sigma_post = as.numeric(sigma_post),
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
map_col_sigma <- map_col %>%
     left_join(
          estadisticos_mapa,
          by = "nombre_key"
     )

# Revisar departamentos sin emparejar
map_col_sigma %>%
     filter(is.na(sigma_post)) %>%
     select(DPTO_CNMBR, nombre_key)

# Excluir San Andrés y Providencia para evitar espacio vacío a la izquierda
map_col_sigma_cont <- map_col_sigma %>%
     filter(
          !grepl(
               pattern     = "SAN ANDRES|PROVIDENCIA|SANTA CATALINA",
               x           = nombre_key,
               ignore.case = TRUE
          )
     )

# Límites del mapa
bbox <- sf::st_bbox(map_col_sigma_cont)

pad_x <- 0.35
pad_y <- 0.35

# Mapa
p_mapa <- ggplot(map_col_sigma_cont) +
     geom_sf(
          aes(fill = sigma_post),
          color     = "gray40",
          linewidth = 0.15
     ) +
     scale_fill_gradientn(
          colors   = c("#eff6ff", "#bfdbfe", "#60a5fa", "#2563eb", "#1e3a8a"),
          na.value = "gray90",
          name     = expression(hat(sigma)[j])
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
     filename = "matematicas_modelo_jerarquico_normal_heterocedastico_mapa_media_posterior_sigma.pdf",
     plot     = p_mapa,
     width    = 5.8,
     height   = 7.2,
     units    = "in"
)

dev.off()

# Fin --------------------------------------------------------------------------