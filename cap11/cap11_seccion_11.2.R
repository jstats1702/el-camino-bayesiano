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
     file = file.path("Examen_Saber_11_20251.txt"),
     sep = ";",
     stringsAsFactors = FALSE
)

dim(dat)

# Revisión inicial por departamento
table(dat$estu_depto_reside, useNA = "ifany")

round(
     100 * table(dat$estu_depto_reside, useNA = "ifany") / nrow(dat),
     1
)

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

# Muestreo proporcional por departamento ---------------------------------------

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
          sample(
               idx_depto,
               size = n_depto,
               replace = FALSE
          )
     )
}

# Base final -------------------------------------------------------------------

dat <- dat[idx_muestra, , drop = FALSE]

dat <- dat |>
     dplyr::select(
          punt_matematicas,
          estu_genero,
          punt_lectura_critica,
          estu_depto_reside
     ) |>
     dplyr::filter(
          !is.na(punt_matematicas),
          !is.na(estu_genero),
          !is.na(punt_lectura_critica),
          !is.na(estu_depto_reside)
     ) |>
     dplyr::mutate(
          estu_genero       = factor(estu_genero),
          estu_depto_reside = factor(estu_depto_reside)
     )

dat$estu_genero <- factor(dat$estu_genero)

# Dimensión base de datos
dim(dat)

# Número de estudiantes
(n <- nrow(dat))

# Número de departamentos
(m <- length(table(dat$estu_depto_reside)))

# Número de estudiantes por departamento
table(dat$estu_depto_reside)

# Resumen global de las variables
summary(dat$punt_matematicas)

summary(dat$punt_lectura_critica)

round(sd(dat$punt_matematicas), 3)

round(sd(dat$punt_lectura_critica), 3)

# Resumen por departamento -----------------------------------------------------

tab_depto <- dat |>
     dplyr::group_by(estu_depto_reside) |>
     dplyr::summarise(
          n              = dplyr::n(),
          media_mate     = mean(punt_matematicas),
          sd_mate        = sd(punt_matematicas),
          media_lectura  = mean(punt_lectura_critica),
          sd_lectura     = sd(punt_lectura_critica),
          corr_mate_lect = cor(punt_matematicas, punt_lectura_critica),
          .groups        = "drop"
     ) |>
     dplyr::arrange(dplyr::desc(media_mate))

tab_depto

summary(tab_depto$n)

summary(tab_depto$media_mate)

summary(tab_depto$media_lectura)

summary(tab_depto$corr_mate_lect)

# Histograma del puntaje de Matemáticas ----------------------------------------

h_matematicas <- hist(
     x      = dat$punt_matematicas,
     breaks = 15,
     plot   = FALSE
)

h_lectura <- hist(
     x      = dat$punt_lectura_critica,
     breaks = 15,
     plot   = FALSE
)

pdf(
     file      = "matematicas_modelo_efectos_mixtos_normal_hist_mate.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

hist(
     x      = dat$punt_matematicas,
     freq   = FALSE,
     breaks = 15,
     col    = "#d8b4fe",
     border = "white",
     xlab   = "Puntaje de Matemáticas",
     ylab   = "Densidad",
     ylim   = c(0, max(h_matematicas$density, h_lectura$density)),
     main   = ""
)

abline(
     v   = mean(dat$punt_matematicas, na.rm = TRUE),
     col = "black",
     lwd = 3
)

dev.off()

# Histograma del puntaje de Lectura Crítica ------------------------------------

pdf(
     file      = "matematicas_modelo_efectos_mixtos_normal_hist_lect.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

hist(
     x      = dat$punt_lectura_critica,
     freq   = FALSE,
     breaks = 15,
     col    = "#93c5fd",
     border = "white",
     xlab   = "Puntaje de Lectura Crítica",
     ylab   = "Densidad",
     ylim   = c(0, max(h_matematicas$density, h_lectura$density)),
     main   = ""
)

abline(
     v   = mean(dat$punt_lectura_critica, na.rm = TRUE),
     col = "black",
     lwd = 3
)

dev.off()

# Mapa puntaje de Matemáticas --------------------------------------------------

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
estadisticos_mapa <- dat %>%
     filter(
          !is.na(estu_depto_reside),
          !is.na(punt_matematicas)
     ) %>%
     group_by(estu_depto_reside) %>%
     summarise(
          yb      = mean(punt_matematicas),
          .groups = "drop"
     ) %>%
     mutate(
          nombre_key = estu_depto_reside %>%
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
          colors = c(
               "#f3e8ff",
               "#d8b4fe",
               "#c084fc",
               "#a855f7",
               "#7e22ce"
          ),
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

dev.off()

# Guardar el mapa
ggsave(
     filename = "matematicas_modelo_efectos_mixtos_normal_mapa_mate.pdf",
     plot     = p_mapa,
     width    = 5.8,
     height   = 7.2,
     units    = "in"
)

# Mapa puntaje de Lectura Crítica ----------------------------------------------

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
estadisticos_mapa <- dat %>%
     filter(
          !is.na(estu_depto_reside),
          !is.na(punt_lectura_critica)
     ) %>%
     group_by(estu_depto_reside) %>%
     summarise(
          yb      = mean(punt_lectura_critica),
          .groups = "drop"
     ) %>%
     mutate(
          nombre_key = estu_depto_reside %>%
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
          colors = c(
               "#dbeafe",
               "#93c5fd",
               "#60a5fa",
               "#2563eb",
               "#1e3a8a"
          ),
          na.value = "gray90",
          name     = expression(bar(x)[j])
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

dev.off()

# Guardar el mapa
ggsave(
     filename = "matematicas_modelo_efectos_mixtos_normal_mapa_lect.pdf",
     plot     = p_mapa,
     width    = 5.8,
     height   = 7.2,
     units    = "in"
)

# Puntaje de Matemáticas por departamento --------------------------------------

# Ordenar departamentos por promedio del puntaje de Matemáticas
media_depto <- tapply(
     X     = dat$punt_matematicas,
     INDEX = dat$estu_depto_reside,
     FUN   = mean,
     na.rm = TRUE
)

depto_orden <- names(sort(media_depto))

dat$estu_depto_reside <- factor(
     dat$estu_depto_reside,
     levels = depto_orden
)

pdf(
     file      = "matematicas_modelo_efectos_mixtos_normal_observaciones_depto_mate.pdf",
     width     = 5.8,
     height    = 7.2,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(4, 6, 1.4, 1.4),
     mgp   = c(2.2, 0.75, 0)
)

# Gráfico de observaciones por departamento
plot(
     x    = dat$punt_matematicas,
     y    = as.numeric(dat$estu_depto_reside),
     pch  = 16,
     cex  = 0.45,
     col  = adjustcolor("#7e22ce", alpha.f = 0.25),
     xlim = c(0, 100),
     ylim = c(0.5, length(depto_orden) + 0.5),
     xlab = "Puntaje de Matemáticas",
     ylab = "",
     axes = FALSE,
     main = "",
     bty  = "n"
)

# Líneas de referencia
abline(
     v   = seq(0, 100, by = 10),
     col = "gray90",
     lty = "dotted"
)

abline(
     h   = seq_along(depto_orden),
     col = "gray95",
     lty = "solid"
)

abline(
     v   = mean(dat$punt_matematicas, na.rm = TRUE),
     col = "black",
     lwd = 2
)

# Redibujar observaciones sobre la grilla
points(
     x   = dat$punt_matematicas,
     y   = as.numeric(dat$estu_depto_reside),
     pch = 16,
     cex = 0.45,
     col = adjustcolor("#7e22ce", alpha.f = 0.25)
)

# Promedios departamentales
points(
     x   = media_depto[depto_orden],
     y   = seq_along(depto_orden),
     pch = 19,
     cex = 0.8,
     col = "#7e22ce"
)

# Ejes
axis(
     side = 1,
     at   = seq(0, 100, by = 10)
)

axis(
     side     = 2,
     at       = seq_along(depto_orden),
     labels   = depto_orden,
     las      = 1,
     cex.axis = 0.65
)

dev.off()

# Puntaje de Lectura Crítica por departamento ----------------------------------

# Ordenar departamentos por promedio del puntaje de Lectura Crítica
media_depto <- tapply(
     X     = dat$punt_lectura_critica,
     INDEX = dat$estu_depto_reside,
     FUN   = mean,
     na.rm = TRUE
)

depto_orden <- names(sort(media_depto))

dat$estu_depto_reside <- factor(
     dat$estu_depto_reside,
     levels = depto_orden
)

pdf(
     file      = "matematicas_modelo_efectos_mixtos_normal_observaciones_depto_lect.pdf",
     width     = 5.8,
     height    = 7.2,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(4, 6, 1.4, 1.4),
     mgp   = c(2.2, 0.75, 0)
)

# Gráfico de observaciones por departamento
plot(
     x    = dat$punt_lectura_critica,
     y    = as.numeric(dat$estu_depto_reside),
     pch  = 16,
     cex  = 0.45,
     col  = adjustcolor("#1e3a8a", alpha.f = 0.25),
     xlim = c(0, 100),
     ylim = c(0.5, length(depto_orden) + 0.5),
     xlab = "Puntaje de Lectura Crítica",
     ylab = "",
     axes = FALSE,
     main = "",
     bty  = "n"
)

# Líneas de referencia
abline(
     v   = seq(0, 100, by = 10),
     col = "gray90",
     lty = "dotted"
)

abline(
     h   = seq_along(depto_orden),
     col = "gray95",
     lty = "solid"
)

abline(
     v   = mean(dat$punt_lectura_critica, na.rm = TRUE),
     col = "black",
     lwd = 2
)

# Redibujar observaciones sobre la grilla
points(
     x   = dat$punt_lectura_critica,
     y   = as.numeric(dat$estu_depto_reside),
     pch = 16,
     cex = 0.45,
     col = adjustcolor("#1e3a8a", alpha.f = 0.25)
)

# Promedios departamentales
points(
     x   = media_depto[depto_orden],
     y   = seq_along(depto_orden),
     pch = 19,
     cex = 0.8,
     col = "#1e3a8a"
)

# Ejes
axis(
     side = 1,
     at   = seq(0, 100, by = 10)
)

axis(
     side     = 2,
     at       = seq_along(depto_orden),
     labels   = depto_orden,
     las      = 1,
     cex.axis = 0.65
)

dev.off()

# Relación entre Matemáticas y sexo --------------------------------------------

tab_genero <- dat |>
     dplyr::group_by(estu_genero) |>
     dplyr::summarise(
          n             = dplyr::n(),
          media_mate    = mean(punt_matematicas),
          sd_mate       = sd(punt_matematicas),
          media_lectura = mean(punt_lectura_critica),
          sd_lectura    = sd(punt_lectura_critica),
          .groups       = "drop"
     )

tab_genero

# Puntaje de Matemáticas por sexo ----------------------------------------------

pdf(
     file      = "matematicas_modelo_efectos_mixtos_normal_boxplots_sexo_mate.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Diagrama de cajas por sexo
boxplot(
     punt_matematicas ~ estu_genero,
     data    = dat,
     col     = adjustcolor("gray", alpha.f = 0.45),
     border  = "gray30",
     outline = FALSE,
     boxwex  = 0.55,
     xlab    = "Sexo",
     ylab    = "Puntaje de Matemáticas",
     main    = "",
     ylim    = range(dat$punt_matematicas, na.rm = TRUE)
)

# Puntos individuales con dispersión horizontal
stripchart(
     punt_matematicas ~ estu_genero,
     data     = dat,
     vertical = TRUE,
     method   = "jitter",
     jitter   = 0.18,
     pch      = 16,
     cex      = 0.3,
     col      = adjustcolor("black", alpha.f = 0.12),
     add      = TRUE
)

# Promedios por sexo
media_genero <- tapply(
     X     = dat$punt_matematicas,
     INDEX = dat$estu_genero,
     FUN   = mean,
     na.rm = TRUE
)

points(
     x   = seq_along(media_genero),
     y   = media_genero,
     pch = 19,
     cex = 1.3,
     col = "black"
)

dev.off()

# Relación entre Matemáticas y Lectura Crítica ---------------------------------

# Correlación
cor_lect_mate <- cor(
     dat$punt_lectura_critica,
     dat$punt_matematicas
)

# Modelo lineal auxiliar
fit_lm <- lm(
     punt_matematicas ~ punt_lectura_critica,
     data = dat
)

summary(fit_lm)

coef(fit_lm)

pdf(
     file      = "matematicas_modelo_efectos_mixtos_normal_scatter_lect_mate.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Gráfico de dispersión
plot(
     x    = dat$punt_lectura_critica,
     y    = dat$punt_matematicas,
     pch  = 16,
     cex  = 0.65,
     col  = adjustcolor("black", alpha.f = 0.2),
     xlab = "Puntaje de Lectura Crítica",
     ylab = "Puntaje de Matemáticas",
     main = ""
)

# Recta de regresión lineal
abline(
     fit_lm,
     col = "black",
     lwd = 2
)

# Promedios marginales
abline(
     v   = mean(dat$punt_lectura_critica),
     col = "gray60",
     lwd = 2,
     lty = 3
)

abline(
     h   = mean(dat$punt_matematicas),
     col = "gray60",
     lwd = 2,
     lty = 3
)

dev.off()

# Tratamiento de datos para ajustar el modelo ----------------------------------

# Codificar sexo: 1 = hombre, 0 = mujer
dat$sexo <- ifelse(
     toupper(trimws(dat$estu_genero)) == "M",
     1,
     0
)

# Centrar Lectura Crítica en el valor de referencia de la prueba
dat$punt_lectura_critica <- dat$punt_lectura_critica - 50

# Extraer los IDs únicos de los departamentos e inicializar las estructuras
depto_ids <- sort(unique(dat$estu_depto_reside))
m         <- length(depto_ids)

Y <- vector("list", m)
X <- vector("list", m)
N <- integer(m)

# Recorrer cada departamento
for (j in seq_len(m)) {
     indices <- dat$estu_depto_reside == depto_ids[j]
     
     Y[[j]] <- as.matrix(dat$punt_matematicas[indices])
     N[j]   <- sum(indices)
     
     # Covariables dentro de cada departamento
     sexo_j <- dat$sexo[indices]
     
     lectura_j <- dat$punt_lectura_critica[indices]
     
     # Matriz de diseño con intercepto, sexo y Lectura Crítica centrada
     X[[j]] <- as.matrix(cbind(1, sexo_j, lectura_j))
}

# Matriz de diseño para efectos aleatorios
Z <- X

# Ajustar regresiones separadas por departamento -------------------------------

BETA_LS <- matrix(
     NA,
     nrow = m,
     ncol = 3
)

colnames(BETA_LS) <- c(
     "Intercepto",
     "Sexo",
     "Lectura"
)

S2_LS <- numeric(m)

for (j in seq_len(m)) {
     fit <- lm(Y[[j]] ~ -1 + X[[j]])
     
     BETA_LS[j, ] <- coef(fit)
     S2_LS[j]     <- summary(fit)$sigma^2
}

BETA_MLS <- colMeans(BETA_LS)

# Líneas de regresión por departamento -----------------------------------------

x_lectura <- unlist(
     lapply(
          X,
          function(x) x[, 3]
     )
)

y_mate <- unlist(Y)

xlim_lectura <- range(x_lectura, na.rm = TRUE)
ylim_mate    <- range(y_mate, na.rm = TRUE)

pdf(
     file      = "matematicas_modelo_efectos_mixtos_normal_lineas_mate_vs_lect.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = xlim_lectura,
     y    = ylim_mate,
     type = "n",
     xlab = "Lectura Crítica centrada",
     ylab = "Puntaje de Matemáticas",
     main = "",
     bty  = "n"
)

for (j in seq_len(m)) {
     abline(
          a   = BETA_LS[j, 1],
          b   = BETA_LS[j, 3],
          col = adjustcolor("gray", alpha.f = 0.65),
          lwd = 1
     )
}

abline(
     a   = BETA_MLS[1],
     b   = BETA_MLS[3],
     col = "black",
     lwd = 3
)

legend(
     x      = "topleft",
     legend = c("Departamentos", "Promedio OLS"),
     lty    = c(1, 1),
     lwd    = c(2, 2),
     col    = c(adjustcolor("gray", alpha.f = 0.65), "black"),
     bty    = "n"
)

box()

dev.off()

# Interceptos vs tamaño muestral -----------------------------------------------

pdf(
     file      = "matematicas_modelo_efectos_mixtos_normal_intercepto_vs_tamano.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = N,
     y    = BETA_LS[, 1],
     pch  = 16,
     cex  = 1.1,
     col  = adjustcolor("gray", alpha.f = 0.85),
     xlab = "Tamaño muestral departamental",
     ylab = "Intercepto",
     main = "",
     bty  = "n"
)

points(
     x   = N,
     y   = BETA_LS[, 1],
     pch = 16,
     cex = 1.1,
     col = adjustcolor("gray", alpha.f = 0.85)
)

abline(
     h   = 0,
     col = "gray50",
     lwd = 2,
     lty = 2
)

abline(
     h   = BETA_MLS[1],
     col = "black",
     lwd = 2.5
)

box()

dev.off()

# Efectos de sexo vs tamaño muestral -------------------------------------------

pdf(
     file      = "matematicas_modelo_efectos_mixtos_normal_sexo_vs_tamano.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = N,
     y    = BETA_LS[, 2],
     pch  = 16,
     cex  = 1.1,
     col  = adjustcolor("gray", alpha.f = 0.85),
     xlab = "Tamaño muestral departamental",
     ylab = "Efecto de sexo",
     main = "",
     bty  = "n"
)

points(
     x   = N,
     y   = BETA_LS[, 2],
     pch = 16,
     cex = 1.1,
     col = adjustcolor("gray", alpha.f = 0.85)
)

abline(
     h   = 0,
     col = "gray50",
     lwd = 2,
     lty = 2
)

abline(
     h   = BETA_MLS[2],
     col = "black",
     lwd = 2.5
)

box()

dev.off()

# Pendientes de Lectura Crítica vs tamaño muestral -----------------------------

pdf(
     file      = "matematicas_modelo_efectos_mixtos_normal_lect_vs_tamano.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mar = c(3.2, 3.2, 1.0, 1.0),
     mgp = c(1.9, 0.7, 0)
)

plot(
     x    = N,
     y    = BETA_LS[, 3],
     pch  = 16,
     cex  = 1.1,
     col  = adjustcolor("gray", alpha.f = 0.85),
     xlab = "Tamaño muestral departamental",
     ylab = "Pendiente de Lectura Crítica",
     main = "",
     bty  = "n"
)

points(
     x   = N,
     y   = BETA_LS[, 3],
     pch = 16,
     cex = 1.1,
     col = adjustcolor("gray", alpha.f = 0.85)
)

abline(
     h   = BETA_MLS[3],
     col = "black",
     lwd = 2.5
)

box()

dev.off()

# Muestreador de Gibbs ---------------------------------------------------------

sample_gamma <- function(yj, Xj, Zj, theta, Sigma, sigma2) {
     Vj <- solve(
          solve(Sigma) +
               t(Zj) %*% Zj / sigma2
     )
     
     mj <- Vj %*% (
          t(Zj) %*% (yj - Xj %*% theta)
     ) / sigma2
     
     as.matrix(
          c(
               mvtnorm::rmvnorm(
                    1,
                    mean  = mj,
                    sigma = Vj
               )
          )
     )
}

sample_theta <- function(y, X, Z, G, sigma2, mu0, Lambda0) {
     p <- length(mu0)
     
     precision <- solve(Lambda0)
     score     <- precision %*% mu0
     
     for (j in seq_along(y)) {
          yj <- y[[j]]
          Xj <- X[[j]]
          Zj <- Z[[j]]
          gj <- G[[j]]
          
          precision <- precision + (t(Xj) %*% Xj) / sigma2
          score <- score +
               (t(Xj) %*% (yj - Zj %*% gj)) / sigma2
     }
     
     V_theta <- solve(precision)
     m_theta <- V_theta %*% score
     
     as.matrix(
          c(
               mvtnorm::rmvnorm(
                    1,
                    mean  = c(m_theta),
                    sigma = V_theta
               )
          )
     )
}

sample_Sigma <- function(G, nu0, S0) {
     q  <- nrow(G[[1]])
     Sn <- matrix(0, q, q)
     
     for (j in seq_along(G)) {
          gj <- G[[j]]
          Sn <- Sn + gj %*% t(gj)
     }
     
     nun <- nu0 + length(G)
     Sn  <- solve(S0 + Sn)
     
     solve(
          stats::rWishart(
               n     = 1,
               df    = nun,
               Sigma = Sn
          )[, , 1]
     )
}

sample_sigma2 <- function(Y, X, Z, theta, G, eta0, sigma20) {
     SSR <- 0
     
     for (j in seq_along(Y)) {
          SSR <- SSR + sum(
               (
                    Y[[j]] -
                         X[[j]] %*% theta -
                         Z[[j]] %*% G[[j]]
               )^2
          )
     }
     
     shape <- 0.5 * (eta0 + sum(sapply(Y, length)))
     rate  <- 0.5 * (eta0 * sigma20 + SSR)
     
     1 / rgamma(
          1,
          shape = shape,
          rate  = rate
     )
}

log_likelihood <- function(Y, X, Z, G, theta, sigma2) {
     out <- 0
     
     for (j in seq_along(Y)) {
          muj <- X[[j]] %*% theta + Z[[j]] %*% G[[j]]
          
          out <- out + sum(
               dnorm(
                    Y[[j]],
                    mean = muj,
                    sd   = sqrt(sigma2),
                    log  = TRUE
               )
          )
     }
     
     out
}

gibbs_sampler <- function(
          Y,
          X,
          Z,
          n_iter,
          n_burn,
          n_skip,
          mu0,
          Lambda0,
          nu0,
          S0,
          eta0,
          sigma20
) {
     # Configuración
     m <- length(Y)
     p <- ncol(X[[1]])
     q <- ncol(Z[[1]])
     
     B <- floor((n_iter - n_burn) / n_skip)
     
     gamma_out  <- vector("list", B)
     Sigma_out  <- vector("list", B)
     theta_out  <- matrix(NA, nrow = B, ncol = p)
     sigma2_out <- numeric(B)
     loglik_out <- numeric(B)
     
     # Inicializar parámetros
     sigma2 <- 1 / rgamma(
          1,
          shape = eta0 / 2,
          rate  = (eta0 * sigma20) / 2
     )
     
     Sigma <- solve(
          stats::rWishart(
               n     = 1,
               df    = nu0,
               Sigma = solve(S0)
          )[, , 1]
     )
     
     theta <- as.matrix(
          c(
               mvtnorm::rmvnorm(
                    1,
                    mean  = mu0,
                    sigma = Lambda0
               )
          )
     )
     
     G <- lapply(
          1:m,
          function(j) {
               as.matrix(
                    c(
                         mvtnorm::rmvnorm(
                              1,
                              mean  = rep(0, q),
                              sigma = Sigma
                         )
                    )
               )
          }
     )
     
     for (b in 1:n_iter) {
          # Actualizar parámetros
          for (j in 1:m) {
               G[[j]] <- sample_gamma(
                    Y[[j]],
                    X[[j]],
                    Z[[j]],
                    theta,
                    Sigma,
                    sigma2
               )
          }
          
          theta  <- sample_theta(Y, X, Z, G, sigma2, mu0, Lambda0)
          Sigma  <- sample_Sigma(G, nu0, S0)
          sigma2 <- sample_sigma2(Y, X, Z, theta, G, eta0, sigma20)
          
          # Mensaje de progreso
          if (b %% ceiling(n_iter / 10) == 0) {
               cat(
                    "Progreso: ",
                    round(100 * b / n_iter),
                    "%\n",
                    sep = ""
               )
          }
          
          # Guardar muestras después de n_burn, aplicando adelgazamiento
          if (b > n_burn && (b - n_burn) %% n_skip == 0) {
               i <- (b - n_burn) / n_skip
               
               gamma_out[[i]] <- G
               Sigma_out[[i]] <- Sigma
               theta_out[i, ] <- theta
               sigma2_out[i]  <- sigma2
               loglik_out[i]  <- log_likelihood(
                    Y,
                    X,
                    Z,
                    G,
                    theta,
                    sigma2
               )
          }
     }
     
     return(
          list(
               gamma  = gamma_out,
               Sigma  = Sigma_out,
               theta  = theta_out,
               sigma2 = sigma2_out,
               loglik = loglik_out
          )
     )
}

# Hiperparámetros --------------------------------------------------------------

# mu0     : promedio de las estimaciones OLS en todos los grupos
# Lambda0 : matriz de covarianzas muestral de las estimaciones OLS en todos los grupos
# nu0     : p + 2
# S0      : Lambda0
# eta0    : 1
# sigma20 : promedio de la varianza OLS muestral dentro de los grupos

mu0     <- as.matrix(c(colMeans(BETA_LS)))
Lambda0 <- cov(BETA_LS)
nu0     <- ncol(X[[1]]) + 2
S0      <- cov(BETA_LS)
eta0    <- 2
sigma20 <- mean(S2_LS)

# Ajuste del modelo ------------------------------------------------------------

# Iteraciones
n_iter <- 110000
n_burn <- 10000
n_skip <- 10

set.seed(123)

# Ajuste del modelo
samples <- gibbs_sampler(
     Y       = Y,
     X       = X,
     Z       = Z,
     n_iter  = n_iter,
     n_burn  = n_burn,
     n_skip  = n_skip,
     mu0     = mu0,
     Lambda0 = Lambda0,
     nu0     = nu0,
     S0      = S0,
     eta0    = eta0,
     sigma20 = sigma20
)

save(
     samples,
     file = "matematicas_muestras_modelo_efectos_mixtos.RData"
)

load(
     file = "matematicas_muestras_modelo_efectos_mixtos.RData"
)

# Traza de la log-verosimilitud ------------------------------------------------

plot(
     x    = samples$loglik,
     type = "p",
     pch  = 16,
     cex  = 0.3,
     xlab = "Iteración",
     ylab = "Log-verosimilitud",
     main = "Log-verosimilitud"
)

dev.off()

# Resumen posterior para theta y sigma -----------------------------------------

post_summary_simple <- function(x) {
     c(
          media = mean(x),
          de    = sd(x),
          q025  = unname(quantile(x, 0.025)),
          q975  = unname(quantile(x, 0.975))
     )
}

# Componentes de theta
tab_theta <- t(
     apply(
          samples$theta,
          2,
          post_summary_simple
     )
)

# Desviación estándar residual sigma
sigma <- sqrt(samples$sigma2)

tab_sigma <- t(
     post_summary_simple(sigma)
)

rownames(tab_sigma) <- "sigma"

# Resultados
round(tab_theta, 3)
round(tab_sigma, 3)

# Media posterior de Sigma -----------------------------------------------------

param_names <- c(
     "Intercepto",
     "Sexo",
     "Lectura Crítica"
)

Sigma_mean <- Reduce(
     f = "+",
     x = samples$Sigma
) / length(samples$Sigma)

rownames(Sigma_mean) <- param_names
colnames(Sigma_mean) <- param_names

round(Sigma_mean, 3)

# Media posterior de la matriz de correlación ----------------------------------

cor_samples <- lapply(
     X   = samples$Sigma,
     FUN = stats::cov2cor
)

cor_mean <- Reduce(
     f = "+",
     x = cor_samples
) / length(cor_samples)

rownames(cor_mean) <- param_names
colnames(cor_mean) <- param_names

round(cor_mean, 3)

# Inferencia posterior para gamma_j --------------------------------------------

B <- nrow(samples$theta)
m <- length(samples$gamma[[1]])
q <- ncol(samples$theta)

# Arreglo: iteración x departamento x parámetro
gamma_array <- array(
     NA,
     dim = c(B, m, q),
     dimnames = list(
          iter  = seq_len(B),
          depto = depto_ids,
          param = param_names
     )
)

for (b in seq_len(B)) {
     for (j in seq_len(m)) {
          gamma_array[b, j, ] <- as.numeric(
               samples$gamma[[b]][[j]]
          )
     }
}

# Resumen posterior para gamma_j -----------------------------------------------

post_summary_gamma <- function(x) {
     c(
          media = mean(x),
          q025  = unname(quantile(x, 0.025)),
          q975  = unname(quantile(x, 0.975))
     )
}

tab_gamma <- data.frame()

for (k in seq_len(q)) {
     aux <- t(
          apply(
               gamma_array[, , k],
               2,
               post_summary_gamma
          )
     )
     
     aux <- data.frame(
          departamento = rownames(aux),
          parametro    = param_names[k],
          media        = aux[, "media"],
          q025         = aux[, "q025"],
          q975         = aux[, "q975"],
          row.names    = NULL
     )
     
     tab_gamma <- rbind(tab_gamma, aux)
}

tab_gamma$clasificacion <- ifelse(
     tab_gamma$q025 > 0,
     "positivo",
     ifelse(
          tab_gamma$q975 < 0,
          "negativo",
          "contiene_cero"
     )
)

tab_gamma$color <- ifelse(
     tab_gamma$clasificacion == "positivo",
     "forestgreen",
     ifelse(
          tab_gamma$clasificacion == "negativo",
          "red3",
          "black"
     )
)

# Función para graficar intervalos por parámetro -------------------------------

plot_gamma_departamental <- function(tabla, parametro, file) {
     tabla_k <- tabla[
          tabla$parametro == parametro,
     ]
     
     tabla_k <- tabla_k[
          order(tabla_k$media),
     ]
     
     pdf(
          file   = file,
          width  = 6.8,
          height = 7.6
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(4, 6, 1.2, 1.2),
          mgp   = c(2.2, 0.75, 0)
     )
     
     plot(
          x    = tabla_k$media,
          y    = seq_len(nrow(tabla_k)),
          xlim = range(tabla_k[, c("q025", "q975")]),
          ylim = c(0.5, nrow(tabla_k) + 0.5),
          pch  = 16,
          cex  = 0.9,
          col  = tabla_k$color,
          axes = FALSE,
          xlab = paste0("Desviación departamental: ", parametro),
          ylab = "",
          main = "",
          bty  = "n"
     )
     
     grid(
          nx  = NULL,
          ny  = NA,
          col = "gray90",
          lty = "dotted"
     )
     
     abline(
          v   = 0,
          col = "gray50",
          lwd = 2,
          lty = 2
     )
     
     segments(
          x0  = tabla_k$q025,
          x1  = tabla_k$q975,
          y0  = seq_len(nrow(tabla_k)),
          y1  = seq_len(nrow(tabla_k)),
          col = tabla_k$color,
          lwd = 2
     )
     
     points(
          x   = tabla_k$media,
          y   = seq_len(nrow(tabla_k)),
          pch = 16,
          cex = 0.9,
          col = tabla_k$color
     )
     
     axis(
          side = 1
     )
     
     axis(
          side     = 2,
          at       = seq_len(nrow(tabla_k)),
          labels   = tabla_k$departamento,
          las      = 1,
          cex.axis = 0.62
     )
     
     dev.off()
}

# Gráficos separados para intercepto, sexo y Lectura Crítica -------------------

plot_gamma_departamental(
     tabla     = tab_gamma,
     parametro = "Intercepto",
     file      = "matematicas_modelo_efectos_mixtos_normal_gamma_intercepto.pdf"
)

plot_gamma_departamental(
     tabla     = tab_gamma,
     parametro = "Sexo",
     file      = "matematicas_modelo_efectos_mixtos_normal_gamma_sexo.pdf"
)

plot_gamma_departamental(
     tabla     = tab_gamma,
     parametro = "Lectura Crítica",
     file      = "matematicas_modelo_efectos_mixtos_normal_gamma_lectura.pdf"
)

# Arreglo para almacenar beta_j = theta + gamma_j ------------------------------

# Dimensiones
m <- length(X)
p <- ncol(X[[1]])
B <- nrow(samples$theta)

beta <- array(
     NA,
     dim = c(m, p, B)
)

for (b in 1:B) {
     for (j in 1:m) {
          beta[j, , b] <- samples$theta[b, ] +
               as.numeric(samples$gamma[[b]][[j]])
     }
}

# Gráfico: Matemáticas frente a Lectura Crítica --------------------------------

# Media posterior de beta_j
mean_beta <- apply(
     beta,
     MARGIN = c(1, 2),
     FUN    = mean
)

# Colores según la desviación gamma_j de Lectura Crítica
gamma_lectura <- tab_gamma[
     tab_gamma$parametro == "Lectura Crítica",
]

gamma_lectura <- gamma_lectura[
     match(depto_ids, gamma_lectura$departamento),
]

col_depto <- gamma_lectura$color

# Visualización
pdf(
     file      = "matematicas_modelo_efectos_mixtos_normal_rectas_lectura.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = dat$punt_lectura_critica,
     y    = dat$punt_matematicas,
     type = "n",
     xlim = range(dat$punt_lectura_critica, na.rm = TRUE),
     ylim = range(dat$punt_matematicas, na.rm = TRUE),
     xlab = "Puntaje de Lectura Crítica",
     ylab = "Puntaje de Matemáticas",
     main = "",
     bty  = "n"
)

points(
     x   = dat$punt_lectura_critica,
     y   = dat$punt_matematicas,
     pch = 16,
     cex = 0.45,
     col = adjustcolor("gray", alpha.f = 0.35)
)

# Rectas departamentales para sexo = 0, es decir, mujeres
sexo_ref <- 0

for (j in 1:m) {
     abline(
          a   = mean_beta[j, 1] + mean_beta[j, 2] * sexo_ref,
          b   = mean_beta[j, 3],
          col = adjustcolor(col_depto[j], alpha.f = 0.80),
          lwd = 1
     )
}

# Recta global para sexo = 0, es decir, mujeres
theta_bar <- colMeans(samples$theta)

abline(
     a   = theta_bar[1] + theta_bar[2] * sexo_ref,
     b   = theta_bar[3],
     col = 4,
     lwd = 3
)

box()

dev.off()

# Fin --------------------------------------------------------------------------