# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Librerías
suppressMessages(
     suppressWarnings({
          library(Rcpp)
          library(igraph)
          library(network)
          library(coda)
          library(corrplot)
     })
)

Rcpp::sourceCpp("dist_functions.cpp")

# Funciones auxiliares ---------------------------------------------------------

# Calcula el número de actores a partir del número de díadas.
get_n <- function(m) {
     n <- (1 + sqrt(1 + 8 * m)) / 2
     
     if (abs(n - round(n)) > sqrt(.Machine$double.eps)) {
          stop("La longitud del vector de díadas no es válida.")
     }
     
     as.integer(round(n))
}

# Convierte una matriz de adyacencia en un vector ordenado por filas.
adjacency_to_dyads <- function(A) {
     n <- nrow(A)
     y <- numeric(n * (n - 1) / 2)
     k <- 1L
     
     for (i in seq_len(n - 1L)) {
          for (j in (i + 1L):n) {
               y[k] <- A[i, j]
               k    <- k + 1L
          }
     }
     
     y
}

# Convierte un vector de díadas en una matriz simétrica.
dyads_to_adjacency <- function(y, n = get_n(length(y))) {
     A <- matrix(
          data = 0,
          nrow = n,
          ncol = n
     )
     
     k <- 1L
     
     for (i in seq_len(n - 1L)) {
          for (j in (i + 1L):n) {
               A[i, j] <- y[k]
               A[j, i] <- y[k]
               k       <- k + 1L
          }
     }
     
     A
}

# Calcula las distancias entre todas las díadas en el orden del muestreador.
pairwise_distances <- function(U) {
     n <- nrow(U)
     d <- numeric(n * (n - 1) / 2)
     k <- 1L
     
     for (i in seq_len(n - 1L)) {
          for (j in (i + 1L):n) {
               d[k] <- sqrt(sum((U[i, ] - U[j, ])^2))
               
               k <- k + 1L
          }
     }
     
     d
}

# Genera valores iniciales para el muestreador.
initialize_distance_model <- function(y, n, K, hyperparameters) {
     U <- matrix(
          data = rnorm(n * K),
          nrow = n,
          ncol = K
     )
     
     U <- scale(
          x      = U,
          center = TRUE,
          scale  = FALSE
     )
     
     p_obs <- mean(y, na.rm = TRUE)
     p_obs <- min(max(p_obs, 0.01), 0.99)
     
     mean_distance <- mean(pairwise_distances(U))
     
     mu <- qlogis(p_obs) + mean_distance
     
     sigma2 <- mean(U^2)
     
     omega2 <- max(
          hyperparameters$b_omega /(hyperparameters$a_omega - 1),
          mu^2
     )
     
     list(
          mu     = mu,
          U      = U,
          sigma2 = sigma2,
          omega2 = omega2
     )
}

# Calcula la log-verosimilitud de las díadas observadas.
log_likelihood_observed <- function(y, mu, U) {
     observed <- !is.na(y)
     
     eta <- mu - pairwise_distances(U)
     
     sum(
          dbinom(
               x    = y[observed],
               size = 1,
               prob = plogis(eta[observed]),
               log  = TRUE
          )
     )
}

# Ajusta el modelo de espacio latente de distancia.
fit_distance_model <- function(
          y,
          K,
          B,
          burn,
          thin,
          hyperparameters,
          seed = 42,
          verbose = TRUE,
          progress_every = 10000L
) {
     n       <- get_n(length(y))
     n_keep  <- B
     B_total <- burn + thin * B
     
     initial <- initialize_distance_model(
          y               = y,
          n               = n,
          K               = K,
          hyperparameters = hyperparameters
     )
     
     mu     <- initial$mu
     U      <- initial$U
     sigma2 <- initial$sigma2
     omega2 <- initial$omega2
     
     missing_dyads <- is.na(y)
     missing_index <- as.numeric(missing_dyads)
     
     y_work <- y
     
     if (any(missing_dyads)) {
          y_work[missing_dyads] <- rbinom(
               n    = sum(missing_dyads),
               size = 1,
               prob = mean(y, na.rm = TRUE)
          )
     }
     
     mu_chain <- numeric(n_keep)
     
     U_chain <- matrix(
          data = NA_real_,
          nrow = n_keep,
          ncol = n * K
     )
     
     sigma2_chain <- numeric(n_keep)
     omega2_chain <- numeric(n_keep)
     loglik_chain <- numeric(n_keep)
     
     del2_U  <- 0.1
     n_U     <- 0L
     n_tun_U <- 100L
     
     del2_mu  <- 0.1
     n_mu     <- 0L
     n_tun_mu <- 100L
     
     keep <- 0L
     
     for (b in seq_len(B_total)) {
          if (any(missing_dyads)) {
               y_work <- sample_Y(
                    I          = n,
                    zeta       = mu,
                    U          = U,
                    na_indices = missing_index,
                    Yna        = y_work
               )
          }
          
          update_U <- sample_U(
               b       = b,
               n_tun_U = n_tun_U,
               del2_U  = del2_U,
               n_U     = n_U,
               n_burn  = burn,
               I       = n,
               K       = K,
               sigsq   = sigma2,
               zeta    = mu,
               U       = U,
               Y       = y_work
          )
          
          U       <- update_U$U
          del2_U  <- update_U$del2_U
          n_U     <- update_U$n_U
          n_tun_U <- update_U$n_tun_U
          
          update_mu <- sample_zeta(
               b          = b,
               n_tun_zeta = n_tun_mu,
               del2_zeta  = del2_mu,
               n_zeta     = n_mu,
               n_burn     = burn,
               I          = n,
               omesq      = omega2,
               zeta       = mu,
               U          = U,
               Y          = y_work
          )
          
          mu       <- update_mu$zeta
          del2_mu  <- update_mu$del2_zeta
          n_mu     <- update_mu$n_zeta
          n_tun_mu <- update_mu$n_tun_zeta
          
          sigma2 <- sample_sigsq(
               I     = n,
               K     = K,
               a_sig = hyperparameters$a_sigma,
               b_sig = hyperparameters$b_sigma,
               U     = U
          )
          
          omega2 <- sample_omesq(
               a_ome = hyperparameters$a_omega,
               b_ome = hyperparameters$b_omega,
               zeta  = mu
          )
          
          if (
               b > burn &&
               (b - burn) %% thin == 0L
          ) {
               keep <- keep + 1L
               
               mu_chain[keep]     <- mu
               U_chain[keep, ]    <- c(U)
               sigma2_chain[keep] <- sigma2
               omega2_chain[keep] <- omega2
               
               loglik_chain[keep] <- log_likelihood_observed(
                    y  = y,
                    mu = mu,
                    U  = U
               )
          }
          
          if (
               verbose &&
               (
                    b == 1L ||
                    b %% progress_every == 0L ||
                    b == B_total
               )
          ) {
               cat(
                    sprintf(
                         paste0(
                              "Iteración %d de %d | ",
                              "aceptación U = %.3f | ",
                              "aceptación mu = %.3f\n"
                         ),
                         b,
                         B_total,
                         n_U / (b * n),
                         n_mu / b
                    )
               )
          }
     }
     
     result <- list(
          y = y,
          settings = list(
               n                = n,
               K                = K,
               B                = B,
               burn             = burn,
               thin             = thin,
               total_iterations = B_total,
               seed             = seed
          ),
          hyperparameters = hyperparameters,
          mu_chain = mu_chain,
          U_chain = U_chain,
          sigma2_chain = sigma2_chain,
          omega2_chain = omega2_chain,
          loglik_chain = loglik_chain,
          acceptance = c(
               U  = n_U / (B_total * n),
               mu = n_mu / B_total
          ),
          proposal_variance = c(
               U  = del2_U,
               mu = del2_mu
          )
     )
     
     result
}

# Convierte la cadena de posiciones en un arreglo.
extract_position_array <- function(fit) {
     n      <- fit$settings$n
     K      <- fit$settings$K
     B_post <- nrow(fit$U_chain)
     
     U_array <- array(
          data = NA_real_,
          dim  = c(B_post, n, K)
     )
     
     for (b in seq_len(B_post)) {
          U_array[b, , ] <- matrix(
               data = fit$U_chain[b, ],
               nrow = n,
               ncol = K
          )
     }
     
     U_array
}

# Alinea una configuración con una referencia mediante Procrustes.
align_procrustes <- function(U, U_reference) {
     U_centered <- scale(
          x      = U,
          center = TRUE,
          scale  = TRUE
     )
     
     reference_centered <- scale(
          x      = U_reference,
          center = TRUE,
          scale  = TRUE
     )
     
     decomposition <- svd(
          t(U_centered) %*% reference_centered
     )
     
     Q <- decomposition$u %*% t(decomposition$v)
     
     U_centered %*% Q
}

# Alinea todas las configuraciones posteriores.
align_position_array <- function(U_array) {
     B_post <- dim(U_array)[1]
     n      <- dim(U_array)[2]
     K      <- dim(U_array)[3]
     
     U_reference <- scale(
          x      = U_array[1, , ],
          center = TRUE,
          scale  = TRUE
     )
     
     U_aligned <- array(
          data = NA_real_,
          dim  = c(B_post, n, K)
     )
     
     for (b in seq_len(B_post)) {
          U_aligned[b, , ] <- align_procrustes(
               U           = U_array[b, , ],
               U_reference = U_reference
          )
     }
     
     U_aligned
}

# Calcula diagnósticos para una cadena.
diagnostic_summary <- function(x) {
     ess <- as.numeric(
          coda::effectiveSize(
               coda::mcmc(x)
          )
     )
     
     c(
          media = mean(x),
          de    = sd(x),
          ess   = ess,
          mcse  = sd(x) / sqrt(ess)
     )
}

# Calcula un resumen posterior.
posterior_summary <- function(x) {
     c(
          media = mean(x),
          de    = sd(x),
          q025  = unname(quantile(x, 0.025)),
          q975  = unname(quantile(x, 0.975))
     )
}

# Calcula las matrices posteriores de probabilidades y distancias.
posterior_relational_matrices <- function(fit) {
     n      <- fit$settings$n
     K      <- fit$settings$K
     B_post <- length(fit$mu_chain)
     
     probability_vector <- interaction_probs0(
          I          = n,
          K          = K,
          B          = B_post,
          zeta_chain = fit$mu_chain,
          U_chain    = fit$U_chain
     )
     
     probability_matrix <- dyads_to_adjacency(
          y = probability_vector,
          n = n
     )
     
     distance_matrix <- matrix(
          data = 0,
          nrow = n,
          ncol = n
     )
     
     for (b in seq_len(B_post)) {
          U <- matrix(
               data = fit$U_chain[b, ],
               nrow = n,
               ncol = K
          )
          
          distance_matrix <- distance_matrix +
               as.matrix(dist(U)) / B_post
     }
     
     list(
          probability = probability_matrix,
          distance    = distance_matrix
     )
}

# Calcula estadísticos de una red.
network_statistics <- function(A) {
     g <- igraph::graph_from_adjacency_matrix(
          adjmatrix = A,
          mode      = "undirected",
          diag      = FALSE
     )
     
     c(
          densidad = igraph::edge_density(
               graph = g,
               loops = FALSE
          ),
          transitividad = igraph::transitivity(
               graph = g,
               type  = "global"
          ),
          asortatividad = igraph::assortativity_degree(
               graph    = g,
               directed = FALSE
          ),
          distancia = igraph::mean_distance(
               graph       = g,
               directed    = FALSE,
               unconnected = TRUE
          ),
          grado_promedio = mean(
               igraph::degree(g)
          ),
          desviacion_grado = sd(
               igraph::degree(g)
          )
     )
}

# Realiza las comprobaciones predictivas posteriores.
posterior_predictive_check <- function(
          fit,
          A_observed,
          B_ppc = 5000L,
          seed = NULL
) {
     if (!is.null(seed)) {
          set.seed(seed)
     }
     
     n      <- fit$settings$n
     K      <- fit$settings$K
     B_post <- length(fit$mu_chain)
     
     B_ppc <- min(
          as.integer(B_ppc),
          B_post
     )
     
     index <- sample(
          x       = seq_len(B_post),
          size    = B_ppc,
          replace = FALSE
     )
     
     replicated_statistics <- matrix(
          data = NA_real_,
          nrow = B_ppc,
          ncol = 6
     )
     
     colnames(replicated_statistics) <- names(
          network_statistics(A_observed)
     )
     
     for (r in seq_len(B_ppc)) {
          b <- index[r]
          
          U <- matrix(
               data = fit$U_chain[b, ],
               nrow = n,
               ncol = K
          )
          
          A_rep <- simulate_data(
               I    = n,
               zeta = fit$mu_chain[b],
               U    = U
          )
          
          replicated_statistics[r, ] <- network_statistics(
               A = A_rep
          )
     }
     
     observed_statistics <- network_statistics(
          A = A_observed
     )
     
     ppp <- vapply(
          X = seq_along(observed_statistics),
          FUN = function(k) {
               mean(
                    replicated_statistics[, k] <
                         observed_statistics[k],
                    na.rm = TRUE
               )
          },
          FUN.VALUE = numeric(1)
     )
     
     names(ppp) <- names(observed_statistics)
     
     list(
          observed   = observed_statistics,
          replicated = replicated_statistics,
          ppp        = ppp
     )
}

# Grafica una comprobación predictiva posterior.
plot_ppc <- function(
          replicated,
          observed,
          ppp,
          file,
          xlab
) {
     pdf(
          file = file,
          width = 5,
          height = 5,
          pointsize = 17
     )
     
     par(
          mfrow = c(1, 1),
          mar = c(3, 3, 1.4, 1.4),
          mgp = c(1.75, 0.75, 0)
     )
     
     hist(
          replicated,
          breaks = 20,
          probability = TRUE,
          col = "gray85",
          border = "white",
          xlab = xlab,
          ylab = "Densidad",
          main = ""
     )
     
     abline(
          v = observed,
          lwd = 2,
          lty = 2
     )
     
     legend(
          x = "topright",
          legend = paste0(
               "ppp = ",
               formatC(
                    ppp,
                    format = "f",
                    digits = 3
               )
          ),
          bty = "n",
          cex = 1.1
     )
     
     box()
     
     dev.off()
}

# Datos ------------------------------------------------------------------------

data("florentine", package = "ergm")

A <- network::as.matrix.network.adjacency(x = flomarriage)

family_names <- flomarriage %v% "vertex.names"

isolated <- which(rowSums(A) == 0)

A <- A[
     -isolated,
     -isolated,
     drop = FALSE
]

family_names <- family_names[-isolated]

rownames(A) <- family_names
colnames(A) <- family_names

n <- nrow(A)
y <- adjacency_to_dyads(A)

g <- igraph::graph_from_adjacency_matrix(
     adjmatrix = A,
     mode      = "undirected",
     diag      = FALSE
)

s       <- igraph::gsize(g)
m       <- n * (n - 1) / 2
density <- s / m

c(
     actores  = n,
     aristas  = s,
     diadas   = m,
     densidad = density
)

# Visualización grafo ----------------------------------------------------------

set.seed(1702)

graph_layout <- igraph::layout_with_kk(g)

pdf(
     file = "florencia_distancia_grafo.pdf",
     width = 5,
     height = 5,
     pointsize = 15
)

par(
     mar = c(1, 1, 1, 1)
)

plot(
     g,
     layout = graph_layout,
     vertex.label = family_names,
     vertex.label.color = "black",
     vertex.label.cex = 0.9,
     vertex.size = 12,
     vertex.color = "white",
     vertex.frame.color = "black",
     edge.color = adjustcolor("blue4", 1),
     edge.width = 0.5,
     edge.curved = 0,
     main = ""
)

dev.off()

# Visualización matriz de adyacencia -------------------------------------------

pdf(
     file      = "florencia_distancia_matriz.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mar = c(0, 0, 0, 0)
)

corrplot::corrplot(
     corr = A,
     is.corr = FALSE,
     method = "color",
     type = "full",
     col = c("white", "blue4"),
     col.lim = c(0, 1),
     tl.col = "black",
     tl.cex = 0.9,
     tl.srt = 90,
     addgrid.col = NA,
     cl.pos = "n",
     diag = TRUE,
     mar = c(0, 0, 0, 0),
     title = ""
)

rect(
     xleft = 0.5,
     ybottom = 0.5,
     xright = ncol(A) + 0.5,
     ytop = nrow(A) + 0.5,
     border = "black",
     lwd = 1,
     xpd = NA
)

dev.off()

# Ajuste del modelo ------------------------------------------------------------

K <- 2L

hyperparameters <- list(
     a_omega = 3,
     b_omega = 2 * 10^2,
     a_sigma = 3,
     b_sigma = 2 * 8 * pi * n^(2/K)
)

fit_distance <- fit_distance_model(
     y               = y,
     K               = K,
     B               = 10000,
     burn            = 10000,
     thin            = 100,
     hyperparameters = hyperparameters,
     seed            = 1702,
     verbose         = TRUE,
     progress_every  = 100000L
)

save(
     fit_distance,
     file = "familias_florentinas_modelo_distancia.RData"
)

load(file = "familias_florentinas_modelo_distancia.RData")

fit_distance$acceptance
fit_distance$proposal_variance

# Diagnósticos -----------------------------------------------------------------

diagnostic_table <- rbind(
     mu = diagnostic_summary(
          fit_distance$mu_chain
     ),
     sigma2 = diagnostic_summary(
          fit_distance$sigma2_chain
     ),
     omega2 = diagnostic_summary(
          fit_distance$omega2_chain
     ),
     logL = diagnostic_summary(
          fit_distance$loglik_chain
     )
)

round(
     diagnostic_table,
     digits = 3
)

# Trazas
par(
     mfrow = c(2, 2),
     mar   = c(3.5, 3.5, 1, 1),
     mgp   = c(2.2, 0.7, 0)
)

plot(
     fit_distance$mu_chain,
     type = "p",
     pch  = 16,
     cex  = 0.2,
     xlab = "Iteración",
     ylab = expression(mu)
)

plot(
     fit_distance$sigma2_chain,
     type = "p",
     pch  = 16,
     cex  = 0.2,
     xlab = "Iteración",
     ylab = expression(sigma^2)
)

plot(
     fit_distance$omega2_chain,
     type = "p",
     pch  = 16,
     cex  = 0.2,
     xlab = "Iteración",
     ylab = expression(omega^2)
)

plot(
     fit_distance$loglik_chain,
     type = "p",
     pch  = 16,
     cex  = 0.2,
     xlab = "Iteración",
     ylab = "Log-verosimilitud"
)

dev.off()

# Inferencia posterior ---------------------------------------------------------

posterior_table <- rbind(
     mu = posterior_summary(
          fit_distance$mu_chain
     ),
     sigma2 = posterior_summary(
          fit_distance$sigma2_chain
     ),
     omega2 = posterior_summary(
          fit_distance$omega2_chain
     )
)

round(
     posterior_table,
     digits = 3
)

# Posiciones latentes ----------------------------------------------------------

U_array <- extract_position_array(
     fit = fit_distance
)

U_aligned <- align_position_array(
     U_array = U_array
)

U_mean <- apply(
     X = U_aligned,
     MARGIN = c(2, 3),
     FUN = mean
)

B_post <- dim(U_aligned)[1]
n      <- dim(U_aligned)[2]

# Colores
U_reference <- scale(
     x = U_array[1, , ],
     center = TRUE,
     scale = FALSE
)

rr <- atan2(
     y = U_reference[, 2],
     x = U_reference[, 1]
)

rr <- (rr - min(rr)) / (max(rr) - min(rr))
gg <- 1 - rr

bb <- U_reference[, 1]^2 + U_reference[, 2]^2
bb <- (bb - min(bb)) / (max(bb) - min(bb))

position_colors <- rgb(
     red = rr,
     green = gg,
     blue = bb
)

position_colors_alpha <- adjustcolor(
     col = position_colors,
     alpha.f = 0.20
)

# Adelgazamiento para la visualización
nthin <- 10L

index_thin <- seq(
     from = nthin,
     to = B_post,
     by = nthin
)

x_range <- range(
     U_aligned[index_thin, , 1],
     U_mean[, 1]
)

y_range <- range(
     U_aligned[index_thin, , 2],
     U_mean[, 2]
)

# Posiciones posteriores y medias
pdf(
     file = "florencia_distancia_posiciones_posteriores.pdf",
     width = 5,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

plot(
     NA,
     xlim = range(x_range, y_range),
     ylim = range(x_range, y_range),
     xlab = "Dimensión 1",
     ylab = "Dimensión 2",
     main = ""
)

for (i in seq_len(n)) {
     points(
          x = U_aligned[index_thin, i, 1],
          y = U_aligned[index_thin, i, 2],
          pch = 15,
          cex = 0.5,
          col = position_colors_alpha[i]
     )
}

text(
     x = U_mean[, 1],
     y = U_mean[, 2],
     labels = seq_len(n),
     cex = 0.9,
     font = 2,
     col = "black"
)

dev.off()

# Posiciones posteriores promedio
pdf(
     file = "florencia_distancia_posiciones_promedio.pdf",
     width = 5,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

plot(
     x = U_mean[, 1],
     y = U_mean[, 2],
     xlim = range(x_range, y_range),
     ylim = range(x_range, y_range),
     xlab = "Dimensión 1",
     ylab = "Dimensión 2",
     pch = 19,
     cex = 1.2,
     col = "blue4",
     main = ""
)

text(
     x = U_mean[, 1],
     y = U_mean[, 2],
     labels = seq_len(n),
     pos = 3,
     cex = 0.9,
     col = "black"
)

dev.off()


# Probabilidades de interacción ------------------------------------------------

relational_matrices <- posterior_relational_matrices(
     fit = fit_distance
)

probability_matrix <- relational_matrices$probability
distance_matrix    <- relational_matrices$distance

rownames(probability_matrix) <- family_names
colnames(probability_matrix) <- family_names

rownames(distance_matrix) <- family_names
colnames(distance_matrix) <- family_names

family_order <- hclust(
     d = as.dist(distance_matrix)
)$order

A_ordered <- A[
     family_order,
     family_order
]

probability_ordered <- probability_matrix[
     family_order,
     family_order
]

# Matriz de adyacencia ordenada
pdf(
     file = "florencia_distancia_adyacencia.pdf",
     width = 5,
     height = 5,
     pointsize = 15
)

par(
     mar = c(0, 0, 0, 0)
)

corrplot::corrplot(
     corr = A_ordered,
     is.corr = FALSE,
     method = "color",
     type = "full",
     col = colorRampPalette(
          c("white", "yellow", "red")
     )(100),
     col.lim = c(0, 1),
     tl.col = "black",
     tl.cex = 0.9,
     tl.srt = 90,
     addgrid.col = NA,
     cl.pos = "n",
     diag = TRUE,
     mar = c(0, 0, 0, 0),
     title = ""
)

rect(
     xleft = 0.5,
     ybottom = 0.5,
     xright = ncol(A_ordered) + 0.5,
     ytop = nrow(A_ordered) + 0.5,
     border = "black",
     lwd = 1,
     xpd = NA
)

dev.off()

# Probabilidades de interacción ------------------------------------------------

pdf(
     file = "florencia_distancia_probabilidades.pdf",
     width = 5,
     height = 5,
     pointsize = 15
)

par(
     mar = c(0, 0, 0, 0)
)

corrplot::corrplot(
     corr = probability_ordered,
     is.corr = FALSE,
     method = "color",
     type = "full",
     col = colorRampPalette(
          c("white", "yellow", "red")
     )(100),
     col.lim = c(0, 1),
     tl.col = "black",
     tl.cex = 0.9,
     tl.srt = 90,
     addgrid.col = NA,
     cl.pos = "n",
     diag = TRUE,
     mar = c(0, 0, 0, 0),
     title = ""
)

rect(
     xleft = 0.5,
     ybottom = 0.5,
     xright = ncol(probability_ordered) + 0.5,
     ytop = nrow(probability_ordered) + 0.5,
     border = "black",
     lwd = 1,
     xpd = NA
)

dev.off()

# Díadas con mayores probabilidades posteriores
dyad_table <- data.frame(
     familia_1 = character(),
     familia_2 = character(),
     observada = integer(),
     probabilidad = numeric()
)

for (i in seq_len(n - 1L)) {
     for (j in (i + 1L):n) {
          dyad_table <- rbind(
               dyad_table,
               data.frame(
                    familia_1   = family_names[i],
                    familia_2   = family_names[j],
                    observada   = A[i, j],
                    probabilidad = probability_matrix[i, j]
               )
          )
     }
}

dyad_table <- dyad_table[
     order(
          dyad_table$probabilidad,
          decreasing = TRUE
     ),
]

head(
     dyad_table,
     n = 15
)

# Evaluación predictiva posterior ----------------------------------------------

ppc <- posterior_predictive_check(
     fit        = fit_distance,
     A_observed = A,
     B_ppc      = 10000,
     seed       = 2468
)

ppc_table <- data.frame(
     estadistico = names(ppc$observed),
     observado   = as.numeric(ppc$observed),
     ppp         = as.numeric(ppc$ppp)
)

ppc_table

plot_ppc(
     replicated = ppc$replicated[, "grado_promedio"],
     observed   = ppc$observed["grado_promedio"],
     ppp        = ppc$ppp["grado_promedio"],
     file       = "florencia_distancia_ppc_grado_promedio.pdf",
     xlab       = "Grado promedio"
)

plot_ppc(
     replicated = ppc$replicated[, "desviacion_grado"],
     observed   = ppc$observed["desviacion_grado"],
     ppp        = ppc$ppp["desviacion_grado"],
     file       = "florencia_distancia_ppc_desviacion_grado.pdf",
     xlab       = "DE del grado"
)

plot_ppc(
     replicated = ppc$replicated[, "densidad"],
     observed   = ppc$observed["densidad"],
     ppp        = ppc$ppp["densidad"],
     file       = "florencia_distancia_ppc_densidad.pdf",
     xlab       = "Densidad"
)

plot_ppc(
     replicated = ppc$replicated[, "transitividad"],
     observed   = ppc$observed["transitividad"],
     ppp        = ppc$ppp["transitividad"],
     file       = "florencia_distancia_ppc_transitividad.pdf",
     xlab       = "Transitividad"
)

plot_ppc(
     replicated = ppc$replicated[, "asortatividad"],
     observed   = ppc$observed["asortatividad"],
     ppp        = ppc$ppp["asortatividad"],
     file       = "florencia_distancia_ppc_asortatividad.pdf",
     xlab       = "Asortatividad por grado"
)

plot_ppc(
     replicated = ppc$replicated[, "distancia"],
     observed   = ppc$observed["distancia"],
     ppp        = ppc$ppp["distancia"],
     file       = "florencia_distancia_ppc_distancia.pdf",
     xlab       = "Distancia promedio"
)

# Fin --------------------------------------------------------------------------