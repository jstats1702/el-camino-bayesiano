# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Librerías
suppressMessages(
     suppressWarnings({
          library(igraph)
          library(sand)
          library(truncnorm)
          library(coda)
          library(cluster)
     })
)

# Datos ------------------------------------------------------------------------

data("elist.lazega", package = "sand")

g <- graph_from_data_frame(
     d = elist.lazega,
     directed = FALSE
)

n <- vcount(g)
s <- ecount(g)
m <- choose(n, 2)

Y <- as.matrix(
     as_adjacency_matrix(
          graph = g,
          sparse = FALSE
     )
)

node_labels <- V(g)$name

c(
     Actors = n,
     Edges = s,
     Dyads = m,
     Density = edge_density(g, loops = FALSE)
)

# Funciones auxiliares ---------------------------------------------------------

# Muestrear las variables auxiliares
sample_z <- function(y, mu, delta, z, pairs) {
     eta <- mu + delta[pairs[, 1]] + delta[pairs[, 2]]
     y_upper <- y[cbind(pairs[, 1], pairs[, 2])]
     
     z_upper <- numeric(nrow(pairs))
     
     index_one <- y_upper == 1
     index_zero <- y_upper == 0
     
     z_upper[index_one] <- truncnorm::rtruncnorm(
          n = sum(index_one),
          a = 0,
          b = Inf,
          mean = eta[index_one],
          sd = 1
     )
     
     z_upper[index_zero] <- truncnorm::rtruncnorm(
          n = sum(index_zero),
          a = -Inf,
          b = 0,
          mean = eta[index_zero],
          sd = 1
     )
     
     z[cbind(pairs[, 1], pairs[, 2])] <- z_upper
     z[cbind(pairs[, 2], pairs[, 1])] <- z_upper
     
     z
}

# Muestrear el efecto global
sample_mu <- function(z, delta, sigma2, pairs) {
     variance_mu <- 1 / (1 / sigma2 + nrow(pairs))
     
     mean_mu <- variance_mu * sum(
          z[cbind(pairs[, 1], pairs[, 2])] -
               delta[pairs[, 1]] -
               delta[pairs[, 2]]
     )
     
     rnorm(
          n = 1,
          mean = mean_mu,
          sd = sqrt(variance_mu)
     )
}

# Muestrear los efectos de sociabilidad
sample_delta <- function(z, mu, tau2, delta) {
     n <- length(delta)
     variance_delta <- 1 / (1 / tau2 + n - 1)
     
     for (i in seq_len(n)) {
          index_other <- setdiff(seq_len(n), i)
          
          mean_delta <- variance_delta * sum(
               z[i, index_other] -
                    mu -
                    delta[index_other]
          )
          
          delta[i] <- rnorm(
               n = 1,
               mean = mean_delta,
               sd = sqrt(variance_delta)
          )
     }
     
     delta - mean(delta)
}

# Muestrear la varianza del efecto global
sample_sigma2 <- function(mu, a_sigma, b_sigma) {
     1 / rgamma(
          n = 1,
          shape = a_sigma + 0.5,
          rate = b_sigma + 0.5 * mu^2
     )
}

# Muestrear la varianza de los efectos de sociabilidad
sample_tau2 <- function(delta, a_tau, b_tau) {
     1 / rgamma(
          n = 1,
          shape = a_tau + length(delta) / 2,
          rate = b_tau + 0.5 * sum(delta^2)
     )
}

# Calcular la log-verosimilitud
log_likelihood_sociality <- function(y, mu, delta, pairs) {
     eta <- mu + delta[pairs[, 1]] + delta[pairs[, 2]]
     
     probabilities <- pnorm(eta)
     probabilities <- pmin(pmax(probabilities, 1e-12), 1 - 1e-12)
     
     y_upper <- y[cbind(pairs[, 1], pairs[, 2])]
     
     sum(
          dbinom(
               x = y_upper,
               size = 1,
               prob = probabilities,
               log = TRUE
          )
     )
}

# Ajustar el modelo mediante Gibbs
gibbs_sampler_sociality <- function(
          y,
          B,
          burn,
          thin,
          a_sigma,
          b_sigma,
          a_tau,
          b_tau,
          progress_every
) {
     n <- nrow(y)
     pairs <- which(upper.tri(y), arr.ind = TRUE)
     total_iterations <- burn + B * thin
     
     mu     <- qnorm(mean(y[upper.tri(y)]))
     delta  <- rep(0, n)
     sigma2 <- 1
     tau2   <- 1
     z      <- matrix(0, nrow = n, ncol = n)
     
     samples <- list(
          mu = numeric(B),
          delta = matrix(0, nrow = B, ncol = n),
          sigma2 = numeric(B),
          tau2 = numeric(B),
          log_likelihood = numeric(B)
     )
     
     stored <- 0L
     
     for (iteration in seq_len(total_iterations)) {
          z <- sample_z(
               y = y,
               mu = mu,
               delta = delta,
               z = z,
               pairs = pairs
          )
          
          mu <- sample_mu(
               z = z,
               delta = delta,
               sigma2 = sigma2,
               pairs = pairs
          )
          
          delta <- sample_delta(
               z = z,
               mu = mu,
               tau2 = tau2,
               delta = delta
          )
          
          sigma2 <- sample_sigma2(
               mu = mu,
               a_sigma = a_sigma,
               b_sigma = b_sigma
          )
          
          tau2 <- sample_tau2(
               delta = delta,
               a_tau = a_tau,
               b_tau = b_tau
          )
          
          if (
               iteration > burn &&
               (iteration - burn) %% thin == 0
          ) {
               stored <- stored + 1L
               
               samples$mu[stored]      <- mu
               samples$delta[stored, ] <- delta
               samples$sigma2[stored]  <- sigma2
               samples$tau2[stored]    <- tau2
               
               samples$log_likelihood[stored] <- log_likelihood_sociality(
                    y = y,
                    mu = mu,
                    delta = delta,
                    pairs = pairs
               )
          }
          
          if (
               progress_every > 0L &&
               iteration %% progress_every == 0L
          ) {
               cat(
                    sprintf(
                         "Iteración %d de %d\n",
                         iteration,
                         total_iterations
                    )
               )
          }
     }
     
     samples
}

# Resumir una muestra posterior
posterior_summary <- function(x) {
     c(
          Mean = mean(x),
          SD = sd(x),
          Lower = unname(quantile(x, 0.025)),
          Upper = unname(quantile(x, 0.975))
     )
}

# Graficar una matriz
plot_matrix <- function(
          matrix_data,
          title,
          fill_title,
          low = "white",
          high = "black"
) {
     matrix_df <- as.data.frame(as.table(matrix_data))
     names(matrix_df) <- c("Actor_i", "Actor_j", "Value")
     
     matrix_df$Actor_i <- as.numeric(matrix_df$Actor_i)
     matrix_df$Actor_j <- as.numeric(matrix_df$Actor_j)
     
     ggplot(
          matrix_df,
          aes(
               x = Actor_i,
               y = Actor_j,
               fill = Value
          )
     ) +
          geom_tile() +
          scale_fill_gradient(
               low = low,
               high = high
          ) +
          scale_y_reverse() +
          coord_equal() +
          labs(
               title = title,
               x = "Actor $i$",
               y = "Actor $j$",
               fill = fill_title
          ) +
          theme_minimal(base_size = 11) +
          theme(
               panel.grid = element_blank(),
               plot.title = element_text(hjust = 0.5)
          )
}

# Calcular estadísticas de una red
network_statistics <- function(graph) {
     values <- c(
          Density = edge_density(
               graph = graph,
               loops = FALSE
          ),
          Transitivity = transitivity(
               graph = graph,
               type = "global"
          ),
          Assortativity = assortativity_degree(
               graph = graph,
               directed = FALSE
          ),
          AveragePath = mean_distance(
               graph = graph,
               directed = FALSE,
               unconnected = TRUE
          ),
          AverageDegree = mean(degree(graph)),
          SDDegree = sd(degree(graph))
     )
     
     values[!is.finite(values)] <- NA_real_
     
     values
}

# Ajuste del modelo ------------------------------------------------------------

a_sigma <- 2
b_sigma <- 1 / 3
a_tau   <- 2
b_tau   <- 1 / 3

B    <- 10000
burn <- 10000
thin <- 10

fit_file <- "lazega_sociabilidad_mcmc.rds"

set.seed(1702)

samples <- gibbs_sampler_sociality(
     y = Y,
     B = B,
     burn = burn,
     thin = thin,
     a_sigma = a_sigma,
     b_sigma = b_sigma,
     a_tau = a_tau,
     b_tau = b_tau,
     progress_every = 10000
)

saveRDS(
     object = samples,
     file = fit_file
)

samples <- readRDS(fit_file)

# Diagnósticos de convergencia -------------------------------------------------

global_draws <- cbind(
     mu = samples$mu,
     sigma2 = samples$sigma2,
     tau2 = samples$tau2,
     log_likelihood = samples$log_likelihood
)

global_effective_size <- coda::effectiveSize(
     coda::as.mcmc(global_draws)
)

global_mcse <- apply(
     global_draws,
     2,
     sd
) / sqrt(global_effective_size)

delta_effective_size <- coda::effectiveSize(
     coda::as.mcmc(samples$delta)
)

delta_mcse <- apply(
     samples$delta,
     2,
     sd
) / sqrt(delta_effective_size)

diagnostic_summary <- data.frame(
     Parameter = c(
          "mu",
          "delta",
          "sigma^2",
          "tau^2",
          "log L"
     ),
     Mean = c(
          mean(samples$mu),
          mean(colMeans(samples$delta)),
          mean(samples$sigma2),
          mean(samples$tau2),
          mean(samples$log_likelihood)
     ),
     SD = c(
          sd(samples$mu),
          mean(apply(samples$delta, 2, sd)),
          sd(samples$sigma2),
          sd(samples$tau2),
          sd(samples$log_likelihood)
     ),
     ESS = c(
          global_effective_size["mu"],
          mean(delta_effective_size),
          global_effective_size["sigma2"],
          global_effective_size["tau2"],
          global_effective_size["log_likelihood"]
     ),
     MCSE = c(
          global_mcse["mu"],
          mean(delta_mcse),
          global_mcse["sigma2"],
          global_mcse["tau2"],
          global_mcse["log_likelihood"]
     ),
     row.names = NULL
)

diagnostic_summary

# Inferencia sobre los parámetros globales -------------------------------------

global_summary <- rbind(
     posterior_summary(samples$mu),
     posterior_summary(samples$sigma2),
     posterior_summary(samples$tau2)
)

global_summary <- data.frame(
     Parameter = c(
          "mu",
          "sigma^2",
          "tau^2"
     ),
     global_summary,
     row.names = NULL,
     check.names = FALSE
)

global_summary

# Inferencia sobre los efectos de sociabilidad ---------------------------------

delta_mean <- colMeans(samples$delta)

delta_interval <- apply(
     samples$delta,
     2,
     quantile,
     probs = c(0.025, 0.975)
)

delta_df <- data.frame(
     Actor = node_labels,
     Mean = delta_mean,
     Lower = delta_interval[1, ],
     Upper = delta_interval[2, ]
)

delta_df$Classification <- ifelse(
     delta_df$Upper < 0,
     "Inferior a cero",
     ifelse(
          delta_df$Lower > 0,
          "Superior a cero",
          "Incluye el cero"
     )
)

delta_df$Classification <- factor(
     delta_df$Classification,
     levels = c(
          "Inferior a cero",
          "Incluye el cero",
          "Superior a cero"
     )
)

delta_df <- delta_df[
     order(delta_df$Mean),
]

delta_df$Order <- seq_len(nrow(delta_df))

label_offset <- 0.04 * diff(
     range(
          c(
               delta_df$Lower,
               delta_df$Upper
          )
     )
)

classification_colors <- c(
     "Inferior a cero" = "firebrick",
     "Incluye el cero" = "gray55",
     "Superior a cero" = "darkgreen"
)

# Visualización
point_colors <- unname(
     classification_colors[
          as.character(delta_df$Classification)
     ]
)

plot_limits <- range(
     c(
          delta_df$Lower,
          delta_df$Upper + label_offset
     )
)

pdf(
     file = "lazega_sociabilidad_efectos_base.pdf",
     width = 7,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(1.4, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

plot(
     x = delta_df$Order,
     y = delta_df$Mean,
     type = "n",
     xlim = c(
          0.5,
          nrow(delta_df) + 0.5
     ),
     ylim = plot_limits,
     xaxt = "n",
     xlab = "",
     ylab = expression(delta[i]),
     bty = "l"
)

abline(
     h = 0,
     lty = 2
)

arrows(
     x0 = delta_df$Order,
     y0 = delta_df$Lower,
     x1 = delta_df$Order,
     y1 = delta_df$Upper,
     angle = 90,
     code = 3,
     length = 0,
     lwd = 1.2,
     col = point_colors
)

points(
     x = delta_df$Order,
     y = delta_df$Mean,
     pch = 16,
     cex = 1,
     col = point_colors
)

delta_df$Actor <- sub(
     pattern = "^V",
     replacement = "",
     x = delta_df$Actor
)

text(
     x = delta_df$Order,
     y = delta_df$Upper + label_offset,
     labels = delta_df$Actor,
     srt = 90,
     adj = c(0, 0.5),
     cex = 0.7
)

legend(
     x = "bottomright",
     legend = names(classification_colors),
     col = classification_colors,
     pch = 16,
     pt.cex = 0.8,
     horiz = FALSE,
     bty = "n",
     xpd = NA,
     cex = 1
)

dev.off()

# Probabilidades posteriores y grados esperados --------------------------------

n_samples <- length(samples$mu)

theta_mean <- matrix(
     0,
     nrow = n,
     ncol = n
)

expected_degree_draws <- matrix(
     0,
     nrow = n_samples,
     ncol = n
)

for (b in seq_len(n_samples)) {
     theta_b <- pnorm(
          samples$mu[b] +
               outer(
                    samples$delta[b, ],
                    samples$delta[b, ],
                    "+"
               )
     )
     
     diag(theta_b) <- 0
     
     theta_mean <- theta_mean + theta_b / n_samples
     expected_degree_draws[b, ] <- rowSums(theta_b)
}

expected_degree_mean <- colMeans(expected_degree_draws)

expected_degree_interval <- apply(
     expected_degree_draws,
     2,
     quantile,
     probs = c(0.025, 0.975)
)

degree_df <- data.frame(
     Actor = node_labels,
     Observed = degree(g),
     Mean = expected_degree_mean,
     Lower = expected_degree_interval[1, ],
     Upper = expected_degree_interval[2, ]
)

degree_df <- degree_df[
     order(degree_df$Observed),
]

degree_df$Order <- seq_len(nrow(degree_df))

# Visualización
plot_limits <- range(
     c(
          degree_df$Observed,
          degree_df$Lower,
          degree_df$Upper
     )
)

pdf(
     file = "lazega_sociabilidad_grados.pdf",
     width = 7,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

plot(
     x = degree_df$Order,
     y = degree_df$Mean,
     type = "n",
     xlim = c(
          0.5,
          nrow(degree_df) + 0.5
     ),
     ylim = plot_limits,
     xaxt = "n",
     xlab = "Actor",
     ylab = "Grado",
     bty = "l"
)

arrows(
     x0 = degree_df$Order,
     y0 = degree_df$Lower,
     x1 = degree_df$Order,
     y1 = degree_df$Upper,
     angle = 90,
     code = 3,
     length = 0,
     lwd = 1
)

points(
     x = degree_df$Order,
     y = degree_df$Mean,
     pch = 16,
     cex = 1
)

points(
     x = degree_df$Order,
     y = degree_df$Observed,
     pch = 4,
     col = 2,
     cex = 1,
     lwd = 1.2
)

axis(
     side = 1,
     at = degree_df$Order,
     labels = sub(
          pattern = "^V",
          replacement = "",
          x = degree_df$Actor
     ),
     las = 2,
     cex.axis = 0.7
)

legend(
     x = "topleft",
     legend = c(
          "Media posterior",
          "Grado observado"
     ),
     col = c(1, 2),
     pch = c(16, 4),
     pt.cex = c(1, 1),
     bty = "n",
     cex = 1
)

dev.off()

# Probabilidades de interacción ------------------------------------------------

probability_order <- order(delta_mean)

Y_ordered <- Y[
     probability_order,
     probability_order
]

theta_ordered <- theta_mean[
     probability_order,
     probability_order
]

actor_ordered <- sub(
     pattern = "^V",
     replacement = "",
     x = node_labels[probability_order]
)

dimnames(Y_ordered) <- list(
     actor_ordered,
     actor_ordered
)

dimnames(theta_ordered) <- list(
     actor_ordered,
     actor_ordered
)

probability_colors <- colorRampPalette(
     c(
          "white",
          "yellow",
          "red"
     )
)(200)

# Matriz de adyacencia
pdf(
     file = "lazega_sociabilidad_adyacencia.pdf",
     width = 6,
     height = 5.8,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(1, 1, 1, 1)
)

corrplot::corrplot(
     corr = Y_ordered,
     is.corr = FALSE,
     method = "color",
     type = "full",
     order = "original",
     col = c(
          "white",
          "red"
     ),
     col.lim = c(0, 1),
     tl.col = "black",
     tl.cex = 0.7,
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
     xright = ncol(Y_ordered) + 0.5,
     ytop = nrow(Y_ordered) + 0.5,
     border = "black",
     lwd = 1,
     xpd = NA
)

dev.off()

# Probabilidades posteriores de interacción
pdf(
     file = "lazega_sociabilidad_probabilidades.pdf",
     width = 6,
     height = 5.8,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(1, 1, 1, 1)
)

corrplot::corrplot(
     corr = theta_ordered,
     is.corr = FALSE,
     method = "color",
     type = "full",
     col = colorRampPalette(
          c(
               "white",
               "yellow",
               "red"
          )
     )(100),
     col.lim = c(0, 1),
     tl.col = "black",
     tl.cex = 0.7,
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
     xright = ncol(theta_ordered) + 0.5,
     ytop = nrow(theta_ordered) + 0.5,
     border = "black",
     lwd = 1,
     xpd = NA
)

dev.off()

# Agrupamiento por niveles de sociabilidad -------------------------------------

n_samples <- nrow(samples$delta)
n <- ncol(samples$delta)

delta_mean <- colMeans(samples$delta)

delta_interval <- apply(
     samples$delta,
     2,
     quantile,
     probs = c(0.025, 0.975)
)

actor_labels <- sub(
     pattern = "^V",
     replacement = "",
     x = node_labels
)

# Reetiqueta una partición con enteros consecutivos
relabel_partition <- function(xi) {
     as.integer(
          match(
               xi,
               unique(xi)
          )
     )
}

# Ordena las etiquetas según la media de los valores de cada grupo
order_partition <- function(xi, values) {
     group_means <- tapply(
          values,
          xi,
          mean
     )
     
     ordered_labels <- names(
          sort(group_means)
     )
     
     label_map <- setNames(
          seq_along(ordered_labels),
          ordered_labels
     )
     
     as.integer(
          label_map[
               as.character(xi)
          ]
     )
}

# Determina el número de grupos mediante el ancho promedio de silueta
find_optimal_k <- function(
          delta,
          max_k = 8L,
          nstart = 20L
) {
     delta <- as.numeric(delta)
     n <- length(delta)
     
     max_k <- min(
          as.integer(max_k),
          n - 1L,
          length(unique(delta))
     )
     
     if (max_k < 2L) {
          return(1L)
     }
     
     delta_matrix <- matrix(
          delta,
          ncol = 1
     )
     
     delta_distance <- dist(
          delta_matrix
     )
     
     k_values <- seq.int(
          from = 2L,
          to = max_k
     )
     
     silhouette_width <- vapply(
          X = k_values,
          FUN = function(k) {
               fit <- tryCatch(
                    kmeans(
                         x = delta_matrix,
                         centers = k,
                         nstart = nstart
                    ),
                    error = function(e) {
                         NULL
                    }
               )
               
               if (is.null(fit)) {
                    return(-Inf)
               }
               
               silhouette_values <- cluster::silhouette(
                    x = fit$cluster,
                    dist = delta_distance
               )
               
               mean(
                    silhouette_values[, "sil_width"]
               )
          },
          FUN.VALUE = numeric(1)
     )
     
     if (all(!is.finite(silhouette_width))) {
          return(1L)
     }
     
     k_values[
          which.max(silhouette_width)
     ]
}

# Particiones posteriores ------------------------------------------------------

n_draws <- n_samples

xi_draws <- matrix(
     NA_integer_,
     nrow = n_draws,
     ncol = n
)

K_draws <- integer(n_draws)

for (b in seq_len(n_draws)) {
     delta_b <- samples$delta[b, ]
     
     K_draws[b] <- find_optimal_k(
          delta = delta_b,
          max_k = min(8L, n - 1L),
          nstart = 20L
     )
     
     if (K_draws[b] == 1L) {
          xi_b <- rep.int(
               1L,
               n
          )
     } else {
          xi_b <- kmeans(
               x = matrix(
                    delta_b,
                    ncol = 1
               ),
               centers = K_draws[b],
               nstart = 20L
          )$cluster
     }
     
     xi_b <- relabel_partition(
          xi_b
     )
     
     xi_draws[b, ] <- order_partition(
          xi = xi_b,
          values = delta_b
     )
}

# Matriz posterior de coagrupamiento -------------------------------------------

posterior_similarity <- matrix(
     0,
     nrow = n,
     ncol = n
)

for (b in seq_len(n_draws)) {
     incidence_b <- outer(
          xi_draws[b, ],
          xi_draws[b, ],
          FUN = "=="
     )
     
     posterior_similarity <- posterior_similarity +
          incidence_b
}

posterior_similarity <- posterior_similarity / n_draws

diag(posterior_similarity) <- 1

# Estimación puntual mediante la pérdida de Binder -----------------------------

upper_indices <- upper.tri(
     posterior_similarity
)

pairwise_probabilities <- posterior_similarity[
     upper_indices
]

binder_loss <- numeric(n_draws)

for (b in seq_len(n_draws)) {
     incidence_b <- outer(
          xi_draws[b, ],
          xi_draws[b, ],
          FUN = "=="
     )
     
     candidate_pairs <- incidence_b[
          upper_indices
     ]
     
     binder_loss[b] <- sum(
          candidate_pairs *
               (1 - pairwise_probabilities) +
               (1 - candidate_pairs) *
               pairwise_probabilities
     )
}

binder_index <- which.min(
     binder_loss
)

xi_hat <- xi_draws[
     binder_index,
]

clusters <- order_partition(
     xi = xi_hat,
     values = delta_mean
)

K_hat <- length(
     unique(clusters)
)

# Orden de los actores para las visualizaciones --------------------------------

clustering_order <- order(
     clusters,
     delta_mean
)

ordered_actor_labels <- actor_labels[
     clustering_order
]

coclustering_ordered <- posterior_similarity[
     clustering_order,
     clustering_order
]

dimnames(coclustering_ordered) <- list(
     ordered_actor_labels,
     ordered_actor_labels
)

# Matriz posterior de coagrupamiento -------------------------------------------

coclustering_colors <- colorRampPalette(
     c(
          "white",
          "yellow",
          "red"
     )
)(200)

pdf(
     file = "lazega_sociabilidad_coagrupamiento.pdf",
     width = 6,
     height = 5.8,
     pointsize = 12
)

par(
     mfrow = c(1, 1),
     mar = c(1, 1, 1, 1)
)

corrplot::corrplot(
     corr = coclustering_ordered,
     is.corr = FALSE,
     method = "color",
     order = "original",
     col = coclustering_colors,
     col.lim = c(0, 1),
     cl.pos = "r",
     cl.ratio = 0.12,
     cl.cex = 0.8,
     tl.pos = "n",
     addgrid.col = NA,
     diag = TRUE,
     mar = c(1, 1, 1, 1)
)

dev.off()

# Efectos de sociabilidad por grupo --------------------------------------------

cluster_delta_df <- data.frame(
     Actor = actor_labels,
     Mean = delta_mean,
     Lower = delta_interval[1, ],
     Upper = delta_interval[2, ],
     Cluster = clusters
)

cluster_delta_df <- cluster_delta_df[
     order(
          cluster_delta_df$Cluster,
          cluster_delta_df$Mean
     ),
]

cluster_delta_df$Order <- seq_len(
     nrow(cluster_delta_df)
)

cluster_colors <- hcl.colors(
     n = K_hat,
     palette = "Dark 3"
)

point_colors <- cluster_colors[
     cluster_delta_df$Cluster
]

label_offset <- 0.04 * diff(
     range(
          c(
               cluster_delta_df$Lower,
               cluster_delta_df$Upper
          )
     )
)

plot_limits <- range(
     c(
          cluster_delta_df$Lower,
          cluster_delta_df$Upper + label_offset
     )
)

cluster_sizes <- table(
     factor(
          cluster_delta_df$Cluster,
          levels = seq_len(K_hat)
     )
)

if (K_hat > 1L) {
     cluster_boundaries <- cumsum(
          cluster_sizes
     )[-K_hat] + 0.5
} else {
     cluster_boundaries <- numeric(0)
}

pdf(
     file = "lazega_sociabilidad_grupos.pdf",
     width = 7,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(1.4, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

plot(
     x = cluster_delta_df$Order,
     y = cluster_delta_df$Mean,
     type = "n",
     xlim = c(
          0.5,
          nrow(cluster_delta_df) + 0.5
     ),
     ylim = plot_limits,
     xaxt = "n",
     xlab = "",
     ylab = expression(delta[i]),
     bty = "l"
)

abline(
     h = 0,
     lty = 2
)

if (length(cluster_boundaries) > 0L) {
     abline(
          v = cluster_boundaries,
          lty = 3,
          col = "gray75"
     )
}

arrows(
     x0 = cluster_delta_df$Order,
     y0 = cluster_delta_df$Lower,
     x1 = cluster_delta_df$Order,
     y1 = cluster_delta_df$Upper,
     angle = 90,
     code = 3,
     length = 0,
     lwd = 1.2,
     col = point_colors
)

points(
     x = cluster_delta_df$Order,
     y = cluster_delta_df$Mean,
     pch = 16,
     cex = 0.8,
     col = point_colors
)

text(
     x = cluster_delta_df$Order,
     y = cluster_delta_df$Upper + label_offset,
     labels = cluster_delta_df$Actor,
     srt = 90,
     adj = c(0, 0.5),
     cex = 0.8
)

dev.off()

# Evaluación predictiva posterior ----------------------------------------------

pairs <- which(
     upper.tri(Y),
     arr.ind = TRUE
)

observed_statistics <- network_statistics(g)

B_ppc <- n_samples

replicated_statistics <- matrix(
     NA_real_,
     nrow = B_ppc,
     ncol = length(observed_statistics)
)

colnames(replicated_statistics) <- names(observed_statistics)

set.seed(1702)

for (r in seq_len(B_ppc)) {
     probabilities <- pnorm(
          samples$mu[r] +
               samples$delta[r, pairs[, 1]] +
               samples$delta[r, pairs[, 2]]
     )
     
     y_replicated_upper <- rbinom(
          n = nrow(pairs),
          size = 1,
          prob = probabilities
     )
     
     Y_replicated <- matrix(
          0,
          nrow = n,
          ncol = n
     )
     
     Y_replicated[
          cbind(
               pairs[, 1],
               pairs[, 2]
          )
     ] <- y_replicated_upper
     
     Y_replicated[
          cbind(
               pairs[, 2],
               pairs[, 1]
          )
     ] <- y_replicated_upper
     
     g_replicated <- graph_from_adjacency_matrix(
          adjmatrix = Y_replicated,
          mode = "undirected",
          diag = FALSE
     )
     
     replicated_statistics[r, ] <- network_statistics(
          g_replicated
     )
}

statistic_labels <- c(
     Density = "Densidad",
     Transitivity = "Transitividad",
     Assortativity = "Asortatividad",
     AveragePath = "Distancia geodésica",
     AverageDegree = "Grado promedio",
     SDDegree = "Desviación del grado"
)

ppp_values <- vapply(
     X = seq_along(observed_statistics),
     FUN = function(j) {
          mean(
               replicated_statistics[, j] <
                    observed_statistics[j],
               na.rm = TRUE
          )
     },
     FUN.VALUE = numeric(1)
)

ppc_summary <- data.frame(
     Statistic = unname(
          statistic_labels[
               names(observed_statistics)
          ]
     ),
     Observed = as.numeric(observed_statistics),
     Mean = apply(
          replicated_statistics,
          2,
          mean,
          na.rm = TRUE
     ),
     Lower = unname(
          apply(
               replicated_statistics,
               2,
               quantile,
               probs = 0.025,
               na.rm = TRUE
          )
     ),
     Upper = unname(
          apply(
               replicated_statistics,
               2,
               quantile,
               probs = 0.975,
               na.rm = TRUE
          )
     ),
     PPP = ppp_values,
     row.names = NULL
)

ppc_summary

statistic_files <- c(
     Density = "lazega_sociabilidad_ppc_densidad.pdf",
     Transitivity = "lazega_sociabilidad_ppc_transitividad.pdf",
     Assortativity = "lazega_sociabilidad_ppc_asortatividad.pdf",
     AveragePath = "lazega_sociabilidad_ppc_distancia.pdf",
     AverageDegree = "lazega_sociabilidad_ppc_grado_promedio.pdf",
     SDDegree = "lazega_sociabilidad_ppc_desviacion_grado.pdf"
)

for (j in seq_along(observed_statistics)) {
     statistic_name <- names(observed_statistics)[j]
     
     replicated_values <- replicated_statistics[, j]
     replicated_values <- replicated_values[is.finite(replicated_values)]
     
     observed_value <- observed_statistics[j]
     ppp_value <- mean(replicated_values < observed_value)
     
     pdf(
          file = statistic_files[statistic_name],
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
          replicated_values,
          breaks = 20,
          probability = TRUE,
          col = "gray85",
          border = "white",
          xlab = unname(statistic_labels[statistic_name]),
          ylab = "Densidad",
          main = ""
     )
     
     abline(
          v = observed_value,
          lwd = 2,
          lty = 2
     )
     
     legend(
          x = "topright",
          legend = paste0(
               "ppp = ",
               formatC(
                    ppp_value,
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

# Fin --------------------------------------------------------------------------