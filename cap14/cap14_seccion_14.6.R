# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Librerías
suppressMessages(
     suppressWarnings({
          library(igraph)
          library(coda)
          library(corrplot)
     })
)

# Funciones auxiliares ---------------------------------------------------------

# Genera una muestra de una distribución Dirichlet.
sample_dirichlet <- function(alpha) {
     values <- rgamma(
          n = length(alpha),
          shape = alpha,
          rate = 1
     )
     
     values / sum(values)
}

# Calcula el número esperado de bloques ocupados.
expected_occupied_blocks <- function(alpha, K, n) {
     probability_empty <- exp(
          lbeta(
               alpha / K,
               alpha - alpha / K + n
          ) -
               lbeta(
                    alpha / K,
                    alpha - alpha / K
               )
     )
     
     K * (1 - probability_empty)
}

# Calibra alpha a partir del número esperado de bloques ocupados.
calibrate_alpha <- function(target, K, n) {
     objective <- function(log_alpha) {
          expected_occupied_blocks(
               alpha = exp(log_alpha),
               K = K,
               n = n
          ) - target
     }
     
     exp(
          uniroot(
               f = objective,
               interval = c(
                    log(1e-6),
                    log(1e6)
               )
          )$root
     )
}

# Construye la lista de díadas de una red no dirigida.
get_pairs <- function(n) {
     which(
          upper.tri(
               matrix(
                    data = 0,
                    nrow = n,
                    ncol = n
               )
          ),
          arr.ind = TRUE
     )
}

# Inicializa la partición mediante agrupamiento jerárquico.
initialize_partition <- function(Y, K) {
     dissimilarity <- dist(
          x = Y,
          method = "euclidean"
     )
     
     clustering <- hclust(
          d = dissimilarity,
          method = "ward.D2"
     )
     
     as.integer(
          cutree(
               tree = clustering,
               k = min(K, nrow(Y))
          )
     )
}

# Calcula el número de díadas y relaciones entre bloques.
block_counts <- function(Y, xi, K, pairs) {
     m <- matrix(
          data = 0,
          nrow = K,
          ncol = K
     )
     
     s <- matrix(
          data = 0,
          nrow = K,
          ncol = K
     )
     
     for (r in seq_len(nrow(pairs))) {
          i <- pairs[r, 1]
          j <- pairs[r, 2]
          
          k <- min(
               xi[i],
               xi[j]
          )
          
          ell <- max(
               xi[i],
               xi[j]
          )
          
          m[k, ell] <- m[k, ell] + 1
          s[k, ell] <- s[k, ell] + Y[i, j]
     }
     
     list(
          m = m,
          s = s
     )
}

# Calcula las probabilidades diádicas para una iteración.
dyadic_probabilities <- function(xi, theta, pairs) {
     block_1 <- pmin(
          xi[pairs[, 1]],
          xi[pairs[, 2]]
     )
     
     block_2 <- pmax(
          xi[pairs[, 1]],
          xi[pairs[, 2]]
     )
     
     theta[
          cbind(
               block_1,
               block_2
          )
     ]
}

# Calcula la log-verosimilitud condicional.
log_likelihood_sbm <- function(Y, xi, theta, pairs) {
     probabilities <- dyadic_probabilities(
          xi = xi,
          theta = theta,
          pairs = pairs
     )
     
     probabilities <- pmin(
          pmax(
               probabilities,
               1e-12
          ),
          1 - 1e-12
     )
     
     sum(
          dbinom(
               x = Y[
                    cbind(
                         pairs[, 1],
                         pairs[, 2]
                    )
               ],
               size = 1,
               prob = probabilities,
               log = TRUE
          )
     )
}

# Ajusta el modelo de bloques estocásticos mediante Gibbs.
gibbs_sampler_sbm <- function(
          Y,
          K,
          B,
          burn,
          thin,
          a_theta,
          b_theta,
          alpha,
          seed = 42,
          verbose = TRUE,
          progress_every = 10000L
) {
     set.seed(seed)
     
     n <- nrow(Y)
     pairs <- get_pairs(n)
     total_iterations <- burn + B * thin
     
     xi <- initialize_partition(
          Y = Y,
          K = K
     )
     
     omega <- sample_dirichlet(
          alpha = alpha / K +
               tabulate(
                    xi,
                    nbins = K
               )
     )
     
     theta <- matrix(
          data = 0,
          nrow = K,
          ncol = K
     )
     
     counts <- block_counts(
          Y = Y,
          xi = xi,
          K = K,
          pairs = pairs
     )
     
     for (k in seq_len(K)) {
          for (ell in k:K) {
               theta[k, ell] <- rbeta(
                    n = 1,
                    shape1 = a_theta + counts$s[k, ell],
                    shape2 = b_theta +
                         counts$m[k, ell] -
                         counts$s[k, ell]
               )
               
               theta[ell, k] <- theta[k, ell]
          }
     }
     
     samples <- list(
          xi = matrix(
               data = NA_integer_,
               nrow = B,
               ncol = n
          ),
          theta = array(
               data = NA_real_,
               dim = c(
                    B,
                    K,
                    K
               )
          ),
          omega = matrix(
               data = NA_real_,
               nrow = B,
               ncol = K
          ),
          K_star = integer(B),
          mean_probability = numeric(B),
          log_likelihood = numeric(B)
     )
     
     stored <- 0L
     
     for (iteration in seq_len(total_iterations)) {
          counts <- block_counts(
               Y = Y,
               xi = xi,
               K = K,
               pairs = pairs
          )
          
          for (k in seq_len(K)) {
               for (ell in k:K) {
                    theta[k, ell] <- rbeta(
                         n = 1,
                         shape1 = a_theta + counts$s[k, ell],
                         shape2 = b_theta +
                              counts$m[k, ell] -
                              counts$s[k, ell]
                    )
                    
                    theta[ell, k] <- theta[k, ell]
               }
          }
          
          for (i in seq_len(n)) {
               other <- setdiff(
                    seq_len(n),
                    i
               )
               
               log_weight <- numeric(K)
               
               for (k in seq_len(K)) {
                    block_1 <- pmin(
                         k,
                         xi[other]
                    )
                    
                    block_2 <- pmax(
                         k,
                         xi[other]
                    )
                    
                    probabilities <- theta[
                         cbind(
                              block_1,
                              block_2
                         )
                    ]
                    
                    probabilities <- pmin(
                         pmax(
                              probabilities,
                              1e-12
                         ),
                         1 - 1e-12
                    )
                    
                    log_weight[k] <- log(omega[k]) +
                         sum(
                              dbinom(
                                   x = Y[i, other],
                                   size = 1,
                                   prob = probabilities,
                                   log = TRUE
                              )
                         )
               }
               
               log_weight <- log_weight - max(log_weight)
               
               assignment_probability <- exp(log_weight)
               assignment_probability <- assignment_probability /
                    sum(assignment_probability)
               
               xi[i] <- sample.int(
                    n = K,
                    size = 1,
                    prob = assignment_probability
               )
          }
          
          block_sizes <- tabulate(
               xi,
               nbins = K
          )
          
          omega <- sample_dirichlet(
               alpha = alpha / K + block_sizes
          )
          
          if (
               iteration > burn &&
               (iteration - burn) %% thin == 0L
          ) {
               stored <- stored + 1L
               
               probabilities <- dyadic_probabilities(
                    xi = xi,
                    theta = theta,
                    pairs = pairs
               )
               
               samples$xi[stored, ] <- xi
               samples$theta[stored, , ] <- theta
               samples$omega[stored, ] <- omega
               samples$K_star[stored] <- sum(block_sizes > 0)
               samples$mean_probability[stored] <- mean(probabilities)
               
               samples$log_likelihood[stored] <- log_likelihood_sbm(
                    Y = Y,
                    xi = xi,
                    theta = theta,
                    pairs = pairs
               )
          }
          
          if (
               verbose &&
               (
                    iteration == 1L ||
                    iteration %% progress_every == 0L ||
                    iteration == total_iterations
               )
          ) {
               cat(
                    sprintf(
                         "Iteración %d de %d | bloques ocupados = %d\n",
                         iteration,
                         total_iterations,
                         sum(block_sizes > 0)
                    )
               )
          }
     }
     
     list(
          Y = Y,
          settings = list(
               n = n,
               K = K,
               B = B,
               burn = burn,
               thin = thin,
               total_iterations = total_iterations,
               seed = seed
          ),
          hyperparameters = list(
               a_theta = a_theta,
               b_theta = b_theta,
               alpha = alpha
          ),
          samples = samples
     )
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
          de = sd(x),
          ess = ess,
          mcse = sd(x) / sqrt(ess)
     )
}

# Reetiqueta una partición según el menor índice de cada bloque.
relabel_partition <- function(xi) {
     group_minimum <- tapply(
          X = seq_along(xi),
          INDEX = xi,
          FUN = min
     )
     
     ordered_labels <- names(
          sort(group_minimum)
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

# Calcula la matriz posterior de coagrupamiento.
posterior_similarity_matrix <- function(xi_draws) {
     B <- nrow(xi_draws)
     n <- ncol(xi_draws)
     
     similarity <- matrix(
          data = 0,
          nrow = n,
          ncol = n
     )
     
     for (b in seq_len(B)) {
          similarity <- similarity +
               outer(
                    xi_draws[b, ],
                    xi_draws[b, ],
                    FUN = "=="
               )
     }
     
     similarity / B
}

# Obtiene una partición puntual mediante la pérdida de Binder.
binder_partition <- function(xi_draws, similarity) {
     B <- nrow(xi_draws)
     upper_indices <- upper.tri(similarity)
     
     pairwise_probability <- similarity[
          upper_indices
     ]
     
     loss <- numeric(B)
     
     for (b in seq_len(B)) {
          incidence <- outer(
               xi_draws[b, ],
               xi_draws[b, ],
               FUN = "=="
          )
          
          candidate <- incidence[
               upper_indices
          ]
          
          loss[b] <- sum(
               candidate *
                    (1 - pairwise_probability) +
                    (1 - candidate) *
                    pairwise_probability
          )
     }
     
     index <- which.min(loss)
     
     list(
          partition = relabel_partition(
               xi_draws[index, ]
          ),
          index = index,
          loss = loss
     )
}

# Calcula el índice de Rand ajustado.
adjusted_rand_index <- function(partition_1, partition_2) {
     contingency <- table(
          partition_1,
          partition_2
     )
     
     choose_2 <- function(x) {
          x * (x - 1) / 2
     }
     
     index <- sum(
          choose_2(contingency)
     )
     
     row_sum <- rowSums(contingency)
     column_sum <- colSums(contingency)
     
     expected <- sum(choose_2(row_sum)) *
          sum(choose_2(column_sum)) /
          choose_2(sum(contingency))
     
     maximum <- 0.5 * (
          sum(choose_2(row_sum)) +
               sum(choose_2(column_sum))
     )
     
     (index - expected) / (maximum - expected)
}

# Calcula las medias posteriores de las probabilidades diádicas.
posterior_probability_matrix <- function(fit) {
     samples <- fit$samples
     B <- nrow(samples$xi)
     n <- fit$settings$n
     pairs <- get_pairs(n)
     
     probability_matrix <- matrix(
          data = 0,
          nrow = n,
          ncol = n
     )
     
     for (b in seq_len(B)) {
          probabilities <- dyadic_probabilities(
               xi = samples$xi[b, ],
               theta = samples$theta[b, , ],
               pairs = pairs
          )
          
          probability_matrix[
               cbind(
                    pairs[, 1],
                    pairs[, 2]
               )
          ] <- probability_matrix[
               cbind(
                    pairs[, 1],
                    pairs[, 2]
               )
          ] + probabilities / B
          
          probability_matrix[
               cbind(
                    pairs[, 2],
                    pairs[, 1]
               )
          ] <- probability_matrix[
               cbind(
                    pairs[, 2],
                    pairs[, 1]
               )
          ] + probabilities / B
     }
     
     probability_matrix
}

# Resume las probabilidades según una partición puntual.
block_probability_summary <- function(probability_matrix, partition) {
     K_hat <- length(
          unique(partition)
     )
     
     result <- matrix(
          data = NA_real_,
          nrow = K_hat,
          ncol = K_hat
     )
     
     for (k in seq_len(K_hat)) {
          for (ell in k:K_hat) {
               index_k <- which(partition == k)
               index_ell <- which(partition == ell)
               
               if (k == ell) {
                    block_values <- probability_matrix[
                         index_k,
                         index_k,
                         drop = FALSE
                    ]
                    
                    block_values <- block_values[
                         upper.tri(block_values)
                    ]
               } else {
                    block_values <- probability_matrix[
                         index_k,
                         index_ell,
                         drop = FALSE
                    ]
               }
               
               result[k, ell] <- mean(block_values)
               result[ell, k] <- result[k, ell]
          }
     }
     
     dimnames(result) <- list(
          paste0(
               "Bloque ",
               seq_len(K_hat)
          ),
          paste0(
               "Bloque ",
               seq_len(K_hat)
          )
     )
     
     result
}

# Calcula estadísticos de una red.
network_statistics <- function(A) {
     graph <- igraph::graph_from_adjacency_matrix(
          adjmatrix = A,
          mode = "undirected",
          diag = FALSE
     )
     
     values <- c(
          densidad = igraph::edge_density(
               graph = graph,
               loops = FALSE
          ),
          transitividad = igraph::transitivity(
               graph = graph,
               type = "global"
          ),
          asortatividad = igraph::assortativity_degree(
               graph = graph,
               directed = FALSE
          ),
          distancia = igraph::mean_distance(
               graph = graph,
               directed = FALSE,
               unconnected = TRUE
          ),
          grado_promedio = mean(
               igraph::degree(graph)
          ),
          desviacion_grado = sd(
               igraph::degree(graph)
          )
     )
     
     values[!is.finite(values)] <- NA_real_
     
     values
}

# Realiza comprobaciones predictivas posteriores.
posterior_predictive_check_sbm <- function(
          fit,
          A_observed,
          B_ppc = 5000L,
          seed = NULL
) {
     if (!is.null(seed)) {
          set.seed(seed)
     }
     
     samples <- fit$samples
     B_post <- nrow(samples$xi)
     n <- fit$settings$n
     pairs <- get_pairs(n)
     
     B_ppc <- min(
          as.integer(B_ppc),
          B_post
     )
     
     index <- sample(
          x = seq_len(B_post),
          size = B_ppc,
          replace = FALSE
     )
     
     observed_statistics <- network_statistics(
          A = A_observed
     )
     
     replicated_statistics <- matrix(
          data = NA_real_,
          nrow = B_ppc,
          ncol = length(observed_statistics)
     )
     
     colnames(replicated_statistics) <- names(
          observed_statistics
     )
     
     for (r in seq_len(B_ppc)) {
          b <- index[r]
          
          probabilities <- dyadic_probabilities(
               xi = samples$xi[b, ],
               theta = samples$theta[b, , ],
               pairs = pairs
          )
          
          y_replicated <- rbinom(
               n = nrow(pairs),
               size = 1,
               prob = probabilities
          )
          
          A_replicated <- matrix(
               data = 0,
               nrow = n,
               ncol = n
          )
          
          A_replicated[
               cbind(
                    pairs[, 1],
                    pairs[, 2]
               )
          ] <- y_replicated
          
          A_replicated[
               cbind(
                    pairs[, 2],
                    pairs[, 1]
               )
          ] <- y_replicated
          
          replicated_statistics[r, ] <- network_statistics(
               A = A_replicated
          )
     }
     
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
     
     names(ppp) <- names(
          observed_statistics
     )
     
     list(
          observed = observed_statistics,
          replicated = replicated_statistics,
          ppp = ppp
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
     replicated <- replicated[
          is.finite(replicated)
     ]
     
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

g <- igraph::make_graph(
     edges = "Zachary"
)

g <- igraph::simplify(
     graph = igraph::as_undirected(g)
)

n <- igraph::vcount(g)

igraph::V(g)$name <- as.character(
     seq_len(n)
)

A <- as.matrix(
     igraph::as_adjacency_matrix(
          graph = g,
          sparse = FALSE
     )
)

storage.mode(A) <- "numeric"

node_labels <- igraph::V(g)$name

faction <- rep(
     "Administrador",
     n
)

faction[
     c(
          1, 2, 3, 4, 5, 6, 7, 8, 9,
          11, 12, 13, 14, 17, 18, 20, 22
     )
] <- "Instructor"

s <- igraph::ecount(g)
m <- choose(n, 2)
density <- igraph::edge_density(
     graph = g,
     loops = FALSE
)

c(
     actores = n,
     aristas = s,
     diadas = m,
     densidad = density
)

# Visualización de los datos ---------------------------------------------------

set.seed(1702)

graph_layout <- igraph::layout_with_fr(
     graph = g
)

# Grafo
pdf(
     file = "karate_sbm_grafo.pdf",
     width = 5,
     height = 5,
     pointsize = 17
)

par(
     mar = c(1, 1, 1, 1)
)

plot(
     g,
     layout = graph_layout,
     vertex.label = node_labels,
     vertex.label.color = "black",
     vertex.label.cex = 0.8,
     vertex.size = 15,
     vertex.color = "white",
     vertex.frame.color = "black",
     vertex.shape = ifelse(
          faction == "Instructor",
          "circle",
          "square"
     ),
     edge.color = adjustcolor(
          "blue4",
          alpha.f = 0.8
     ),
     edge.width = 0.7,
     edge.curved = 0,
     main = ""
)

legend(
     x = "topright",
     legend = c(
          "Instructor",
          "Administrador"
     ),
     pch = c(
          21,
          22
     ),
     pt.bg = "white",
     col = "black",
     bty = "n",
     cex = 0.9
)

dev.off()

# Matriz de adyacencia
pdf(
     file = "karate_sbm_matriz.pdf",
     width = 5,
     height = 5,
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
     col = c(
          "white",
          "blue4"
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
     xright = ncol(A) + 0.5,
     ytop = nrow(A) + 0.5,
     border = "black",
     lwd = 1,
     xpd = NA
)

dev.off()

# Ajuste del modelo ------------------------------------------------------------

K <- 8L
K_target <- 4

a_theta <- 1
b_theta <- 1

alpha <- calibrate_alpha(
     target = K_target,
     K = K,
     n = n
)

c(
     K = K,
     bloques_esperados = expected_occupied_blocks(
          alpha = alpha,
          K = K,
          n = n
     ),
     alpha = alpha
)

fit_sbm <- gibbs_sampler_sbm(
     Y = A,
     K = K,
     B    = 10000,
     burn = 10000,
     thin = 10,
     a_theta = a_theta,
     b_theta = b_theta,
     alpha = alpha,
     seed = 1702,
     verbose = TRUE,
     progress_every = 10000L
)

saveRDS(
     object = fit_sbm,
     file = "karate_sbm_modelo.rds"
)

fit_sbm <- readRDS(
     file = "karate_sbm_modelo.rds"
)

samples <- fit_sbm$samples

# Diagnósticos -----------------------------------------------------------------

diagnostic_table <- rbind(
     K_star = diagnostic_summary(
          samples$K_star
     ),
     theta_bar = diagnostic_summary(
          samples$mean_probability
     ),
     logL = diagnostic_summary(
          samples$log_likelihood
     )
)

round(
     diagnostic_table,
     digits = 3
)

# Trazas
par(
     mfrow = c(3, 1),
     mar = c(3, 3, 1.2, 1.2),
     mgp = c(1.75, 0.75, 0)
)

plot(
     samples$K_star,
     type = "p",
     pch  = 16,
     cex  = 0.2,
     xlab = "Iteración",
     ylab = expression(K^"*"),
     main = ""
)

plot(
     samples$mean_probability,
     type = "p",
     pch  = 16,
     cex  = 0.2,
     xlab = "Iteración",
     ylab = expression(bar(theta)),
     main = ""
)

plot(
     samples$log_likelihood,
     type = "p",
     pch  = 16,
     cex  = 0.2,
     xlab = "Iteración",
     ylab = "Log-verosimilitud",
     main = ""
)

dev.off()

# Inferencia sobre el número de bloques ----------------------------------------

K_star_probability <- prop.table(
     table(
          factor(
               samples$K_star,
               levels = seq_len(K)
          )
     )
)

K_star_table <- data.frame(
     K_star = seq_len(K),
     probabilidad = as.numeric(
          K_star_probability
     )
)

K_star_table <- K_star_table[
     K_star_table$probabilidad > 0,
]

K_star_table

K_star_mode <- K_star_table$K_star[
     which.max(
          K_star_table$probabilidad
     )
]

K_star_mean <- mean(
     samples$K_star
)

c(
     media = K_star_mean,
     moda = K_star_mode
)

pdf(
     file = "karate_sbm_numero_bloques.pdf",
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
     x = K_star_table$K_star,
     y = K_star_table$probabilidad,
     type = "h",
     lwd = 5,
     lend = 1,
     xlab = "Número de bloques ocupados",
     ylab = "Probabilidad posterior",
     xaxt = "n",
     ylim = c(
          0,
          max(K_star_table$probabilidad) * 1.05
     ),
     col = "gray40"
)

axis(
     side = 1,
     at = K_star_table$K_star
)

box()

dev.off()

# Inferencia sobre la partición ------------------------------------------------

posterior_similarity <- posterior_similarity_matrix(
     xi_draws = samples$xi
)

binder_result <- binder_partition(
     xi_draws = samples$xi,
     similarity = posterior_similarity
)

partition_hat <- binder_result$partition

K_hat <- length(
     unique(partition_hat)
)

partition_table <- data.frame(
     actor = node_labels,
     bloque = partition_hat,
     faccion = faction
)

partition_table <- partition_table[
     order(
          partition_table$bloque,
          as.integer(partition_table$actor)
     ),
]

partition_table

block_sizes <- table(
     partition_hat
)

block_sizes

faction_numeric <- ifelse(
     faction == "Instructor",
     1L,
     2L
)

ari <- adjusted_rand_index(
     partition_1 = partition_hat,
     partition_2 = faction_numeric
)

ari

table(
     bloque = partition_hat,
     faccion = faction
)

partition_order <- order(
     partition_hat,
     -igraph::degree(g)
)

ordered_labels <- node_labels[
     partition_order
]

similarity_ordered <- posterior_similarity[
     partition_order,
     partition_order
]

dimnames(similarity_ordered) <- list(
     ordered_labels,
     ordered_labels
)

cluster_colors <- grDevices::hcl.colors(
     n = max(
          K_hat,
          3L
     ),
     palette = "Dark 3"
)[
     seq_len(K_hat)
]

block_sizes_ordered <- as.integer(
     block_sizes
)

block_boundaries <- cumsum(
     block_sizes_ordered
)

block_boundaries <- block_boundaries[
     -length(block_boundaries)
]

n_nodes <- nrow(
     similarity_ordered
)

# Grafo con comunidades
pdf(
     file = "karate_sbm_particion.pdf",
     width = 5,
     height = 5,
     pointsize = 17
)

par(
     mar = c(1, 1, 1, 1)
)

plot(
     g,
     layout = graph_layout,
     vertex.label = node_labels,
     vertex.label.color = "black",
     vertex.label.cex = 0.8,
     vertex.size = 12,
     vertex.color = cluster_colors[
          partition_hat
     ],
     vertex.frame.color = "black",
     vertex.shape = ifelse(
          faction == "Instructor",
          "circle",
          "square"
     ),
     edge.color = "gray70",
     edge.width = 0.7,
     edge.curved = 0,
     main = ""
)

legend(
     x = "topright",
     legend = paste(
          "Bloque",
          seq_len(K_hat)
     ),
     pch = 21,
     pt.bg = cluster_colors,
     col = "black",
     bty = "n",
     cex = 0.8
)

legend(
     x = "top",
     legend = c(
          "Instructor",
          "Administrador"
     ),
     pch = c(
          21,
          22
     ),
     pt.bg = "white",
     col = "black",
     bty = "n",
     cex = 0.8
)

dev.off()

# Matriz de coagrupamiento
coagrupamiento_colors <- colorRampPalette(
     c(
          "white",
          "yellow",
          "red"
     )
)(200)

pdf(
     file = "karate_sbm_coagrupamiento.pdf",
     width = 5,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = 0 * c(1, 1, 1, 1)
)

corrplot::corrplot(
     corr = similarity_ordered,
     is.corr = FALSE,
     method = "color",
     type = "full",
     col = coagrupamiento_colors,
     col.lim = c(0, 1),
     cl.pos = "n",
     tl.col = "black",
     tl.cex = 0.7,
     tl.srt = 90,
     addgrid.col = NA,
     diag = TRUE,
     mar = 0 * c(1, 1, 1, 1),
     title = ""
)

for (boundary in block_boundaries) {
     segments(
          x0 = boundary + 0.5,
          y0 = 0.5,
          x1 = boundary + 0.5,
          y1 = n_nodes + 0.5,
          col = "black",
          lwd = 1
     )
     
     segments(
          x0 = 0.5,
          y0 = n_nodes - boundary + 0.5,
          x1 = n_nodes + 0.5,
          y1 = n_nodes - boundary + 0.5,
          col = "black",
          lwd = 1
     )
}

rect(
     xleft = 0.5,
     ybottom = 0.5,
     xright = ncol(similarity_ordered) + 0.5,
     ytop = nrow(similarity_ordered) + 0.5,
     border = "black",
     lwd = 1,
     xpd = NA
)

dev.off()

# Probabilidades de interacción ------------------------------------------------

probability_matrix <- posterior_probability_matrix(
     fit = fit_sbm
)

rownames(probability_matrix) <- node_labels
colnames(probability_matrix) <- node_labels

block_probability_matrix <- block_probability_summary(
     probability_matrix = probability_matrix,
     partition = partition_hat
)

round(
     block_probability_matrix,
     digits = 3
)

A_ordered <- A[
     partition_order,
     partition_order
]

probability_ordered <- probability_matrix[
     partition_order,
     partition_order
]

dimnames(A_ordered) <- list(
     ordered_labels,
     ordered_labels
)

dimnames(probability_ordered) <- list(
     ordered_labels,
     ordered_labels
)

# Matriz de adyacencia
pdf(
     file = "karate_sbm_adyacencia.pdf",
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

for (boundary in block_boundaries) {
     segments(
          x0 = boundary + 0.5,
          y0 = 0.5,
          x1 = boundary + 0.5,
          y1 = n_nodes + 0.5,
          col = "black",
          lwd = 1
     )
     
     segments(
          x0 = 0.5,
          y0 = n_nodes - boundary + 0.5,
          x1 = n_nodes + 0.5,
          y1 = n_nodes - boundary + 0.5,
          col = "black",
          lwd = 1
     )
}

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

# Matriz de probabilidades de interacción
probability_colors <- colorRampPalette(
     c(
          "white",
          "yellow",
          "red"
     )
)(200)

pdf(
     file = "karate_sbm_probabilidades.pdf",
     width = 5,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = 0 * c(1, 1, 1, 1)
)

corrplot::corrplot(
     corr = probability_ordered,
     is.corr = FALSE,
     method = "color",
     type = "full",
     col = probability_colors,
     col.lim = c(0, 1),
     cl.pos = "n",
     tl.col = "black",
     tl.cex = 0.7,
     tl.srt = 90,
     addgrid.col = NA,
     diag = TRUE,
     mar = 0 * c(1, 1, 1, 1),
     title = ""
)

for (boundary in block_boundaries) {
     segments(
          x0 = boundary + 0.5,
          y0 = 0.5,
          x1 = boundary + 0.5,
          y1 = n_nodes + 0.5,
          col = "black",
          lwd = 1
     )
     
     segments(
          x0 = 0.5,
          y0 = n_nodes - boundary + 0.5,
          x1 = n_nodes + 0.5,
          y1 = n_nodes - boundary + 0.5,
          col = "black",
          lwd = 1
     )
}

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

# Evaluación predictiva posterior ----------------------------------------------

ppc <- posterior_predictive_check_sbm(
     fit = fit_sbm,
     A_observed = A,
     B_ppc = 5000,
     seed = 2468
)

ppc_table <- data.frame(
     estadistico = names(
          ppc$observed
     ),
     observado = as.numeric(
          ppc$observed
     ),
     ppp = as.numeric(
          ppc$ppp
     )
)

ppc_table

plot_ppc(
     replicated = ppc$replicated[, "grado_promedio"],
     observed = ppc$observed["grado_promedio"],
     ppp = ppc$ppp["grado_promedio"],
     file = "karate_sbm_ppc_grado_promedio.pdf",
     xlab = "Grado promedio"
)

plot_ppc(
     replicated = ppc$replicated[, "desviacion_grado"],
     observed = ppc$observed["desviacion_grado"],
     ppp = ppc$ppp["desviacion_grado"],
     file = "karate_sbm_ppc_desviacion_grado.pdf",
     xlab = "DE del grado"
)

plot_ppc(
     replicated = ppc$replicated[, "densidad"],
     observed = ppc$observed["densidad"],
     ppp = ppc$ppp["densidad"],
     file = "karate_sbm_ppc_densidad.pdf",
     xlab = "Densidad"
)

plot_ppc(
     replicated = ppc$replicated[, "transitividad"],
     observed = ppc$observed["transitividad"],
     ppp = ppc$ppp["transitividad"],
     file = "karate_sbm_ppc_transitividad.pdf",
     xlab = "Transitividad"
)

plot_ppc(
     replicated = ppc$replicated[, "asortatividad"],
     observed = ppc$observed["asortatividad"],
     ppp = ppc$ppp["asortatividad"],
     file = "karate_sbm_ppc_asortatividad.pdf",
     xlab = "Asortatividad por grado"
)

plot_ppc(
     replicated = ppc$replicated[, "distancia"],
     observed = ppc$observed["distancia"],
     ppp = ppc$ppp["distancia"],
     file = "karate_sbm_ppc_distancia.pdf",
     xlab = "Distancia promedio"
)

# Fin --------------------------------------------------------------------------