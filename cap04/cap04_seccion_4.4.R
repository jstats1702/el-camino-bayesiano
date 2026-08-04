# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Librerías
suppressMessages(suppressWarnings(library(igraph)))
suppressMessages(suppressWarnings(library(ggraph)))

# Grafo ------------------------------------------------------------------------

lazega <- graph_from_data_frame(d = elist.lazega, directed = FALSE)

# Clase de objeto
class(lazega)

# ¿Dirigida?
is_directed(lazega)

# Orden
n <- vcount(lazega)
n

# Tamaño
s <- ecount(lazega)
s

# Número de díadas
m <- n * (n - 1) / 2
m

# Visualización del grafo ------------------------------------------------------

pdf(file = "lazega_binomial_viz_grafo.pdf", width = 5, height = 5, pointsize = 17)

par(mfrow = c(1, 1), mar = 0 * c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

set.seed(123)

graph_layout <- layout_with_kk(lazega)

plot(
     lazega,
     layout             = graph_layout,
     vertex.label       = 1:n,
     vertex.label.color = "black",
     vertex.label.cex   = 0.8,
     vertex.size        = 15,
     vertex.color       = 0,
     vertex.frame.color = "black",
     edge.color         = adjustcolor("blue4", alpha.f = 0.5),
     edge.curved        = 0,
     main               = ""
)

dev.off()

# Visualización de la matriz de adyacencia -------------------------------------

pdf(file = "lazega_binomial_viz_matriz.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

# Matriz de adyacencia
Y <- as.matrix(as_adjacency_matrix(graph = lazega, names = FALSE))

corrplot::corrplot(
     corr        = Y,
     is.corr     = FALSE,
     method      = "color",
     type        = "full",
     col         = c("white", "blue4"),
     col.lim     = c(0, 1),
     tl.pos      = "lt",
     tl.col      = "black",
     tl.cex      = 0.55,
     tl.srt      = 90,
     addgrid.col = NA,
     cl.pos      = "n",
     diag        = TRUE,
     mar         = c(0, 0, 0, 0)
)

rect(
     xleft   = 0.5,
     ybottom = 0.5,
     xright  = ncol(Y) + 0.5,
     ytop    = nrow(Y) + 0.5,
     border  = "black",
     lwd     = 1,
     xpd     = NA
)

dev.off()

# Modelo Beta-Binomial ---------------------------------------------------------

# Distribución previa
a <- 1
b <- 1

# Distribución posterior
ap <- a + s
ap

bp <- b + m - s
bp

# Número de muestras de Monte Carlo
B <- 10000

# Simulación de Monte Carlo
set.seed(123)
theta_mc <- rbeta(n = B, shape1 = ap, shape2 = bp)

# Inferencia posterior
out <- c(
     mean(theta_mc),
     sd(theta_mc) / abs(mean(theta_mc)),
     quantile(x = theta_mc, probs = c(0.025, 0.975))
)

names(out) <- c("Estimación", "CV", "Q2.5%", "Q97.5%")
round(out, 4)

# Estadísticos observados ------------------------------------------------------

t0 <- as.matrix(c(
     mean(igraph::degree(lazega)),
     sd(igraph::degree(lazega)),
     igraph::edge_density(lazega),
     igraph::transitivity(lazega),
     igraph::assortativity_degree(lazega, directed = FALSE),
     igraph::mean_distance(lazega, directed = FALSE, unconnected = TRUE)
))

colnames(t0) <- c("Lazega")

rownames(t0) <- c(
     "grado_prom",
     "grado_desv",
     "densidad",
     "transitividad",
     "asortatividad",
     "distancia"
)

round(t0, 3)

# Distribución predictiva posterior -------------------------------------------

t_mc <- NULL

set.seed(123)

for (i in 1:B) {
     # Datos replicados
     Y <- matrix(data = 0, nrow = n, ncol = n)
     Y[lower.tri(Y)] <- rbinom(n = n * (n - 1) / 2, size = 1, prob = theta_mc[i])
     
     g <- igraph::graph_from_adjacency_matrix(
          adjmatrix = Y + t(Y),
          mode      = "undirected"
     )
     
     # Estadísticos de prueba
     t_mc <- rbind(
          t_mc,
          c(
               mean(igraph::degree(g)),
               sd(igraph::degree(g)),
               igraph::edge_density(g),
               igraph::transitivity(g),
               igraph::assortativity_degree(g, directed = FALSE),
               igraph::mean_distance(g, directed = FALSE, unconnected = TRUE)
          )
     )
}

colnames(t_mc) <- rownames(t0)
t_mc[is.nan(t_mc)] <- NA

# Visualizaciones --------------------------------------------------------------

nombres <- rownames(t0)

for (k in 1:6) {
     pdf(
          file      = paste0("lazega_binomial_chequeo_", nombres[k], ".pdf"),
          width     = 5,
          height    = 5,
          pointsize = 20
     )
     
     par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))
     
     # Estadístico observado
     t_obs <- t0[k, 1]
     
     # Valores predictivos simulados
     t_rep <- t_mc[, k]
     t_rep <- t_rep[is.finite(t_rep)]
     
     # Valor p predictivo posterior
     ppp <- round(mean(t_rep <= t_obs), 3)
     cat(nombres[k], "ppp = ", round(ppp, 4), "\n", sep = " ")
     
     hist(
          x      = t_rep,
          freq   = FALSE,
          col    = "gold2",
          border = "gold2",
          xlab   = "t",
          ylab   = expression(p(t ~ "|" ~ bold(y))),
          main   = ""
     )
     
     abline(v = t_obs, col = 1, lwd = 2, lty = 1)
     
     dev.off()
}

# Fin --------------------------------------------------------------------------