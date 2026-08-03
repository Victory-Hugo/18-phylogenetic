library(saasi)
library(diversitree)
library(ape)
library(phytools)

tree <- ape::read.tree(
  system.file(
    "extdata",
    "nexstrain_ebola_ebov-2013_smaller.nwk",
    package = "saasi"
  )
)

metadata <- readr::read_tsv(
  system.file(
    "extdata",
    "nextstrain_ebola_ebov-2013_metadata.tsv",
    package = "saasi"
  )
)

tip_data <- data.frame(
  tip_label = metadata$strain,
  state = metadata$country
)

ebola_tree <- prepare_tree_for_saasi(tree, tip_data)

plot_saasi(
  ebola_tree,
  saasi_result = NULL,
  tip_cex = 1,
  res = 900
)

Q <- estimate_transition_rates(
  ebola_tree,
  method = "fitMk",
  matrix_structure = "SYM"
)

rates <- estimate_bds_parameters(
  ebola_tree,
  mu = 5,
  r0_max = 3,
  r0_min = 1.5,
  psi_max = 15,
  infectious_period_min = 20 / 365,
  infectious_period_max = 40 / 365,
  n_starts = 100
)

pars <- data.frame(
  state = colnames(Q),
  lambda = rates$lambda,
  mu = rates$mu,
  psi = rates$psi,
  row.names = NULL
)

pars[pars$state == "Liberia", "psi"] <-
  pars[pars$state == "Liberia", "psi"] / 2

saasi_ebola <- saasi(ebola_tree, Q, pars)

p1 <- plot_saasi(
  ebola_tree,
  saasi_ebola,
  tip_cex = 1,
  node_cex = 0.5,
  res = 900
)

pars <- data.frame(
  state = colnames(Q),
  lambda = rates$lambda,
  mu = rates$mu,
  psi = rates$psi,
  row.names = NULL
)

saasi_ebola2 <- saasi(ebola_tree, Q, pars)


p2 <- plot_saasi(
  ebola_tree,
  saasi_ebola2,
  tip_cex = 1,
  node_cex = 0.5,
  res = 900
)
