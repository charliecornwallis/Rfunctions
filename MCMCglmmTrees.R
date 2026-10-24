# prior=list(R=list(R1=list(V=diag(6), nu=5.002)),
#            G=list(G1=list(V=diag(6), nu=5.002),
#                   G2=list(V=diag(6), nu=5.002),
#                   G3=list(V=diag(6), nu=5.002)))
# treetoburn=0
# data=dat
# samples=10
# trees=trees
# tiplabel=dat$animal
# fixed=as.formula(cbind(Ztemp,Zsqrt_precip,Zsqrt_temp_sd,Zsqrt_precip_sd,Zsqrt_temp_sd_byrs,Zsqrt_precip_sd_byrs)~
#                    at.level(fam_non,'Pair'):trait+
#                    at.level(fam_non,'Family'):trait+
#                    at.level(fam_non,'Nonfamily'):trait-1)
# random=as.formula(~ us(at.level(fam_non,'Pair'):trait):animal+
#                     us(at.level(fam_non,'Family'):trait):animal+
#                     us(at.level(fam_non,'Nonfamily'):trait):animal)
# residual=as.formula(~us(trait):units)
# Gnumber=3
# Rnumber=1
# Geffects=c((6*6),(6*6),(6*6))
# Reffects=(6*6)
# family=c("gaussian","gaussian","gaussian","gaussian","gaussian","gaussian")
# nitts_tree=1
# thins_tree=1
# burns_tree=0

#***************************************
#Function for sampling of across trees using MCMCglmm models#
#***************************************

MCMCglmmTrees <- function(
  prior, treetoburn, data, mev = 0, trees, tiplabel, fixed, random, residual,
  Gnumber, Rnumber, Geffects, Reffects, family = "gaussian", samples = 1000,
  nitts_tree, thins_tree, burns_tree, pr = FALSE, trunc = FALSE,
  phylo_effect = NULL, inverse_tol = .Machine$double.eps^0.5
) {
# Explanation
# prior #model prior
# treetoburn #number of trees to discard before saving iterations
# data #data to be used
# mev=0 #is there a mev effect in the model mev=0 (no), mev=1 (yes)
# trees #trees to sample across
# tiplabel #name of tiplabel in dataset
# fixed #fixed effect formula
# random #G side random effects formula
# residual #R side random effects formula
# Gnumber #number of G side random effects (Number of G in prior list)
# Rnumber #number of R side random effects (Number of R in prior list)
# Geffects #number of parameters per G term (e.g. 2x2 covar per G = 4) given as a vector [e.g c(1,4,4) ] with length = Gnumber [e.g. 3: G1 = single variance, G2=2x2, G3=2x2]
# Reffects #
# family="gaussian" # vector of families e.g. c("gaussian","gaussian") for 2 trait multi-response model
# samples=1000
# nitts_tree
# thins_tree
# burns_tree

  stopifnot(
    length(trees) > 0L, samples > 0, samples == as.integer(samples),
    treetoburn >= 0, treetoburn == as.integer(treetoburn),
    length(tiplabel) == nrow(data), !anyNA(tiplabel),
    is.finite(inverse_tol), inverse_tol > 0
  )

  # The inverse name must match the factor actually used in the random formula.
  # Do not silently attach a tiplabel inverse to an unrelated animal effect.
  random_variables <- all.vars(random)
  if (is.null(phylo_effect)) {
    candidates <- intersect(c("animal", "tiplabel"), random_variables)
    if (length(candidates) != 1L) {
      stop("Specify phylo_effect: no unique animal/tiplabel random effect.")
    }
    phylo_effect <- candidates
  }
  if (length(phylo_effect) != 1L || !phylo_effect %in% random_variables) {
    stop("phylo_effect must name a factor present in the random formula.")
  }
  data$tiplabel <- tiplabel
  if (!phylo_effect %in% names(data) || !identical(
    as.character(data[[phylo_effect]]), as.character(tiplabel)
  )) {
    stop("The phylogenetic random-effect values must match tiplabel by row.")
  }

  # Preserve the legacy order when enough distinct input trees are supplied.
  # For 100 trees and 1,500 updates, shuffle complete equal-weight cycles;
  # 500 warm-up updates and 1,000 retained updates remain unchanged.
  updates <- as.integer(samples + treetoburn)
  if (length(trees) >= updates) {
    tree_indices <- seq_len(updates)
  } else {
    tree_indices <- integer(0)
    for (cycle in seq_len(ceiling(updates / length(trees)))) {
      tree_indices <- c(tree_indices, sample.int(length(trees)))
    }
    tree_indices <- tree_indices[seq_len(updates)]
  }

  tree1 <- MCMCglmm::inverseA(
    trees[[tree_indices[1]]], tol = inverse_tol
  )$Ainv
  if (!all(as.character(tiplabel) %in% rownames(tree1))) {
    stop("The input tree does not contain every analysis tiplabel.")
  }
  inverse_list <- setNames(list(tree1), phylo_effect)
  model <- MCMCglmm::MCMCglmm(
    fixed = fixed, random = random, rcov = residual,
    ginverse = inverse_list, family = family, data = data, prior = prior,
    nitt = samples, burnin = 0, thin = 1, pr = pr, verbose = FALSE,
    pl = TRUE, trunc = FALSE
  )
  modeltmp <- model

  for (i in seq_len(updates)) {
    tree <- trees[[tree_indices[i]]]
    invtree <- MCMCglmm::inverseA(tree, tol = inverse_tol)$Ainv
    if (!all(as.character(tiplabel) %in% rownames(invtree))) {
      stop("A sampled tree is missing analysis tiplabels.")
    }
    inverse_list <- setNames(list(invtree), phylo_effect)

    # The package expects a liability vector, not a nested Liab list.
    liab_start <- as.numeric(modeltmp$Liab[1, ])
    G_start <- list()
    for (j in seq_len(Gnumber)) {
      a <- if (j == 1L) 1L else 1L + sum(Geffects[seq_len(j - 1L)])
      b <- sum(Geffects[seq_len(j)])
      G_start[[paste0("G", j)]] <- if (Geffects[j] > 1L) {
        matrix(modeltmp$VCV[1, a:b], nrow = sqrt(Geffects[j]))
      } else {
        modeltmp$VCV[1, a:b]
      }
    }
    R_start <- list()
    glength <- sum(Geffects)
    for (j in seq_len(Rnumber)) {
      a <- glength + if (j == 1L) 1L else {
        1L + sum(Reffects[seq_len(j - 1L)])
      }
      b <- glength + sum(Reffects[seq_len(j)])
      R_start[[paste0("R", j)]] <- if (Reffects[j] > 1L) {
        matrix(modeltmp$VCV[1, a:b], nrow = sqrt(Reffects[j]))
      } else {
        modeltmp$VCV[1, a:b]
      }
    }
    start <- list(Liab = liab_start, G = G_start, R = R_start)

    modeltmp <- MCMCglmm::MCMCglmm(
      fixed = fixed, random = random, rcov = residual,
      ginverse = inverse_list, family = family, data = data, prior = prior,
      nitt = nitts_tree, thin = thins_tree, burnin = burns_tree,
      start = start, pr = pr, verbose = FALSE, pl = TRUE, trunc = FALSE
    )
    if (i > treetoburn) {
      model$VCV[i - treetoburn, ] <- modeltmp$VCV[1, ]
      model$Sol[i - treetoburn, ] <- modeltmp$Sol[1, ]
      model$Liab[i - treetoburn, ] <- modeltmp$Liab[1, ]
    }
  }
  model$phylogeny_sampling <- list(
    effect = phylo_effect, available_trees = length(trees),
    warmup_updates = treetoburn, retained_updates = samples,
    tree_indices = tree_indices,
    retained_tree_indices = tree_indices[treetoburn + seq_len(samples)],
    inverse_tolerance = inverse_tol
  )
  return(model)
}
