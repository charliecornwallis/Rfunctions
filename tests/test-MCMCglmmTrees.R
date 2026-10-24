# Run from the repository root, or pass a helper path and evidence directory.
# Bounded simulated fits check computation, not convergence or inference.
library(ape)
library(MCMCglmm)
library(digest)
library(jsonlite)

arguments <- commandArgs(trailingOnly = TRUE)
helper_path <- if (length(arguments) > 0L) arguments[1] else {
  "MCMCglmmTrees.R"
}
evidence_dir <- if (length(arguments) > 1L) arguments[2] else NULL
invisible(parse(file = helper_path))
source(helper_path)
checks <- c()

# Existing positional arguments remain valid; additions are optional.
legacy_arguments <- c(
  "prior", "treetoburn", "data", "mev", "trees", "tiplabel", "fixed",
  "random", "residual", "Gnumber", "Rnumber", "Geffects", "Reffects",
  "family", "samples", "nitts_tree", "thins_tree", "burns_tree", "pr",
  "trunc"
)
checks["legacy_argument_order"] <- identical(
  head(names(formals(MCMCglmmTrees)), length(legacy_arguments)),
  legacy_arguments
)
checks["legacy_inverse_tolerance"] <- identical(
  eval(formals(MCMCglmmTrees)$inverse_tol), .Machine$double.eps^0.5
)

set.seed(2024)
test_tree <- rcoal(8)
test_tree$tip.label <- paste0("species_", seq_len(8))
test_data <- data.frame(
  animal = rep(test_tree$tip.label, each = 2), response = rnorm(16)
)
test_data$tiplabel <- test_data$animal
test_data$species <- test_data$animal
original_data <- test_data
model_arguments <- list(
  prior = list(R = list(V = 1, nu = 0.002),
               G = list(G1 = list(V = 1, nu = 0.002))),
  treetoburn = 2, data = test_data,
  trees = rep(list(test_tree), 6), tiplabel = test_data$animal,
  fixed = response ~ 1, random = ~ animal, residual = ~ units,
  Gnumber = 1, Rnumber = 1, Geffects = 1, Reffects = 1,
  samples = 4, nitts_tree = 20, thins_tree = 1, burns_tree = 19
)
test_models <- list()
for (effect in c("animal", "tiplabel", "species")) {
  fit_arguments <- model_arguments
  fit_arguments$random <- as.formula(paste("~", effect))
  if (effect == "species") {
    fit_arguments$phylo_effect <- effect
    fit_arguments$trees <- rep(list(test_tree), 2)
  }
  set.seed(73131)
  test_models[[effect]] <- do.call(MCMCglmmTrees, fit_arguments)
  fit <- test_models[[effect]]
  checks[paste0(effect, "_inverse_associated")] <-
    sum(fit$Random$nat) > 0 &&
    identical(fit$phylogeny_sampling$effect, effect)
  checks[paste0(effect, "_finite_retained_draws")] <-
    nrow(fit$Sol) == 4L && nrow(fit$VCV) == 4L &&
    nrow(fit$Liab) == 4L && all(is.finite(fit$Sol)) &&
    all(is.finite(fit$VCV)) && all(is.finite(fit$Liab))
}
checks["equivalent_factor_names_identical_draws"] <- identical(
  unname(as.matrix(test_models$animal$Sol)),
  unname(as.matrix(test_models$tiplabel$Sol))
) && identical(
  unname(as.matrix(test_models$animal$VCV)),
  unname(as.matrix(test_models$tiplabel$VCV))
)
checks["sufficient_tree_count_preserves_order"] <- identical(
  test_models$animal$phylogeny_sampling$tree_indices, seq_len(6)
)
checks["finite_tree_set_reused_equally"] <- identical(
  tabulate(test_models$species$phylogeny_sampling$tree_indices, nbins = 2),
  c(3L, 3L)
) && identical(
  tabulate(test_models$species$phylogeny_sampling$retained_tree_indices,
           nbins = 2), c(2L, 2L)
)
checks["caller_data_unchanged"] <- identical(test_data, original_data)

# Invalid bindings must fail before sampling rather than silently fit IID.
for (case in c("ambiguous", "missing_effect", "mismatched_labels",
               "missing_tree_tip")) {
  invalid_arguments <- model_arguments
  if (case == "ambiguous") {
    invalid_arguments$random <- ~ animal + tiplabel
  } else if (case == "missing_effect") {
    invalid_arguments$phylo_effect <- "not_in_formula"
  } else if (case == "mismatched_labels") {
    invalid_arguments$data$animal[1] <- "wrong_species"
  } else {
    invalid_arguments$trees <- rep(
      list(drop.tip(test_tree, test_tree$tip.label[1])), 6
    )
  }
  set.seed(73132)
  rng_before <- .Random.seed
  rejected <- tryCatch(
    do.call(MCMCglmmTrees, invalid_arguments), error = identity
  )
  checks[paste0(case, "_rejected_before_sampling")] <-
    inherits(rejected, "error") && identical(.Random.seed, rng_before)
}

stopifnot(all(checks))
if (!is.null(evidence_dir)) {
  dir.create(evidence_dir, recursive = TRUE, showWarnings = FALSE)
  write_json(
    list(
      date = "2026-10-24", review = "self-check",
      helper = normalizePath(helper_path),
      sha256 = digest(helper_path, algo = "sha256", file = TRUE),
      checks = as.list(checks), diagnostic_fits = length(test_models),
      production_fits = 0, inference_validated = FALSE
    ),
    file.path(evidence_dir, "verification.json"),
    pretty = TRUE, auto_unbox = TRUE
  )
  writeLines(
    trimws(capture.output(sessionInfo()), which = "right"),
    file.path(evidence_dir, "sessionInfo.txt")
  )
}
message(sum(checks), " checks passed; three bounded simulated fits.")
