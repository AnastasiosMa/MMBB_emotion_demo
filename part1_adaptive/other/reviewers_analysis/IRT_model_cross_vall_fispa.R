library(mirt)
library(readr)

set.seed(20260930)

data_dir <- "data"
output_dir <- file.path(data_dir, "cross_validation_fi_spa")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_group <- function(group) {
  response_path <- file.path(data_dir, paste0("binary_responses_", group, ".csv"))
  info_path <- file.path(data_dir, paste0("trial_info_", group, ".csv"))

  if (!file.exists(response_path) || !file.exists(info_path)) {
    stop("Could not find group files under ", data_dir,)}

  responses <- read_csv(response_path, show_col_types = FALSE)
  source_item_names <- setdiff(names(responses), c("Participant", "Age", "Gender"))
  info <- read_csv(info_path, show_col_types = FALSE)

  if (length(source_item_names) != nrow(info)) {
    stop(group, ": item columns do not match trial-info rows")
  }

  item_names <- paste0(
    "excerpt_", info$Track110,
    "_target_", info$TargetEmo,
    "_comparison_", info$ComparisonEmo
  )
  item_matrix <- as.matrix(responses[, source_item_names, drop = FALSE])
  storage.mode(item_matrix) <- "double"
  colnames(item_matrix) <- item_names

  if (ncol(item_matrix) != nrow(info)) {
    stop(group, ": item columns do not match trial-info rows")
  }

  if (!all(is.na(item_matrix) | item_matrix %in% c(0, 1))) {
    stop(group, ": response matrix contains values other than 0, 1, or NA")
  }

  list(responses = item_matrix, item_names = item_names, info = info)
}

fi <- read_group("fi")
spa <- read_group("spa")

if (!identical(fi$item_names, spa$item_names)) {
  stop("Finnish and Spanish item columns are not in the same order")
}

item_key <- c("Labels", "TargetEmo", "ComparisonEmo", "Track110")
if (!identical(fi$info[item_key], spa$info[item_key])) {
  stop("Finnish and Spanish trial metadata are not aligned")
}

estimate_eap <- function(responses, item_indices, difficulties,
                         prior_mean, prior_sd, theta_grid) {
  estimates <- rep(NA_real_, nrow(responses))

  for (person in seq_len(nrow(responses))) {
    observed <- !is.na(responses[person, item_indices])
    if (!any(observed)) {
      next
    }

    person_responses <- responses[person, item_indices[observed]]
    person_difficulties <- difficulties[item_indices[observed]]
    log_likelihood <- numeric(length(theta_grid))

    for (item in seq_along(person_responses)) {
      probability <- plogis(theta_grid - person_difficulties[item])
      probability <- pmin(pmax(probability, 1e-12), 1 - 1e-12)
      if (person_responses[item] == 1) {
        log_likelihood <- log_likelihood + log(probability)
      } else {
        log_likelihood <- log_likelihood + log1p(-probability)
      }
    }

    log_posterior <- log_likelihood +
      dnorm(theta_grid, mean = prior_mean, sd = prior_sd, log = TRUE)
    posterior_weights <- exp(log_posterior - max(log_posterior))
    estimates[person] <- sum(theta_grid * posterior_weights) /
      sum(posterior_weights)
  }

  estimates
}

select_excerpt_split <- function(info, difficulties, proportion = 0.8,
                                 seed, attempts = 20000) {
  item_count <- length(difficulties)
  difficulty_order <- order(difficulties)
  difficulty_stratum <- character(item_count)
  difficulty_stratum[difficulty_order] <- as.character(cut(
    seq_len(item_count),
    breaks = c(0, item_count / 3, 2 * item_count / 3, item_count),
    labels = c("Low", "Mid", "High"),
    include.lowest = TRUE
  ))

  excerpt_id <- as.character(info$Track110)
  excerpts <- unique(excerpt_id)
  stratum_levels <- c("Low", "Mid", "High")
  excerpt_index <- match(excerpt_id, excerpts)
  stratum_index <- match(difficulty_stratum, stratum_levels)
  profile <- matrix(0, nrow = length(excerpts), ncol = length(stratum_levels),
                    dimnames = list(excerpts, stratum_levels))

  for (item in seq_len(item_count)) {
    profile[excerpt_index[item], stratum_index[item]] <-
      profile[excerpt_index[item], stratum_index[item]] + 1
  }

  target_by_stratum <- proportion * colSums(profile)
  target_total <- proportion * item_count
  scale_by_stratum <- pmax(target_by_stratum, 1)
  best_score <- Inf
  best_selected <- NULL

  set.seed(seed)
  for (attempt in seq_len(attempts)) {
    selected_excerpts <- runif(length(excerpts)) < proportion
    selected_counts <- colSums(profile[selected_excerpts, , drop = FALSE])
    heldout_counts <- colSums(profile[!selected_excerpts, , drop = FALSE])

    score <- sum(((selected_counts - target_by_stratum) /
                    scale_by_stratum)^2) +
      ((sum(selected_counts) - target_total) / target_total)^2

    if (any(selected_counts == 0) || any(heldout_counts == 0)) {
      score <- score + 1e6
    }

    if (score < best_score) {
      best_score <- score
      best_selected <- selected_excerpts
    }
  }

  selected_item <- excerpt_id %in% excerpts[best_selected]

  list(
    calibration_items = which(selected_item),
    heldout_items = which(!selected_item),
    difficulty_stratum = difficulty_stratum,
    selected_excerpts = excerpts[best_selected]
  )
}

run_direction <- function(training, testing, direction, split_seed) {
  initial_training_fit <- mirt(
    training$responses,
    1,
    itemtype = "Rasch",
    verbose = FALSE
  )
  fit_statistics <- itemfit(initial_training_fit, fit_stats = "infit")
  fit_rows <- match(training$item_names, fit_statistics$item)
  infit <- fit_statistics$infit[fit_rows]
  outfit <- fit_statistics$outfit[fit_rows]
  finite_fit <- is.finite(infit) & is.finite(outfit)
  acceptable_fit <- finite_fit & infit >= 0.5 & infit <= 1.5 &
    outfit >= 0.5 & outfit <= 1.5

  fit_screen <- data.frame(
    direction = direction,
    item = training$item_names,
    excerpt = training$info$Track110,
    infit = infit,
    outfit = outfit,
    retained = acceptable_fit,
    exclusion_reason = ifelse(
      !finite_fit,
      "nonfinite_infit_or_outfit",
      ifelse(acceptable_fit, "retained", "outside_0.5_to_1.5")
    )
  )

  if (!any(acceptable_fit)) {
    stop(direction, ": no items passed the infit/outfit screen")
  }

  training$responses <- training$responses[, acceptable_fit, drop = FALSE]
  testing$responses <- testing$responses[, acceptable_fit, drop = FALSE]
  training$item_names <- training$item_names[acceptable_fit]
  testing$item_names <- testing$item_names[acceptable_fit]
  training$info <- training$info[acceptable_fit, , drop = FALSE]
  testing$info <- testing$info[acceptable_fit, , drop = FALSE]

  training_fit <- mirt(
    training$responses,
    1,
    itemtype = "Rasch",
    verbose = FALSE
  )

  item_parameters <- coef(training_fit, IRTpars = TRUE, simplify = TRUE)$items
  difficulties <- item_parameters[training$item_names, "b"]
  latent_parameters <- coef(training_fit, simplify = TRUE)
  prior_mean <- as.numeric(latent_parameters$means[1])
  prior_sd <- sqrt(as.numeric(latent_parameters$cov[1, 1]))

  split <- select_excerpt_split(
    training$info,
    difficulties,
    proportion = 0.8,
    seed = split_seed
  )
  calibration_items <- split$calibration_items
  heldout_items <- split$heldout_items

  theta_grid <- seq(prior_mean - 6 * prior_sd,
                    prior_mean + 6 * prior_sd,
                    length.out = 1201)
  ability <- estimate_eap(
    testing$responses,
    calibration_items,
    difficulties,
    prior_mean,
    prior_sd,
    theta_grid
  )

  predictions <- vector("list", nrow(testing$responses))
  person_metrics <- vector("list", nrow(testing$responses))

  for (person in seq_len(nrow(testing$responses))) {
    heldout_response <- testing$responses[person, heldout_items]
    observed <- !is.na(heldout_response) & !is.na(ability[person])
    if (!any(observed)) {
      next
    }

    observed_items <- heldout_items[observed]
    actual <- heldout_response[observed]
    predicted_probability <- plogis(
      ability[person] - difficulties[observed_items]
    )
    predicted_class <- as.integer(predicted_probability >= 0.5)
    clipped_probability <- pmin(pmax(predicted_probability, 1e-12), 1 - 1e-12)

    predictions[[person]] <- data.frame(
      direction = direction,
      person_row = person,
      item = testing$item_names[observed_items],
      excerpt = testing$info$Track110[observed_items],
      difficulty_stratum = split$difficulty_stratum[observed_items],
      actual = actual,
      predicted_probability = predicted_probability,
      predicted_class = predicted_class
    )

    person_metrics[[person]] <- data.frame(
      direction = direction,
      person_row = person,
      theta_eap = ability[person],
      n_calibration_responses = sum(
        !is.na(testing$responses[person, calibration_items])
      ),
      n_heldout_responses = length(actual),
      accuracy = mean(predicted_class == actual),
      brier = mean((predicted_probability - actual)^2),
      log_loss = mean(-actual * log(clipped_probability) -
                        (1 - actual) * log1p(-clipped_probability))
    )
  }

  prediction_table <- do.call(
    rbind, predictions[!vapply(predictions, is.null, logical(1))]
  )
  person_table <- do.call(
    rbind, person_metrics[!vapply(person_metrics, is.null, logical(1))]
  )

  split_table <- data.frame(
    direction = direction,
    item = training$item_names,
    excerpt = training$info$Track110,
    training_difficulty = as.numeric(difficulties),
    difficulty_stratum = split$difficulty_stratum,
    role = ifelse(seq_along(difficulties) %in% calibration_items,
                  "ability_estimation", "heldout_prediction")
  )

  heldout_actual <- prediction_table$actual
  clipped_all <- pmin(pmax(prediction_table$predicted_probability, 1e-12),
                      1 - 1e-12)
  pooled_summary <- data.frame(
    direction = direction,
    training_group = training$group,
    testing_group = testing$group,
    n_training_people = nrow(training$responses),
    n_testing_people = nrow(testing$responses),
    n_items = length(difficulties),
    n_items_before_fit_screen = nrow(fit_screen),
    n_items_excluded_fit_screen = sum(!fit_screen$retained),
    n_calibration_items = length(calibration_items),
    calibration_item_percent = 100 * length(calibration_items) /
      length(difficulties),
    n_heldout_items = length(heldout_items),
    heldout_item_percent = 100 * length(heldout_items) /
      length(difficulties),
    n_calibration_excerpts = length(split$selected_excerpts),
    n_heldout_responses = nrow(prediction_table),
    heldout_accuracy = mean(prediction_table$predicted_class == heldout_actual),
    heldout_brier = mean((prediction_table$predicted_probability - heldout_actual)^2),
    heldout_log_loss = mean(
      -heldout_actual * log(clipped_all) -
        (1 - heldout_actual) * log1p(-clipped_all)
    ),
    mean_person_accuracy = mean(person_table$accuracy),
    training_theta_mean = prior_mean,
    training_theta_sd = prior_sd
  )

  list(
    summary = pooled_summary,
    predictions = prediction_table,
    person_metrics = person_table,
    item_split = split_table,
    item_fit = fit_screen
  )
}

fi$group <- "FI"
spa$group <- "SPA"

fi_to_spa <- run_direction(fi, spa, "FI_to_SPA", split_seed = 20260931)
spa_to_fi <- run_direction(spa, fi, "SPA_to_FI", split_seed = 20260932)

write_csv(
  rbind(fi_to_spa$summary, spa_to_fi$summary),
  file.path(output_dir, "summary.csv")
)
write_csv(
  rbind(fi_to_spa$item_split, spa_to_fi$item_split),
  file.path(output_dir, "item_split.csv")
)
write_csv(
  rbind(fi_to_spa$item_fit, spa_to_fi$item_fit),
  file.path(output_dir, "item_fit_screen.csv")
)
write_csv(
  rbind(fi_to_spa$person_metrics, spa_to_fi$person_metrics),
  file.path(output_dir, "person_metrics.csv")
)
write_csv(
  rbind(fi_to_spa$predictions, spa_to_fi$predictions),
  file.path(output_dir, "heldout_predictions.csv")
)

print(rbind(fi_to_spa$summary, spa_to_fi$summary))