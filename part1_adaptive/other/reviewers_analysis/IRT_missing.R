library(mirt)
library(readr)

setwd("Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/")
binary_responses <- read_csv("reviewers_analysis/data/binary_responses_all.csv")
binary_responses=binary_responses[,1:(ncol(binary_responses)-3)]
binary_responses_missing <- read_csv("reviewers_analysis/data/binary_responses_all_missing.csv")
binary_responses_missing=binary_responses_missing[,1:(ncol(binary_responses_missing)-3)]


# Find items with only one observed response category
bad_items <- sapply(binary_responses_missing, function(x) {
  vals <- unique(na.omit(x))
  length(vals) < 2})

# Names of problematic items
bad_item_names <- names(binary_responses_missing)[bad_items]

# Print them
bad_item_names

# Remove from BOTH datasets
binary_responses_missing_filt <- binary_responses_missing[, !bad_items]
binary_responses_filt <- binary_responses[, !bad_items]

# Check dimensions
dim(binary_responses)
dim(binary_responses_filt)

# Initial fits
fitRasch <- mirt(binary_responses_filt,
                 1,
                 itemtype = "Rasch",
                 verbose = TRUE,
                 guess = 0.5)

fitRasch_missing <- mirt(binary_responses_missing_filt,
                         1,
                         itemtype = "Rasch",
                         verbose = TRUE,
                         guess = 0.5)

paramsRasch <- coef(fitRasch, IRTpars = TRUE, simplify = TRUE)
paramsRasch_missing <- coef(fitRasch_missing, IRTpars = TRUE, simplify = TRUE)

b1 <- as.vector(paramsRasch$items[,2])
b2 <- as.vector(paramsRasch_missing$items[,2])
cat(sprintf("Difficulty correlation: r = %.3f; R-squared = %.3f\n",
            cor(b1, b2),
            cor(b1, b2)^2))
plot(paramsRasch$items[,2],
     paramsRasch_missing$items[,2])
abline(0,1,col="red")

# Extract difficulties
diff1 <- coef(fitRasch, IRTpars = TRUE, simplify = TRUE)$items[,2]
diff2 <- coef(fitRasch_missing, IRTpars = TRUE, simplify = TRUE)$items[,2]

# Find items with extreme difficulty
bad_diff <- (diff1 > 4) | (diff2 > 4)

# Item names
bad_diff_names <- names(diff1)[bad_diff]
bad_diff_names_missing <- names(diff2)[bad_diff]

bad_diff_names

# Remove from both datasets
binary_responses_filt2 <- binary_responses_filt[, !colnames(binary_responses_filt) %in% bad_diff_names]
binary_responses_missing_filt2 <- binary_responses_missing_filt[, !colnames(binary_responses_missing_filt) %in% bad_diff_names_missing]

# Item fit statistics
fit_stats <- itemfit(fitRasch,
                     fit_stats = "infit")

fit_stats_missing <- itemfit(fitRasch_missing,
                             fit_stats = "infit")

# Flag items with infit or outfit outside acceptable range
bad_fit <- (
  fit_stats$infit < 0.5 |
    fit_stats$infit > 1.5 |
    fit_stats$outfit < 0.5 |
    fit_stats$outfit > 1.5 |
    fit_stats_missing$infit < 0.5 |
    fit_stats_missing$infit > 1.5 |
    fit_stats_missing$outfit < 0.5 |
    fit_stats_missing$outfit > 1.5
)

# Item names
bad_fit_names <- fit_stats$item[bad_fit]
bad_diff_names_missing <- fit_stats_missing$item[bad_fit]

bad_fit_names

# Remove from both datasets
binary_responses_filt2 <-
  binary_responses_filt2[, !colnames(binary_responses_filt2) %in% bad_fit_names]

binary_responses_missing_filt2 <-
  binary_responses_missing_filt2[, !colnames(binary_responses_missing_filt2) %in% bad_diff_names_missing]

paramsRasch <- coef(fitRasch, IRTpars = TRUE, simplify = TRUE)
paramsRasch_missing <- coef(fitRasch_missing, IRTpars = TRUE, simplify = TRUE)

b1 <- as.vector(paramsRasch$items[,2])
b2 <- as.vector(paramsRasch_missing$items[,2])
cat(sprintf("Difficulty correlation: r = %.3f; R-squared = %.3f\n",
            cor(b1, b2),
            cor(b1, b2)^2))

# Refit
fitRasch <- mirt(binary_responses_filt2,
                 1,
                 itemtype = "Rasch",
                 verbose = TRUE, guess = 0.5)

fitRasch_missing <- mirt(binary_responses_missing_filt2,
                         1,
                         itemtype = "Rasch",
                         verbose = TRUE,guess = 0.5)

paramsRasch <- coef(fitRasch, IRTpars = TRUE, simplify = TRUE)
paramsRasch_missing <- coef(fitRasch_missing, IRTpars = TRUE, simplify = TRUE)

b1 <- as.vector(paramsRasch$items[,2])
b2 <- as.vector(paramsRasch_missing$items[,2])
cat(sprintf("Difficulty correlation: r = %.3f; R-squared = %.3f\n",
            cor(b1, b2),
            cor(b1, b2)^2))

plot(paramsRasch$items[,2],
     paramsRasch_missing$items[,2])
abline(0,1,col="red")

# ---- Percentage of ties effect
# Put both difficulty vectors in the same explicit item order
item_names <- names(binary_responses_filt2)
item_names_missing <- names(binary_responses_missing_filt2)

# Tie rate among originally observed responses, separately for each item
main_matrix <- as.matrix(binary_responses_filt2[, item_names, drop = FALSE])
missing_matrix <- as.matrix(binary_responses_missing_filt2[, item_names_missing, drop = FALSE])

tie_rate <- colSums(!is.na(main_matrix) & is.na(missing_matrix)) /
  colSums(!is.na(main_matrix))

# Item-level difficulty change: negative means lower b in missing condition
difficulty_change <- b2 - b1

item_results <- data.frame(
  item = item_names,
  tie_rate = tie_rate,
  b_main = b1,
  b_missing = b2,
  difficulty_change = difficulty_change
)

cat("Items retained:", nrow(item_results), "\n")
cat(sprintf("Tie rate: min %.1f%%, mean %.1f%%,std %.1f%%,median %.1f%%, max %.1f%%\n",
            100 * min(tie_rate),
            100 * mean(tie_rate),
            100 * sd(tie_rate),
            100 * median(tie_rate),
            100 * max(tie_rate)))
cat(sprintf("Mean difficulty change: %.3f\n",
            mean(difficulty_change)))
cat(sprintf("Correlation of tie rate with difficulty change: r = %.3f\n",
            cor(tie_rate, difficulty_change)))
cat(sprintf("Correlation of initial difficulty and difficulty change: r = %.3f\n",
            cor(b1, difficulty_change)))
cat(sprintf("Difficulty correlation: r = %.3f; R-squared = %.3f\n",
            cor(b1, b2),
            cor(b1, b2)^2))

# Inspect the items with the highest tie rates
print(item_results[order(-item_results$tie_rate), ][1:10, ])

# Plot tie rate against difficulty change
plot(100 * tie_rate, difficulty_change,
     xlab = "Ties among observed responses (%)",
     ylab = "Difficulty change (missing - main)")
abline(lm(difficulty_change ~ tie_rate), col = "red")

difficulties<-data.frame(x=(scale(paramsRasch$items[,2])),y=(scale(paramsRasch_missing$items[,2])))
participant_scores<-data.frame(x=(mlScores),y=(mlScores_missing))

model <- lm(paramsRasch$items[,2] ~ paramsRasch_missing$items[,2])
summary(model)
