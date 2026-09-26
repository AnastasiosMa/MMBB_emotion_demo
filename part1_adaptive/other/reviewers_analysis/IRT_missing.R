library(mirt)
library(readr)

setwd("Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/reviewers_analysis/")
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

# Fit models
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


#best fir for rasch and 1 dimension
mlScores <- fscores(fitRasch, method = 'ML')
mlScores_missing <- fscores(fitRasch_missing, method = 'ML')

paramsRasch <- coef(fitRasch, IRTpars = TRUE, simplify = TRUE)
paramsRasch_missing <- coef(fitRasch_missing, IRTpars = TRUE, simplify = TRUE)

difficulties<-data.frame(x=(scale(paramsRasch$items[,2])),y=(scale(paramsRasch_missing$items[,2])))
participant_scores<-data.frame(x=(mlScores),y=(mlScores_missing))

model <- lm(scale(paramsRasch$items[,2]) ~ scale(paramsRasch_missing$items[,2]))
summary(model)

#write.csv(paramsRasch, 'Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/data/output/binary_responses/irt_models/rasch_mirt.csv', row.names=FALSE)
#write.csv(infit_outfit_Rasch, 'Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/data/output/binary_responses/irt_models/rasch_infit_outfit.csv', row.names=FALSE)
#write.csv(participantScores, 'Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/data/output/binary_responses/irt_models/participantScores.csv', row.names=FALSE)

summary(fitRasch)
