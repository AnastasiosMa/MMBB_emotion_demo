library(mirt)
library(readr)

setwd("Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/")
binary_responses <- read_csv("reviewers_analysis/data/binary_responses.csv")
binary_responses_missing <- read_csv("reviewers_analysis/data/binary_responses_missing.csv")

# Train IRT model
fitRasch_fi <- mirt(binary_responses_fi, 1, itemtype = 'Rasch', verbose = T,guess = 0.5)
print(fitRasch_fi)

fitRasch_spa <- mirt(binary_responses_spa, 1, itemtype = 'Rasch', verbose = T,guess = 0.5)


#best fir for rasch and 1 dimension
mlScores_spa <- fscores(fitRasch_spa, method = 'ML')
mlScores_fi <- fscores(fitRasch_fi, method = 'ML')

paramsRasch_spa <- coef(fitRasch_spa, IRTpars = TRUE, simplify = TRUE)
paramsRasch_fi <- coef(fitRasch_fi, IRTpars = TRUE, simplify = TRUE)

difficulties<-data.frame(x=(scale(paramsRasch_spa$items[,2])),y=(scale(paramsRasch_fi$items[,2])))

model <- lm(scale(paramsRasch_spa$items[,2]) ~ scale(paramsRasch_fi$items[,2]))
summary(model)

write.csv(paramsRasch, 'Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/data/output/binary_responses/irt_models/rasch_mirt.csv', row.names=FALSE)
write.csv(infit_outfit_Rasch, 'Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/data/output/binary_responses/irt_models/rasch_infit_outfit.csv', row.names=FALSE)
write.csv(participantScores, 'Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/data/output/binary_responses/irt_models/participantScores.csv', row.names=FALSE)

summary(fitRasch)
