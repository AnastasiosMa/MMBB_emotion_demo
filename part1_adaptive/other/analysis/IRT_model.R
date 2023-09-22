#install.packages("mirt")

library(mirt)

# Extension for 'mirt' 
# devtools::install_github("masurp/ggmirt")
library(ggmirt)

library(readr)
binary_responses <- read_csv("Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/data/output/binary_responses/fear_binary_responses_only.csv")
View(binary_responses)

fitRasch <- mirt(binary_responses, 1, itemtype = "Rasch", verbose = T,guess = 0.5)
fitRasch

itemfit(fitRasch)


paramsRasch <- coef(fitRasch, IRTpars = TRUE, simplify = TRUE)
write.csv(paramsRasch, 'Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/data/output/binary_responses/irt_models/rasch_mirt.csv', row.names=FALSE)

summary(fitRasch)


tracePlot(fitRasch, theta_range = c(-5, 5), facet = F, legend = T) + 
  scale_color_brewer(palette = "Set3") +
  labs(title = "1PL - Traceplot")

fit3PL <- mirt(data = binary_responses, 
               model = 1,  # alternatively, we could also just specify model = 1 in this case
               itemtype = "3PL", 
               verbose = FALSE,guess = 0.5)

fit3PL
summary(fit3PL)

params3PL <- coef(fit3PL, IRTpars = TRUE, simplify = TRUE)
write.csv(params3PL, 'Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/data/output/binary_responses/irt_models/3pl_mirt.csv', row.names=FALSE)


itemfit(fit3PL)
itempersonMap(fit3PL)

fit2PL <- mirt(binary_responses, 1, itemtype = "2PL", verbose = F,guess = 0.5,method = 'EM',technical = list('NCYCLES',2000)) #technical = list("NCYCLES",1000)
fit2PL
params2PL <- coef(fit2PL, IRTpars = TRUE, simplify = TRUE)
write.csv(params2PL, 'Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/data/output/binary_responses/irt_models/2pl_mirt.csv', row.names=FALSE)

