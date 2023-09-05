#install.packages("mirt")

library(mirt)

# Extension for 'mirt' 
# devtools::install_github("masurp/ggmirt")
library(ggmirt)

library(readr)
binary_responses <- read_csv("Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/data/output/binary_responses/binary_responses.csv")
View(binary_responses)

fitRasch <- mirt(binary_responses, 1, itemtype = "Rasch", verbose = T)

coef(fitRasch, IRTpars = TRUE, simplify = TRUE)

fitRasch

summary(fitRasch)


tracePlot(fitRasch, theta_range = c(-5, 5), facet = F, legend = T) + 
  scale_color_brewer(palette = "Set3") +
  labs(title = "1PL - Traceplot")

fit3PL <- mirt(data = d, 
               model = 1,  # alternatively, we could also just specify model = 1 in this case
               itemtype = "3PL", 
               verbose = FALSE)

fit3PL
summary(fit3PL)

params3PL <- coef(fit3PL, IRTpars = TRUE, simplify = TRUE)
round(params3PL$items, 2)

itempersonMap(fit3PL)

