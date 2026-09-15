library(synthpop)
library(readxl)
data<-read.csv('data/osf_materials/pilot_preprocessed.csv')
syn_data <- syn(data,method='cart',m=1)
syn<-syn_data$syn
compare(syn,data)
write.csv(syn,'pilot_synthetic_data.csv')