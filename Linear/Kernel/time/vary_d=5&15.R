setwd("C:/Users/smm/Desktop/DC-code/Setting2-Linear/Kernel/time")
source("bws_d=5_time.R")
save.image("C:/Users/smm/Desktop/DC-code/Setting2-Linear/Kernel/time/linear_bws_d=5_time.RData")
rm(list = ls())

source("bws_d=15_time.R")
save.image("C:/Users/smm/Desktop/DC-code/Setting2-Linear/Kernel/time/linear_bws_d=15_time.RData")
rm(list = ls())

source("abla_wht_d=5_time.R")
save.image("C:/Users/smm/Desktop/DC-code/Setting2-Linear/Kernel/time/linear_abla_wht_d=5_time.RData")
rm(list = ls())

source("abla_wht_d=15_time.R")
save.image("C:/Users/smm/Desktop/DC-code/Setting2-Linear/Kernel/time/linear_abla_wht_d=15_time.RData")
rm(list = ls())
