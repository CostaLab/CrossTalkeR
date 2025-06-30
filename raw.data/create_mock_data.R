

CTR <- utils::read.csv("CTR_LR.csv")
EXP <- utils::read.csv("EXP_LR.csv")
ctkroutput <- readRDS("LR_data_final.Rds")


usethis::use_data(CTR, overwrite = TRUE)
usethis::use_data(EXP, overwrite = TRUE)
usethis::use_data(ctkroutput, overwrite = TRUE)
