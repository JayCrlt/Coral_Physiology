# Second_part Diego
rm(list=ls())
Results_directory = paste(getwd(),"/Results",sep="") ; Data_directory = paste(getwd(),"/Data",sep="")
setwd(Data_directory) ; options(mc.cores = parallel::detectCores())
library('tidyverse') ; library('parallel') ; library('brms') ; library('abind') ; library('readxl')
library('patchwork') ; library("tidybayes") ; library('ggpubr') ; library('rstatix')
Data_Raw_Metabo <- read_excel("Data_Metabo.xlsx")
Data_Raw_Metabo = Data_Raw_Metabo %>% mutate(Calcif_kg_yr = Corrected_Production/100,
                                             log_area = log(Surface_area_cm2),
                                             log_calcif = log(Calcif_kg_yr),
                                             log_respi = log(NightValue),
                                             log_photo = log(NightValue-Dayvalue)) %>% 
  mutate(Species = fct_recode(Species, 
                              "Acropora hyacinthus" = "Acropora_hyacinthus", "Montipora verilli" = "Montipora_verilli",
                              "Astrea curta" = "Astrea_curta", "Pocillopora verrucosa" = "Pocillopora_verrucosa",
                              "Porites lutea" = "Porites_lutea", "Napopora irregularis" = "Napopora_irregularis")) 

if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() == "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")}

load("Models/mixt_calcif_PB.RData") ; load("Models/mixt_calcif.RData") ; load("Models/mixt_respi_PB.RData")
load("Models/mixt_respi.RData") ; load("Models/mixt_photo_PB.RData") ; load("Models/mixt_photo.RData")

data_1 = data.frame(Species = rep(Data_Raw_Metabo$Species,4),
         Log_Area = rep(Data_Raw_Metabo$log_area,4),
         Log_Metabolism_Rate = c(predict(mixt_calcif, Data_Raw_Metabo)[,1],
                                 predict(mixt_respi, Data_Raw_Metabo)[,1],
                                 predict(mixt_photo, Data_Raw_Metabo)[,1],
                                 predict(mixt_photo, Data_Raw_Metabo)[,1] - predict(mixt_respi, Data_Raw_Metabo)[,1]),
         Function = rep(c("1) Calcifcation","2) Respiration","3) Photosynthesis","4) Net Primary Production"), each = 250))
data_2 = data.frame(Species = Data_Raw_Metabo$Species,
                    Surface_Area = Data_Raw_Metabo$Surface_area_cm2/10000,
                    Log_Calcification = predict(mixt_calcif, Data_Raw_Metabo)[,1],
                    Log_NPP = c(predict(mixt_photo, Data_Raw_Metabo)[,1])) # - predict(mixt_respi, Data_Raw_Metabo)[,1]))

Fig_1_Diego = ggplot(data_1) + facet_wrap(~Function, ncol = 4) + 
  geom_point(aes(x = Log_Area, y = Log_Metabolism_Rate, col = Species)) +
  scale_color_viridis_d() + theme_classic() + scale_x_continuous(name = expression("Log(Area) ("*cm^2*")")) +
  scale_y_continuous(name = expression("Log(Metabolic Rate) (Units)")) 
Fig_2_Diego = ggplot(data_2, aes(x = ((exp(Log_NPP)*24/1000)/14.003)*1000, 
                                 y = ((exp(Log_Calcification)/365*1000)/100.087)*1000, col = Species)) + 
  geom_point() + scale_x_continuous(name = expression("(Net Primary Production) ("*mmol.day^-1*")")) +
  scale_color_viridis_d() + theme_classic() + scale_y_continuous(name = expression("(Calcification) ("*mmol.day^-1*")")) +
  geom_smooth(method = "nls", formula=y~a*x^b, method.args = list(start=c(a=1,b=1)), se=FALSE)
Fig_1_Diego / Fig_2_Diego + plot_layout(guides = "collect")

data_3 = data.frame(NPP = ((exp(data_2$Log_NPP)*24/1000)/14.003)*1000,
                    Calcification = ((exp(data_2$Log_Calcification)/365*1000)/100.087)*1000,
           Ratio = (((exp(data_2$Log_NPP)*24/1000)/14.003)*1000)/(((exp(data_2$Log_Calcification)/365*1000)/100.087)*1000),
           Species = Data_Raw_Metabo$Species,
           Surface_Area = Data_Raw_Metabo$Surface_area_cm2)

h1 = data_3 %>% filter(Species == "Acropora hyacinthus") %>% ggplot(aes(x = Ratio)) + 
  geom_density(fill="#44015400", color="black", alpha=0.85) + theme_classic()
h2 = data_3 %>% filter(Species == "Astrea curta") %>% ggplot(aes(x = Ratio)) + 
  geom_density(fill="#414487FF", color="black", alpha=0.85) + theme_classic()
h3 = data_3 %>% filter(Species == "Montipora verilli") %>% ggplot(aes(x = Ratio)) + 
  geom_density(fill="#2A788EFF", color="black", alpha=0.85) + theme_classic()
h4 = data_3 %>% filter(Species == "Napopora irregularis") %>% ggplot(aes(x = Ratio)) + 
  geom_density(fill="#22A884FF", color="black", alpha=0.85) + theme_classic()
h5 = data_3 %>% filter(Species == "Pocillopora verrucosa") %>% ggplot(aes(x = Ratio)) + 
  geom_density(fill="#7AD151FF", color="black", alpha=0.85) + theme_classic()
h6 = data_3 %>% filter(Species == "Porites lutea") %>% ggplot(aes(x = Ratio)) + 
  geom_density(fill="#FDE725FF", color="black", alpha=0.85) + theme_classic()

Density_Plot = h1|h2|h3|h4|h5|h6

data_ACR = sample(rnorm(1000,mean=300,sd=200)) ; data_POC = data_ACR
Area = data.frame(ID = seq(1,1000,1), Surface_Area_ACR = data_ACR, Surface_Area_POC = data_POC) %>% 
  filter(Surface_Area_ACR > 0, Surface_Area_POC > 0)

ACR_Juv = draw_sample(dist=Area[,c(1,2)], n=50, skew = 0.8, kurts = 0, replacement = T ,output_name = c("sample", "1"))
ACR_Bal = draw_sample(dist=Area[,c(1,2)], n=50, skew = 0.5, kurts = 0, replacement = T ,output_name = c("sample", "1"))
ACR_Adu = draw_sample(dist=Area[,c(1,2)], n=50, skew = 0.2, kurts = 0, replacement = T ,output_name = c("sample", "1"))
POC_Juv = ACR_Juv ; POC_Bal = ACR_Bal ; POC_Adu = ACR_Adu

#####################
### Matrix panels ###
#####################

### ACROPORA HYACINTHUS ###
# Combination between Acropora and Pocillopora
data_4 = data_3 %>% filter(Species %in% c("Acropora hyacinthus","Pocillopora verrucosa")) 
# Build a new dataset according to the random sampling with the 3 skewness conditions
data = data.frame(Dataset = c(rep("Juv", 100), rep("Bal", 100), rep("Adu", 100)),
                  Surface_Area = c(sort(ACR_Juv$sample$x),sort(POC_Juv$sample$x),sort(ACR_Bal$sample$x),
                                   sort(POC_Bal$sample$x),sort(ACR_Adu$sample$x),sort(POC_Adu$sample$x)),
                  Species = rep(rep(c("Acropora hyacinthus","Pocillopora verrucosa"), each = 50),3))  %>% 
  mutate(log_area = log(Surface_Area)) %>% group_split(Species)
## ADULTS ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Adu")
data_ACR = data_ACR %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Adu")
data_POC = data_POC %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_adu = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_adu = data.frame(expand.grid(n_col, n_col), comb_adu) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_adu)

## BALANCED ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Bal")
data_ACR = data_ACR %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Bal")
data_POC = data_POC %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_bal = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_bal = data.frame(expand.grid(n_col, n_col), comb_bal) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_bal)

## JUVENILES ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Juv")
data_ACR = data_ACR %>% mutate(Calcif =((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS =((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Juv")
data_POC = data_POC %>% mutate(Calcif =((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS =((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_juv = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_juv = data.frame(expand.grid(n_col, n_col), comb_juv) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_juv)
# Dataset Acropora - Pocillopora
Comb_dataset_ACR = rbind(comb_juv, comb_bal, comb_adu) %>% as.data.frame() %>% 
  mutate(Sampling = rep(c("1) Juveniles", "2) Balance", "3) Adults"), each = 2500), 
         ID = rep(seq(1,2500,1),3), Pocillopora = Pocillopora *2, Acropora = Acropora *2,
         Sum = Acropora + Pocillopora) %>% filter(Sum <= 100)

###MONTIPORA VERILLI ###
# Combination between Acropora and Pocillopora
data_4 = data_3 %>% filter(Species %in% c("Montipora verilli","Pocillopora verrucosa")) 
# Build a new dataset according to the random sampling with the 3 skewness conditions
data = data.frame(Dataset = c(rep("Juv", 100), rep("Bal", 100), rep("Adu", 100)),
                  Surface_Area = c(sort(ACR_Juv$sample$x),sort(POC_Juv$sample$x),sort(ACR_Bal$sample$x),
                                   sort(POC_Bal$sample$x),sort(ACR_Adu$sample$x),sort(POC_Adu$sample$x)),
                  Species = rep(rep(c("Montipora verilli","Pocillopora verrucosa"), each = 50),3))  %>% 
  mutate(log_area = log(Surface_Area)) %>% group_split(Species)
## ADULTS ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Adu")
data_ACR = data_ACR %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Adu")
data_POC = data_POC %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_adu = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_adu = data.frame(expand.grid(n_col, n_col), comb_adu) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_adu)

## BALANCED ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Bal")
data_ACR = data_ACR %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Bal")
data_POC = data_POC %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_bal = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_bal = data.frame(expand.grid(n_col, n_col), comb_bal) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_bal)

## JUVENILES ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Juv")
data_ACR = data_ACR %>% mutate(Calcif =((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS =((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Juv")
data_POC = data_POC %>% mutate(Calcif =((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS =((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_juv = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_juv = data.frame(expand.grid(n_col, n_col), comb_juv) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_juv)
# Dataset Acropora - Pocillopora
Comb_dataset_MON = rbind(comb_juv, comb_bal, comb_adu) %>% as.data.frame() %>% 
  mutate(Sampling = rep(c("1) Juveniles", "2) Balance", "3) Adults"), each = 2500), 
         ID = rep(seq(1,2500,1),3), Pocillopora = Pocillopora *2, Acropora = Acropora *2,
         Sum = Acropora + Pocillopora) %>% filter(Sum <= 100)

### PORITES LUTEA ###
# Combination between Acropora and Pocillopora
data_4 = data_3 %>% filter(Species %in% c("Porites lutea","Pocillopora verrucosa")) 
# Build a new dataset according to the random sampling with the 3 skewness conditions
data = data.frame(Dataset = c(rep("Juv", 100), rep("Bal", 100), rep("Adu", 100)),
                  Surface_Area = c(sort(ACR_Juv$sample$x),sort(POC_Juv$sample$x),sort(ACR_Bal$sample$x),
                                   sort(POC_Bal$sample$x),sort(ACR_Adu$sample$x),sort(POC_Adu$sample$x)),
                  Species = rep(rep(c("Porites lutea","Pocillopora verrucosa"), each = 50),3))  %>% 
  mutate(log_area = log(Surface_Area)) %>% group_split(Species)
## ADULTS ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Adu")
data_ACR = data_ACR %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Adu")
data_POC = data_POC %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_adu = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_adu = data.frame(expand.grid(n_col, n_col), comb_adu) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_adu)

## BALANCED ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Bal")
data_ACR = data_ACR %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Bal")
data_POC = data_POC %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_bal = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_bal = data.frame(expand.grid(n_col, n_col), comb_bal) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_bal)

## JUVENILES ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Juv")
data_ACR = data_ACR %>% mutate(Calcif =((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS =((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Juv")
data_POC = data_POC %>% mutate(Calcif =((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS =((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_juv = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_juv = data.frame(expand.grid(n_col, n_col), comb_juv) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_juv)
# Dataset Acropora - Pocillopora
Comb_dataset_POR = rbind(comb_juv, comb_bal, comb_adu) %>% as.data.frame() %>% 
  mutate(Sampling = rep(c("1) Juveniles", "2) Balance", "3) Adults"), each = 2500), 
         ID = rep(seq(1,2500,1),3), Pocillopora = Pocillopora *2, Acropora = Acropora *2,
         Sum = Acropora + Pocillopora) %>% filter(Sum <= 100)

### ASTREA CURTA ###
# Combination between Acropora and Pocillopora
data_4 = data_3 %>% filter(Species %in% c("Astrea curta","Pocillopora verrucosa")) 
# Build a new dataset according to the random sampling with the 3 skewness conditions
data = data.frame(Dataset = c(rep("Juv", 100), rep("Bal", 100), rep("Adu", 100)),
                  Surface_Area = c(sort(ACR_Juv$sample$x),sort(POC_Juv$sample$x),sort(ACR_Bal$sample$x),
                                   sort(POC_Bal$sample$x),sort(ACR_Adu$sample$x),sort(POC_Adu$sample$x)),
                  Species = rep(rep(c("Astrea curta","Pocillopora verrucosa"), each = 50),3))  %>% 
  mutate(log_area = log(Surface_Area)) %>% group_split(Species)
## ADULTS ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Adu")
data_ACR = data_ACR %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Adu")
data_POC = data_POC %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_adu = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_adu = data.frame(expand.grid(n_col, n_col), comb_adu) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_adu)

## BALANCED ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Bal")
data_ACR = data_ACR %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Bal")
data_POC = data_POC %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_bal = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_bal = data.frame(expand.grid(n_col, n_col), comb_bal) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_bal)

## JUVENILES ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Juv")
data_ACR = data_ACR %>% mutate(Calcif =((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS =((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Juv")
data_POC = data_POC %>% mutate(Calcif =((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS =((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_juv = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_juv = data.frame(expand.grid(n_col, n_col), comb_juv) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_juv)
# Dataset Acropora - Pocillopora
Comb_dataset_AST = rbind(comb_juv, comb_bal, comb_adu) %>% as.data.frame() %>% 
  mutate(Sampling = rep(c("1) Juveniles", "2) Balance", "3) Adults"), each = 2500), 
         ID = rep(seq(1,2500,1),3), Pocillopora = Pocillopora *2, Acropora = Acropora *2,
         Sum = Acropora + Pocillopora) %>% filter(Sum <= 100)

###NAPOPORA IRREGULARIS ###
# Combination between Acropora and Pocillopora
data_4 = data_3 %>% filter(Species %in% c("Napopora irregularis","Pocillopora verrucosa")) 
# Build a new dataset according to the random sampling with the 3 skewness conditions
data = data.frame(Dataset = c(rep("Juv", 100), rep("Bal", 100), rep("Adu", 100)),
                  Surface_Area = c(sort(ACR_Juv$sample$x),sort(POC_Juv$sample$x),sort(ACR_Bal$sample$x),
                                   sort(POC_Bal$sample$x),sort(ACR_Adu$sample$x),sort(POC_Adu$sample$x)),
                  Species = rep(rep(c("Napopora irregularis","Pocillopora verrucosa"), each = 50),3))  %>% 
  mutate(log_area = log(Surface_Area)) %>% group_split(Species)
## ADULTS ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Adu")
data_ACR = data_ACR %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Adu")
data_POC = data_POC %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_adu = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_adu = data.frame(expand.grid(n_col, n_col), comb_adu) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_adu)

## BALANCED ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Bal")
data_ACR = data_ACR %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Bal")
data_POC = data_POC %>% mutate(Calcif = ((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS=((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_bal = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_bal = data.frame(expand.grid(n_col, n_col), comb_bal) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_bal)

## JUVENILES ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Juv")
data_ACR = data_ACR %>% mutate(Calcif =((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS =((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Juv")
data_POC = data_POC %>% mutate(Calcif =((exp(Calcif.Estimate)/365*1000)/100.087)*1000, PhotoS =((exp(PhotoS.Estimate)*24/1000)/14.003)*1000) %>% 
  select(Dataset, Surface_Area, Species, Calcif, PhotoS) 
col_m_Calcif = cumsum(data_POC$Calcif) ; col_m_PhotoS = cumsum(data_POC$PhotoS)
# Define a matrix for Calcification
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_Calcif + row_m_Calcif[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_c = m
# Define a matrix for Photosynthesis
m = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m[i,] = col_m_PhotoS + row_m_PhotoS[i]} ; n_col = seq(1,50,1)
m_number = matrix(ncol = 50, nrow = 50) ; for (i in 1:50) {m_number[i,] = i + seq(1,50,1)} ; m_heatmap_p = m
# Define a matrix for Photosynthesis/Calcification
m_heatmap = m_heatmap_p/m_heatmap_c
# Aggregate the matrix hotosynthesis/Calcification
comb_juv = c(m_heatmap[1,], m_heatmap[2,], m_heatmap[3,], m_heatmap[4,], m_heatmap[5,], 
             m_heatmap[6,], m_heatmap[7,], m_heatmap[8,], m_heatmap[9,], m_heatmap[10,], 
             m_heatmap[11,], m_heatmap[12,], m_heatmap[13,], m_heatmap[14,], m_heatmap[15,], 
             m_heatmap[16,], m_heatmap[17,], m_heatmap[18,], m_heatmap[19,], m_heatmap[20,], 
             m_heatmap[21,], m_heatmap[22,], m_heatmap[23,], m_heatmap[24,], m_heatmap[25,], 
             m_heatmap[26,], m_heatmap[27,], m_heatmap[28,], m_heatmap[29,], m_heatmap[30,], 
             m_heatmap[31,], m_heatmap[32,], m_heatmap[33,], m_heatmap[34,], m_heatmap[35,], 
             m_heatmap[36,], m_heatmap[37,], m_heatmap[38,], m_heatmap[39,], m_heatmap[40,], 
             m_heatmap[41,], m_heatmap[42,], m_heatmap[43,], m_heatmap[44,], m_heatmap[45,], 
             m_heatmap[46,], m_heatmap[47,], m_heatmap[48,], m_heatmap[49,], m_heatmap[50,])
comb_juv = data.frame(expand.grid(n_col, n_col), comb_juv) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_juv)
# Dataset Acropora - Pocillopora
Comb_dataset_NAP = rbind(comb_juv, comb_bal, comb_adu) %>% as.data.frame() %>% 
  mutate(Sampling = rep(c("1) Juveniles", "2) Balance", "3) Adults"), each = 2500), 
         ID = rep(seq(1,2500,1),3), Pocillopora = Pocillopora *2, Acropora = Acropora *2,
         Sum = Acropora + Pocillopora) %>% filter(Sum <= 100)

#### FINAL DATASET ###
Comb_dataset = rbind(Comb_dataset_ACR, Comb_dataset_AST, Comb_dataset_MON, Comb_dataset_NAP, Comb_dataset_POR) %>% as.data.frame() %>% 
  mutate(Combination = rep(c("Acropora - Pocillopora", "Astrea - Pocillopora", "Montipora  Pocillopora", 
                             "Napopora  Pocillopora", "Porites - Pocillopora"), each = 3675)) %>% 
  ggplot(., aes(x=Pocillopora, y=Acropora, col = Ratio, fill=Ratio)) + geom_tile() + scale_fill_viridis_c(option="inferno") + 
  scale_color_viridis_c(option="inferno") + theme_classic() + xlab("Pocillopora cover (%)") + ylab("X cover (%)") + facet_grid(Sampling~Combination)




###################################
#### Relationship NPP - Calcif ####
###################################

mutate(Calcification = ((exp(log_calcif)/365*1000)/100.087)*1000,
       Net_Photosynthesis = ((exp(log_photo)*24/1000)/14.003)*1000) %>% 
  mutate(Calcification_predicted = (((exp(predict(mixt_calcif, Data_Raw_Metabo))/365*1000)/100.087)*1000)[,1],
         Net_Photosynthesis_predicted = (((exp(predict(mixt_photo, Data_Raw_Metabo))*24/1000)/14.003)*1000)[,1]) %>% 
  select(Species, Concerning_Box, Surface_area_cm2, Calcification, Calcification_predicted, Net_Photosynthesis, Net_Photosynthesis_predicted)
data_Fig_Diego = data.frame(Species = rep(data_Fig_Diego$Species,2),
                            Concerning_Box = rep(data_Fig_Diego$Concerning_Box,2),
                            Surface_Area = rep(data_Fig_Diego$Surface_area_cm2,2),
                            Calcification = c(data_Fig_Diego$Calcification, data_Fig_Diego$Calcification_predicted),
                            Net_Photosynthesis = c(data_Fig_Diego$Net_Photosynthesis, data_Fig_Diego$Net_Photosynthesis_predicted),
                            Value = c(rep("Observed", 250), rep("Predicted", 250)))
data_Fig_Diego = data_Fig_Diego %>%  mutate(Shape = paste(Concerning_Box, Value, sep ="_"))

mixt_Diego = brm(bf(Calcification ~ a*Net_Photosynthesis^b, 
                    a ~ 1 + (1|Species), b ~ 1 + (1|Species), nl = TRUE), iter = 5000,
                 data = data_Fig_Diego, family = gaussian(),
                 prior = c(prior(normal(10,10), nlpar = "a"), prior(normal(0.5,.5), nlpar = "b")),
                 control = list(adapt_delta = 0.999, max_treedepth = 30), chains = 3, cores = 3)

mixt_Diego_df = cbind(fitted(mixt_Diego),mixt_Diego$data) %>% 
  mutate(Shape = data_Fig_Diego$Shape)

ggplot(mixt_Diego_df, aes(x = Net_Photosynthesis, y = Calcification, col = Species)) + 
  geom_ribbon(aes(x = Net_Photosynthesis, ymin = Q2.5, ymax = Q97.5, fill = Species), alpha = .5, show.legend = F) + 
  geom_point(aes(shape = Shape), alpha = .75, show.legend = F, size = 2) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() +
  scale_shape_manual(values = c(15,0,16,1,17,2)) + 
  scale_y_continuous(name = expression("Calcification")) + 
  scale_x_continuous(name = expression("Photosynthesis")) 