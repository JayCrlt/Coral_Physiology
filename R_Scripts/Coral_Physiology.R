# SET UP THE SCRIPT
rm(list=ls());Results_directory=paste(getwd(),"/Results",sep="");Data_directory=paste(getwd(),"/Data",sep="")
setwd(Data_directory) ; options(mc.cores = parallel::detectCores()) ; load("data_Carlot.RData")
library('tidyverse');library('brms');library('ggpubr');library('abind');library('readxl');library('ggmcmc')
library('patchwork');library("tidybayes");library('rstatix');library("drawsample");library('parallel') 

if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() == "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")}

#### FIGURE 1 ####
# Model Calcification

my_priors <- set_prior("normal(0, 5)", class = "b") +
  set_prior("normal(0, 5)", class = "Intercept") +
  set_prior("gamma(2, 0.1)", class = "sigma")
mixt_calcif <- brm(log_calcif ~ log_area + (1 + log_area | Species),
                   iter = 5000, data = Data_Raw_Metabo, family = gaussian(),
                   control = list(adapt_delta = 0.999, max_treedepth = 30),
                   prior = my_priors, chains = 3)
plot(mixt_calcif) # traceplot and posterior distributions, add to online supp
pp_check(mixt_calcif, type = "scatter_avg") # posterior_predictive check, add to online supp
bayes_R2(mixt_calcif) # bayesian R2, add to your main figure and results section?
CaCO3_df = cbind(fitted(mixt_calcif),mixt_calcif$data) 

## Plot the results
Fig_1A = ggplot(CaCO3_df, aes(x = log_area, y = log_calcif, col = Species)) + 
  geom_ribbon(aes(x = log_area, ymin = Q2.5, ymax = Q97.5, fill = Species), alpha = .5, show.legend = F) + 
  theme_classic() + geom_point(alpha = .5, show.legend = F) + scale_color_viridis_d()+scale_fill_viridis_d() + 
  facet_wrap(~Species, ncol = 6, labeller = labeller(Species = c("Acropora hyacinthus" = "A. hyacinthus", 
                                                                 "Montipora verilli" = "M. verilli", 
                                                                 "Napopora irregularis" = "N. irregularis", 
                                                                 "Astrea curta" = "A. curta",
                                                                 "Pocillopora verrucosa" = "P. verrucosa",  
                                                                 "Porites lutea" = "P. lutea"))) +
  scale_y_continuous(name = expression("log(Calcification rate) (kg.yr"^-1*")")) + 
  scale_x_continuous(name = expression("log(Surface Area) (cm"^2*")")) +
  theme(legend.text = element_text(face = "italic")) 
Fig_2A = mixt_calcif %>% spread_draws(c(b_log_area,sd_Species__log_area), r_Species[Species,log_area]) %>%
  as.data.frame() %>% filter(log_area == "log_area") %>% 
  mutate(condition_mean = b_log_area + r_Species,
         condition_sd = sd_Species__log_area + r_Species) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = sd(condition_mean)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
                              "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
                              "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", 
                              "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 1, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the calcification allometric model"))
Fig_3A = mixt_calcif %>% spread_draws(c(b_Intercept,sd_Species__Intercept), r_Species[Species,]) %>%
  mutate(condition_mean = b_Intercept + r_Species,
         condition_sd = sd_Species__Intercept + r_Species) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = sd(condition_mean)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
                              "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
                              "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", 
                              "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(alpha~"coefficient of the calcification allometric model"))

# Model Respiration
mixt_respi = brm(log_respi ~ log_area + (1 + log_area | Species),
    iter = 5000, data = Data_Raw_Metabo, family = gaussian(),
    control = list(adapt_delta = 0.999, max_treedepth = 30),
    prior = my_priors, chains = 3)
plot(mixt_respi) # traceplot and posterior distributions, add to online supp
pp_check(mixt_respi, type = "scatter_avg") # posterior_predictive check, add to online supp
bayes_R2(mixt_respi) # bayesian R2, add to your main figure and results section?
respi_df = cbind(fitted(mixt_respi),mixt_respi$data) 
## Plot the results
Fig_1B = ggplot(respi_df, aes(x = log_area, y = log_respi, col = Species)) + 
  geom_ribbon(aes(x = log_area, ymin = Q2.5, ymax = Q97.5, fill = Species), alpha = .5) + 
  theme_classic() + geom_point(alpha = .5) + scale_color_viridis_d() + scale_fill_viridis_d() + 
  facet_wrap(~Species, ncol = 6, labeller = labeller(Species = c("Acropora hyacinthus" = "A. hyacinthus", 
                                                                 "Montipora verilli" = "M. verilli", 
                                                                 "Napopora irregularis" = "N. irregularis", 
                                                                 "Astrea curta" = "A. curta",
                                                                 "Pocillopora verrucosa" = "P. verrucosa",  
                                                                 "Porites lutea" = "P. lutea"))) +
  scale_y_continuous(name = expression("log(Respiration rate) (mg.h"^-1*")")) + 
  scale_x_continuous(name = expression("log(Surface Area) (cm"^2*")")) +
  theme(legend.text = element_text(face = "italic")) 
Fig_2B = mixt_respi %>% spread_draws(c(b_log_area,sd_Species__log_area), r_Species[Species,log_area]) %>%
  as.data.frame() %>% filter(log_area == "log_area") %>% 
  mutate(condition_mean = b_log_area + r_Species,
         condition_sd = sd_Species__log_area + r_Species) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = sd(condition_mean)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
                              "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
                              "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", 
                              "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 1, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the respiration allometric model"))
Fig_3B = mixt_respi %>% spread_draws(c(b_Intercept,sd_Species__Intercept), r_Species[Species,]) %>%
  mutate(condition_mean = b_Intercept + r_Species,
         condition_sd = sd_Species__Intercept + r_Species) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = sd(condition_mean)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
                              "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
                              "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", 
                              "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(alpha~"coefficient of the respiration allometric model"))

# Model Photosynthesis
mixt_photo = brm(log_photo ~ log_area + (1 + log_area | Species),
                 iter = 5000, data = Data_Raw_Metabo, family = gaussian(),
                 control = list(adapt_delta = 0.999, max_treedepth = 35),
                 prior = my_priors, chains = 3)
plot(mixt_photo) # traceplot and posterior distributions, add to online supp
pp_check(mixt_photo, type = "scatter_avg") # posterior_predictive check, add to online supp
bayes_R2(mixt_photo) # bayesian R2, add to your main figure and results section?
photo_df = cbind(fitted(mixt_photo),mixt_photo$data) 
## Plot the results
Fig_1C = ggplot(photo_df, aes(x = log_area, y = (log_photo), col = Species)) + 
  geom_ribbon(aes(x = (log_area), ymin = (Q2.5), ymax = (Q97.5), fill = Species), alpha = .5, show.legend = F) + 
  theme_classic() + geom_point(alpha = .5, show.legend = F) + scale_color_viridis_d() + scale_fill_viridis_d() + 
  facet_wrap(~Species, ncol = 6, labeller = labeller(Species = c("Acropora hyacinthus" = "A. hyacinthus", 
                                                                 "Montipora verilli" = "M. verilli", 
                                                                 "Napopora irregularis" = "N. irregularis", 
                                                                 "Astrea curta" = "A. curta",
                                                                 "Pocillopora verrucosa" = "P. verrucosa",  
                                                                 "Porites lutea" = "P. lutea"))) +
  scale_y_continuous(name = expression("log(Photosynthesis rate) (mg.h"^-1*")")) + 
  scale_x_continuous(name = expression("log(Surface Area) (cm"^2*")")) +
  theme(legend.text = element_text(face = "italic")) 
Fig_2C = mixt_photo %>% spread_draws(c(b_log_area,sd_Species__log_area), r_Species[Species,log_area]) %>%
  as.data.frame() %>% filter(log_area == "log_area") %>% 
  mutate(condition_mean = b_log_area + r_Species,
         condition_sd = sd_Species__log_area + r_Species) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = sd(condition_mean)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
                              "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
                              "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", 
                              "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 1, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the Photosynthesis allometric model"))
Fig_3C = mixt_photo %>% spread_draws(c(b_Intercept,sd_Species__Intercept), r_Species[Species,]) %>%
  mutate(condition_mean = b_Intercept + r_Species,
         condition_sd = sd_Species__Intercept + r_Species) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = sd(condition_mean)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
                              "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
                              "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", 
                              "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(alpha~"coefficient of the photosynthesis allometric model"))

#Final Plot – Figure 1
Production = ((Fig_1A/Fig_1B/Fig_1C) | (Fig_2A/Fig_2B/Fig_2C)) + plot_layout(guides = "collect", widths = c(4, 1))
ggsave(Production, filename = paste(Results_directory,"Figure_1.eps", sep = "/"), device=cairo_ps, 
       fallback_resolution = 500, width = 40, height = 25, units = "cm")

#### FIGURE 2 ####

data_ACR = sample(rnorm(1000,mean=300,sd=200)) ; data_POC = data_ACR
Area = data.frame(ID = seq(1,1000,1), Surface_Area_ACR = data_ACR, Surface_Area_POC = data_POC) %>% 
  filter(Surface_Area_ACR > 0, Surface_Area_POC > 0)

ACR_Juv = draw_sample(dist=Area[,c(1,2)], n=50, skew = 0.8, kurts = 0, replacement = T ,output_name = c("sample", "1"))
ACR_Bal = draw_sample(dist=Area[,c(1,2)], n=50, skew = 0.5, kurts = 0, replacement = T ,output_name = c("sample", "1"))
ACR_Adu = draw_sample(dist=Area[,c(1,2)], n=50, skew = 0.2, kurts = 0, replacement = T ,output_name = c("sample", "1"))
POC_Juv = ACR_Juv ; POC_Bal = ACR_Bal ; POC_Adu = ACR_Adu

conv_mmol_CaCO3 = 1/365*1000/100.087*1000 #year_to_day ; kg_to_g ; g_CaCO3_to_mol ; mol_to_mmol
conv_mmol_C = 1*24/1000/14.003*1000 #hour_to_day ; mg_to_g ; g_C_to_mol ; mol_to_mmol

# Combination between Acropora and Pocillopora
#### THIS PROCESS SHOULD BE DONE FOR THE 4 OTHERS COMBINATIONS (i.e. Porites, Napopora, Astrea and Montipora)
# Build a new dataset according to the random sampling with the 3 skewness conditions
data = data.frame(Dataset = c(rep("Juv", 100), rep("Bal", 100), rep("Adu", 100)),
                  Surface_Area = c(sort(ACR_Juv$sample$x),sort(POC_Juv$sample$x),sort(ACR_Bal$sample$x),
                                   sort(POC_Bal$sample$x),sort(ACR_Adu$sample$x),sort(POC_Adu$sample$x)),
                  Species = rep(rep(c("Acropora hyacinthus","Pocillopora verrucosa"), each = 50),3))  %>% 
  mutate(log_area = log(Surface_Area)) %>% group_split(Species)
## ADULTS ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), PhotoS = 
                        predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Adu")
data_ACR = data_ACR %>% mutate(Calcif = exp(Calcif.Estimate)*conv_mmol_CaCO3, 
                               PhotoS = exp(PhotoS.Estimate)*conv_mmol_C) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), 
                      PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Adu")
data_POC = data_POC %>% mutate(Calcif = exp(Calcif.Estimate)*conv_mmol_CaCO3, 
                               PhotoS = exp(PhotoS.Estimate)*conv_mmol_C) %>% 
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
comb_adu = as.vector(t(as.matrix(m_heatmap)))
comb_adu = data.frame(expand.grid(n_col, n_col), comb_adu) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_adu)

## BALANCED ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), 
                      PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Bal")
data_ACR = data_ACR %>% mutate(Calcif = exp(Calcif.Estimate)*conv_mmol_CaCO3, 
                               PhotoS = exp(PhotoS.Estimate)*conv_mmol_C) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), 
                      PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Bal")
data_POC = data_POC %>% mutate(Calcif = exp(Calcif.Estimate)*conv_mmol_CaCO3, 
                               PhotoS = exp(PhotoS.Estimate)*conv_mmol_C) %>% 
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
comb_bal = as.vector(t(as.matrix(m_heatmap)))
comb_bal = data.frame(expand.grid(n_col, n_col), comb_bal) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_bal)

## JUVENILES ##
# Predict Acropora data to build a matrix
data_ACR = data.frame(data[[1]], Calcif = predict(mixt_calcif, data[[1]]), 
                      PhotoS = predict(mixt_photo, data[[1]])) %>% filter(Dataset == "Juv")
data_ACR = data_ACR %>% mutate(Calcif = exp(Calcif.Estimate)*conv_mmol_CaCO3, 
                               PhotoS = exp(PhotoS.Estimate)*conv_mmol_C) %>%
  select(Dataset, Surface_Area, Species, Calcif, PhotoS)
row_m_Calcif = cumsum(data_ACR$Calcif) ; row_m_PhotoS = cumsum(data_ACR$PhotoS)
# Predict Pocillopora data to build a matrix
data_POC = data.frame(data[[2]], Calcif = predict(mixt_calcif, data[[2]]), 
                      PhotoS = predict(mixt_photo, data[[2]])) %>% filter(Dataset == "Juv")
data_POC = data_POC %>% mutate(Calcif = exp(Calcif.Estimate)*conv_mmol_CaCO3, 
                               PhotoS = exp(PhotoS.Estimate)*conv_mmol_C) %>% 
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
t = as.vector(t(as.matrix(m_heatmap)))
comb_juv = as.vector(t(as.matrix(m_heatmap)))
comb_juv = data.frame(expand.grid(n_col, n_col), comb_juv) %>% rename(Pocillopora = Var1, Acropora = Var2, Ratio = comb_juv)
# Dataset Acropora - Pocillopora
Left_plot = rbind(comb_juv, comb_bal, comb_adu) %>% as.data.frame() %>% 
  mutate(Sampling = rep(c("1) Juveniles", "2) Balance", "3) Adults"), each = 2500), 
         ID = rep(seq(1,2500,1),3), Pocillopora = Pocillopora *2, Acropora = Acropora *2,
         Sum = Acropora + Pocillopora) %>% filter(Sum <= 100)  %>% 
  ggplot(., aes(x=Pocillopora, y=Acropora, col = Ratio, fill=Ratio)) + geom_tile() + scale_fill_viridis_c(option="inferno") + 
  scale_color_viridis_c(option="inferno") + theme_classic() + xlab("Pocillopora cover (%)") + ylab("Acropora cover (%)") + 
  facet_grid(~Sampling)

data_aggregated = data.frame(Species = rep(Data_Raw_Metabo$Species,2),
                            Concerning_Box = rep(Data_Raw_Metabo$Incubation_chamber,2),
                            Surface_Area = rep(Data_Raw_Metabo$Surface_Area_cm2,2),
                            Calcification = c(Data_Raw_Metabo$Calcification_mmol_day, 
                                              Data_Raw_Metabo$Calcification_mmol_day_estimate),
                            Net_Photosynthesis = c(Data_Raw_Metabo$Photosynthesis_mmol_day, 
                                                   Data_Raw_Metabo$Photosynthesis_mmol_day_estimate),
                            Value = c(rep("Observed", 250), rep("Predicted", 250))) %>% 
  mutate(Shape = paste(Concerning_Box, Value, sep ="_"))

mixt_Energy = brm(bf(Calcification ~ a*Net_Photosynthesis^b, 
                    a ~ 1 + (1|Species), b ~ 1 + (1|Species), nl = TRUE), iter = 5000,
                 data = data_aggregated, family = gaussian(),
                 prior = c(prior(uniform(-10,10), nlpar = "a"), prior(normal(0.5,.5), nlpar = "b")),
                 control = list(adapt_delta = 0.999, max_treedepth = 30), chains = 3, cores = 1)
mixt_Energy_df = cbind(fitted(mixt_Energy),mixt_Energy$data) %>% mutate(Shape = data_aggregated$Shape)

Right_plot = ggplot(mixt_Energy_df, aes(x = Net_Photosynthesis, y = Calcification, col = Species)) + 
  geom_ribbon(aes(x = Net_Photosynthesis, ymin = Q2.5, ymax = Q97.5, fill = Species), alpha = .5, show.legend = F) + 
  geom_point(aes(shape = Shape), alpha = .75, show.legend = F, size = 2) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() +   scale_shape_manual(values = c(15,0,16,1,17,2)) + 
  scale_y_continuous(name = expression("Calcification")) + scale_x_continuous(name = expression("Photosynthesis")) 

#### FIGURE S1 ####
# The figure S1 done on Keynote (macOS application)

#### FIGURE S2 ####

# Model Calcification
Fig_4A = ggplot(CaCO3_df, aes(x = log_area, y = log(exp(log_calcif)/exp(log_area)), col = Species)) + 
  geom_ribbon(aes(x = log_area, ymin = log(exp(Q2.5)/exp(log_area)), ymax = log(exp(Q97.5)/exp(log_area)), fill = Species), 
              alpha = .5, show.legend = F) + 
  theme_classic() + geom_point(alpha = .5, show.legend = F) + scale_color_viridis_d()+scale_fill_viridis_d() + 
  facet_wrap(~Species, ncol = 6, labeller = labeller(Species = c("Acropora hyacinthus" = "A. hyacinthus", 
                                                                 "Montipora verilli" = "M. verilli", 
                                                                 "Napopora irregularis" = "N. irregularis", 
                                                                 "Astrea curta" = "A. curta",
                                                                 "Pocillopora verrucosa" = "P. verrucosa",  
                                                                 "Porites lutea" = "P. lutea"))) +
  scale_y_continuous(name = expression("log(Calcification std. rate) (kg.cm"^-2*"yr"^-1*")")) + 
  scale_x_continuous(name = expression("log(Surface Area) (cm"^2*")")) +
  theme(legend.text = element_text(face = "italic")) 
Fig_5A = mixt_calcif %>% spread_draws(c(b_log_area,sd_Species__log_area), r_Species[Species,log_area]) %>%
  as.data.frame() %>% filter(log_area == "log_area") %>% 
  mutate(condition_mean = b_log_area + r_Species - 1,
         condition_sd = sd_Species__log_area + r_Species - 1) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = sd(condition_mean)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
                              "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
                              "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", 
                              "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the calcification allometric model"))
Fig_6A = data.frame(log_area = rep(0, 6), Species = unique(CaCO3_df$Species)) %>% cbind(predict(mixt_calcif, .)) %>% 
  as.data.frame() %>% mutate(mean_mean = Estimate, mean_sd = Est.Error) %>% 
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora hyacinthus",
                              "M. verilli" = "Montipora verilli", "N. irregularis" = "Napopora irregularis", 
                              "A. curta" = "Astrea curta", "P. verrucosa" = "Pocillopora verrucosa", 
                              "P. lutea" = "Porites lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(alpha~"coefficient of the calcification allometric model"))

# Model Respiration
Fig_4B = ggplot(respi_df, aes(x = log_area, y = log(exp(log_respi)/exp(log_area)), col = Species)) + 
  geom_ribbon(aes(x = log_area, ymin = log(exp(Q2.5)/exp(log_area)), ymax = log(exp(Q97.5)/exp(log_area)), fill = Species), 
              alpha = .5, show.legend = F) + 
  theme_classic() + geom_point(alpha = .5, show.legend = F) + scale_color_viridis_d()+scale_fill_viridis_d() + 
  facet_wrap(~Species, ncol = 6, labeller = labeller(Species = c("Acropora hyacinthus" = "A. hyacinthus", 
                                                                 "Montipora verilli" = "M. verilli", 
                                                                 "Napopora irregularis" = "N. irregularis", 
                                                                 "Astrea curta" = "A. curta",
                                                                 "Pocillopora verrucosa" = "P. verrucosa",  
                                                                 "Porites lutea" = "P. lutea"))) +
  scale_y_continuous(name = expression("log(Respiration rate std.) (mg.cm"^-2*"h"^-1*")")) + 
  scale_x_continuous(name = expression("log(Surface Area) (cm"^2*")")) +
  theme(legend.text = element_text(face = "italic")) 
Fig_5B = mixt_respi %>% spread_draws(c(b_log_area,sd_Species__log_area), r_Species[Species,log_area]) %>%
  as.data.frame() %>% filter(log_area == "log_area") %>% 
  mutate(condition_mean = b_log_area + r_Species - 1,
         condition_sd = sd_Species__log_area + r_Species - 1) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = sd(condition_mean)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
                              "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
                              "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", 
                              "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the respiration allometric model"))
Fig_6B = data.frame(log_area = rep(0, 6), Species = unique(respi_df$Species)) %>% cbind(predict(mixt_respi, .)) %>% 
  as.data.frame() %>% mutate(mean_mean = Estimate, mean_sd = Est.Error) %>% 
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora hyacinthus",
                              "M. verilli" = "Montipora verilli", "N. irregularis" = "Napopora irregularis", 
                              "A. curta" = "Astrea curta", "P. verrucosa" = "Pocillopora verrucosa", 
                              "P. lutea" = "Porites lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(alpha~"coefficient of the respiration allometric model"))

# Model Photosynthesis
Fig_4C = ggplot(photo_df, aes(x = log_area, y = log(exp(log_photo)/exp(log_area)), col = Species)) + 
  geom_ribbon(aes(x = log_area, ymin = log(exp(Q2.5)/exp(log_area)), ymax = log(exp(Q97.5)/exp(log_area)), fill = Species), 
              alpha = .5, show.legend = F) + 
  theme_classic() + geom_point(alpha = .5, show.legend = F) + scale_color_viridis_d()+scale_fill_viridis_d() + 
  facet_wrap(~Species, ncol = 6, labeller = labeller(Species = c("Acropora hyacinthus" = "A. hyacinthus", 
                                                                 "Montipora verilli" = "M. verilli", 
                                                                 "Napopora irregularis" = "N. irregularis", 
                                                                 "Astrea curta" = "A. curta",
                                                                 "Pocillopora verrucosa" = "P. verrucosa",  
                                                                 "Porites lutea" = "P. lutea"))) +
  scale_y_continuous(name = expression("log(Photosynthesis rate std.) (mg.cm"^-2*"h"^-1*")")) + 
  scale_x_continuous(name = expression("log(Surface Area) (cm"^2*")")) +
  theme(legend.text = element_text(face = "italic")) 
Fig_5C = mixt_photo %>% spread_draws(c(b_log_area,sd_Species__log_area), r_Species[Species,log_area]) %>%
  as.data.frame() %>% filter(log_area == "log_area") %>% 
  mutate(condition_mean = b_log_area + r_Species - 1,
         condition_sd = sd_Species__log_area + r_Species - 1) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = sd(condition_mean)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
                              "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
                              "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", 
                              "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the photosynthesis allometric model"))
Fig_6C = data.frame(log_area = rep(0, 6), Species = unique(photo_df$Species)) %>% cbind(predict(mixt_photo, .)) %>% 
  as.data.frame() %>% mutate(mean_mean = Estimate, mean_sd = Est.Error) %>% 
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora hyacinthus",
                              "M. verilli" = "Montipora verilli", "N. irregularis" = "Napopora irregularis", 
                              "A. curta" = "Astrea curta", "P. verrucosa" = "Pocillopora verrucosa", 
                              "P. lutea" = "Porites lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(alpha~"coefficient of the photosynthesis allometric model"))

#Final Plot – Figure S1
Productivity = ((Fig_4A/Fig_4B/Fig_4C) | (Fig_5A/Fig_5B/Fig_5C)) + plot_layout(guides = "collect", widths = c(4, 1))
ggsave(Productivity, filename = paste(Results_directory,"Figure_S1.eps", sep = "/"), device=cairo_ps, 
       fallback_resolution = 500, width = 40, height = 25, units = "cm")

#### FIGURE S3 ####

# Bayesian outputs – Calcification
modeltranformed_calcif <- ggmcmc::ggs(mixt_calcif) %>% filter(Parameter %in% unique(modeltranformed_calcif$Parameter)[1:6]) %>% filter(Iteration >= 2500)
C_calcif = ggplot(modeltranformed_calcif) + geom_line(aes(y = value, x = Iteration, color = as.factor(Chain))) + 
  facet_grid(Parameter ~ ., scale  = 'free_y', switch = 'y') + theme_classic() +
  labs(title = "Caterpillar Plots – Calcification", col   = "Chains") + fishualize::scale_colour_fish_d(option = "Chlorurus_microrhinos", end = .3)
D_calcif = ggplot(modeltranformed_calcif) + geom_density(aes(x = value), fill = "#6600CC60") + 
  facet_grid(Parameter ~ ., scale  = 'free', switch = 'y') + theme_classic() +
  labs(title = "Density Plots – Calcification") 

# Bayesian outputs – Respiration
modeltranformed_respi <- ggmcmc::ggs(mixt_respi) %>% filter(Parameter %in% unique(modeltranformed_calcif$Parameter)[1:6]) %>% filter(Iteration >= 2500)
C_respi = ggplot(modeltranformed_respi) + geom_line(aes(y = value, x = Iteration, color = as.factor(Chain))) + 
  facet_grid(Parameter ~ ., scale  = 'free_y', switch = 'y') + theme_classic() +
  labs(title = "Caterpillar Plots – Respiration", col   = "Chains") + fishualize::scale_colour_fish_d(option = "Chlorurus_microrhinos", end = .3)
D_respi = ggplot(modeltranformed_respi) + geom_density(aes(x = value), fill = "#6600CC60") + 
  facet_grid(Parameter ~ ., scale  = 'free', switch = 'y') + theme_classic() +
  labs(title = "Density Plots – Respiration") 

# Bayesian outputs – Photosynthesis
modeltranformed_photo <- ggmcmc::ggs(mixt_photo) %>% filter(Parameter %in% unique(modeltranformed_calcif$Parameter)[1:6]) %>% filter(Iteration >= 2500)
C_photo = ggplot(modeltranformed_photo) + geom_line(aes(y = value, x = Iteration, color = as.factor(Chain))) + 
  facet_grid(Parameter ~ ., scale  = 'free_y', switch = 'y') + theme_classic() +
  labs(title = "Caterpillar Plots – Photosynthesis", col   = "Chains") + fishualize::scale_colour_fish_d(option = "Chlorurus_microrhinos", end = .3)
D_photo = ggplot(modeltranformed_photo) + geom_density(aes(x = value), fill = "#6600CC60") + 
  facet_grid(Parameter ~ ., scale  = 'free', switch = 'y') + theme_classic() +
  labs(title = "Density Plots – Photosynthesis") 

#Final Plot – Figure S3
Bayes_Outputs_1 = C_calcif + D_calcif + C_respi + D_respi + C_photo + D_photo + plot_layout(guides = "collect", widths = c(3,1,3,1,3,1))
Bayes_Outputs_2 = pp_check(mixt_calcif, type = "scatter_avg") + pp_check(mixt_respi, type = "scatter_avg") + pp_check(mixt_photo, type = "scatter_avg")
Bayes_Outputs = Bayes_Outputs_1 / Bayes_Outputs_2

ggsave(Bayes_Outputs, filename = paste(Results_directory,"Figure_S3.eps", sep = "/"), device=cairo_ps, 
       fallback_resolution = 600, width = 40, height = 25, units = "cm")
