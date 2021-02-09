# SET UP THE SCRIPT
rm(list=ls())
Results_directory = paste(getwd(),"/Results",sep="") ; Data_directory = paste(getwd(),"/Data",sep="")
setwd(Data_directory) ; options(mc.cores = parallel::detectCores())
library('tidyverse') ; library('parallel') ; library('brms') ; library('abind') ; library('readxl')
library('patchwork') ; library("tidybayes") ; library('ggpubr') ; library('rstatix')
Data_Raw_Metabo <- read_excel("Data_Metabo.xlsx")

if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() == "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")}

#############################################################################################################
#############################################################################################################
######################################## CORAL METABOLISM PRODUCTION ######################################## 
#############################################################################################################
#############################################################################################################

# BAYESIAN FRAMEWORK
## Managing the dataset
Data_Raw_Metabo = Data_Raw_Metabo %>% mutate(Calcif_kg_yr = Corrected_Production/100,
                                             log_area = log(Surface_area_cm2),
                                             log_calcif = log(Calcif_kg_yr),
                                             log_respi = log(NightValue),
                                             log_photo = log(NightValue-Dayvalue)) %>% 
  mutate(Species = fct_recode(Species, "Acropora hyacinthus" = "Acropora_hyacinthus",
  "Montipora verilli" = "Montipora_verilli", "Napopora irregularis" = "Napopora_irregularis",
  "Astrea curta" = "Astrea_curta", "Pocillopora verrucosa" = "Pocillopora_verrucosa",
  "Porites lutea" = "Porites_lutea")) 

#############################################################################################################
##############################  Launching the model for the calcification rate ############################## 
#############################################################################################################

# Model
mixt_calcif = brm(bf(log_calcif ~ a*log_area+b,
                     a ~ 1 + (1|Species), 
                     b ~ 1 + (1|Species), 
                     nl = TRUE), iter = 5000,
                  data = Data_Raw_Metabo, family = gaussian(),
                  prior = c(prior(normal(-8,2), nlpar = "b"), prior(normal(0.8,0.2), nlpar = "a")),
                  control = list(adapt_delta = 0.999, max_treedepth = 30), chains = 3, cores = 1)
CaCO3_df = cbind(fitted(mixt_calcif),mixt_calcif$data) 
## Plot the results
Fig_1A = ggplot(CaCO3_df, aes(x = log_area, y = log_calcif, col = Species)) + 
  geom_ribbon(aes(x = log_area, ymin = Q2.5, ymax = Q97.5, fill = Species), alpha = .5, show.legend = F) + 
  theme_classic() + geom_point(alpha = .5, show.legend = F) + scale_color_viridis_d() + scale_fill_viridis_d() + 
  facet_wrap(~Species, ncol = 6, labeller = labeller(Species = c("Acropora hyacinthus" = "A. hyacinthus", 
  "Montipora verilli" = "M. verilli", "Napopora irregularis" = "N. irregularis", "Astrea curta" = "A. curta",
  "Pocillopora verrucosa" = "P. verrucosa",  "Porites lutea" = "P. lutea"))) +
  scale_y_continuous(name = expression("log(Calcification rate) (kg.yr"^-1*")")) + 
  scale_x_continuous(name = expression("log(Surface Area) (cm"^2*")")) +
  theme(legend.text = element_text(face = "italic")) 
Fig_2A = mixt_calcif %>% spread_draws(c(b_a_Intercept,sd_Species__a_Intercept), r_Species__a[Species,]) %>%
  mutate(condition_mean = b_a_Intercept + r_Species__a,
         condition_sd = sd_Species__a_Intercept + r_Species__a) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = mean(condition_sd)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
  "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
  "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 1, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the calcification allometric model"))

#############################################################################################################
###############################  Launching the model for the respiration rate ############################### 
#############################################################################################################

# Model
mixt_respi = brm(bf(log_respi ~ a*log_area+b,
                    a ~ 1 + (1|Species), 
                    b ~ 1 + (1|Species), 
                    nl = TRUE), iter = 5000,
                 data = Data_Raw_Metabo, family = gaussian(),
                 prior = c(prior(normal(-4,3), nlpar = "b"), prior(normal(1,0.3), nlpar = "a")),
                 control = list(adapt_delta = 0.999, max_treedepth = 30), chains = 3, cores = 1)
respi_df = cbind(fitted(mixt_respi),mixt_respi$data) 
## Plot the results
Fig_1B = ggplot(respi_df, aes(x = log_area, y = log_respi, col = Species)) + 
  geom_ribbon(aes(x = log_area, ymin = Q2.5, ymax = Q97.5, fill = Species), alpha = .5) + 
  theme_classic() + geom_point(alpha = .5) + scale_color_viridis_d() + scale_fill_viridis_d() + 
  facet_wrap(~Species, ncol = 6, labeller = labeller(Species = c("Acropora hyacinthus" = "A. hyacinthus", 
  "Montipora verilli" = "M. verilli", "Napopora irregularis" = "N. irregularis", "Astrea curta" = "A. curta",
  "Pocillopora verrucosa" = "P. verrucosa",  "Porites lutea" = "P. lutea"))) +
  scale_y_continuous(name = expression("log(Respiration rate) (mg.h"^-1*")")) + 
  scale_x_continuous(name = expression("log(Surface Area) (cm"^2*")")) +
  theme(legend.text = element_text(face = "italic")) 
Fig_2B = mixt_respi %>% spread_draws(c(b_a_Intercept,sd_Species__a_Intercept), r_Species__a[Species,]) %>%
  mutate(condition_mean = b_a_Intercept + r_Species__a,
         condition_sd = sd_Species__a_Intercept + r_Species__a) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = mean(abs(condition_sd))) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
  "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
  "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 1, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the respiration allometric model"))

#############################################################################################################
##############################  Launching the model for the photosynthesis rate ############################# 
#############################################################################################################

# Model
mixt_photo = brm(bf(log_photo ~ a*log_area+b,
                    a ~ 1 + (1|Species), 
                    b ~ 1 + (1|Species), 
                    nl = TRUE), iter = 5000,
                 data = Data_Raw_Metabo, family = gaussian(),
                 prior = c(prior(normal(-4,3), nlpar = "b"), prior(normal(1,0.3), nlpar = "a")),
                 control = list(adapt_delta = 0.999, max_treedepth = 30), chains = 3, cores = 1)
photo_df = cbind(fitted(mixt_photo),mixt_photo$data) 
## Plot the results
Fig_1C = ggplot(photo_df, aes(x = log_area, y = (log_photo), col = Species)) + 
  geom_ribbon(aes(x = (log_area), ymin = (Q2.5), ymax = (Q97.5), fill = Species), alpha = .5, show.legend = F) + 
  theme_classic() + geom_point(alpha = .5, show.legend = F) + scale_color_viridis_d() + scale_fill_viridis_d() + 
  facet_wrap(~Species, ncol = 6, labeller = labeller(Species = c("Acropora hyacinthus" = "A. hyacinthus", 
  "Montipora verilli" = "M. verilli", "Napopora irregularis" = "N. irregularis", "Astrea curta" = "A. curta",
  "Pocillopora verrucosa" = "P. verrucosa",  "Porites lutea" = "P. lutea"))) +
  scale_y_continuous(name = expression("log(Photosynthesis rate) (mg.h"^-1*")")) + 
  scale_x_continuous(name = expression("log(Surface Area) (cm"^2*")")) +
  theme(legend.text = element_text(face = "italic")) 
Fig_2C = mixt_photo %>% spread_draws(c(b_a_Intercept,sd_Species__a_Intercept), r_Species__a[Species,]) %>%
  mutate(condition_mean = b_a_Intercept + r_Species__a,
         condition_sd = sd_Species__a_Intercept + r_Species__a) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = mean(abs(condition_sd))) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
  "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
  "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 1, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the Photosynthesis allometric model"))

#############################################################################################################
##########################################  Final Production Results ######################################## 
#############################################################################################################

Production = ((Fig_1A/Fig_1B/Fig_1C) | (Fig_2A/Fig_2B/Fig_2C)) + plot_layout(guides = "collect", widths = c(3, 2))
ggsave(Production, filename = paste(Results_directory,"Figure_1.eps", sep = "/"), device=cairo_ps, 
       fallback_resolution = 600, width = 40, height = 25, units = "cm")

#############################################################################################################
#############################################################################################################
####################################### CORAL METABOLISM PRODUCTIVITY ####################################### 
#############################################################################################################
#############################################################################################################

# BAYESIAN FRAMEWORK
## Managing the dataset
Data_Raw_Metabo = Data_Raw_Metabo %>% mutate(productivity_calcif = ((Calcif_kg_yr*1000)/Surface_area_cm2)/10,
                                             productivity_respi  = (NightValue)/Surface_area_cm2,
                                             productivity_photo  = (NightValue-Dayvalue)/Surface_area_cm2)

#############################################################################################################
##############################  Launching the model for the calcification rate ############################## 
#############################################################################################################

# Model
load("Models/mixt_calcif_PB.RData")
mixt_calcif_PB = brm(bf(productivity_calcif ~ a*Surface_area_cm2^b,
                     a ~ 1 + (1|Species), 
                     b ~ 1 + (1|Species), 
                     nl = TRUE), iter = 5000,
                  data = Data_Raw_Metabo, family = gaussian(),
                  prior = c(prior(normal(.25,.25), nlpar = "a"), prior(normal(-0.2,0.1), nlpar = "b")),
                  control = list(adapt_delta = 0.999, max_treedepth = 30), chains = 3, cores = 1)
CaCO3_PB_df = cbind(fitted(mixt_calcif_PB),mixt_calcif_PB$data) 
## Plot the results
Fig_3A = ggplot(CaCO3_PB_df, aes(x = log(Surface_area_cm2), y = log(productivity_calcif), col = Species)) + 
  geom_ribbon(aes(x = log(Surface_area_cm2), ymin = log(Q2.5), ymax = log(Q97.5), fill = Species), alpha = .5, show.legend = F) + 
  theme_classic() + geom_point(alpha = .5, show.legend = F) + scale_color_viridis_d() + scale_fill_viridis_d() + 
  facet_wrap(~Species, ncol = 6, labeller = labeller(Species = c("Acropora hyacinthus" = "A. hyacinthus", 
  "Montipora verilli" = "M. verilli", "Napopora irregularis" = "N. irregularis", "Astrea curta" = "A. curta",
  "Pocillopora verrucosa" = "P. verrucosa",  "Porites lutea" = "P. lutea"))) +
  scale_y_continuous(name = expression("log(Calcification rate std.) (kg.m"^-2*"yr"^-1*")")) + 
  scale_x_continuous(name = expression("log(Surface Area) (cm"^2*")")) +
  theme(legend.text = element_text(face = "italic")) 
Fig_4A = mixt_calcif_PB %>% spread_draws(c(b_b_Intercept,sd_Species__b_Intercept), r_Species__b[Species,]) %>%
  mutate(condition_mean = b_b_Intercept + r_Species__b,
         condition_sd = sd_Species__b_Intercept + r_Species__b) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = mean(condition_sd)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
  "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
  "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the calcification allometric model std."))   
Fig_4D = mixt_calcif_PB %>% spread_draws(c(b_a_Intercept,sd_Species__a_Intercept), r_Species__a[Species,]) %>%
  mutate(condition_mean = b_a_Intercept + r_Species__a,
         condition_sd = sd_Species__a_Intercept + r_Species__a) %>% 
  group_by(Species) %>% summarise(mean_mean = log(mean(condition_mean)), mean_sd = mean(condition_sd)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
  "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
  "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(alpha~"coefficient of the calcification allometric model std.")) 

#############################################################################################################
###############################  Launching the model for the respiration rate ############################### 
#############################################################################################################

# Model
load("Models/mixt_respi_PB.RData")
mixt_respi_PB = brm(bf(productivity_respi ~ a*Surface_area_cm2^b,
                        a ~ 1 + (1|Species), 
                        b ~ 1 + (1|Species), 
                        nl = TRUE), iter = 5000,
                     data = Data_Raw_Metabo, family = gaussian(),
                     prior = c(prior(normal(.2,.2), nlpar = "a"), prior(normal(0,.3), nlpar = "b")),
                     control = list(adapt_delta = 0.999, max_treedepth = 30), chains = 3, cores = 1)
respi_PB_df = cbind(fitted(mixt_respi_PB),mixt_respi_PB$data) 
## Plot the results
Fig_3B = ggplot(respi_PB_df, aes(x = log(Surface_area_cm2), y = log(productivity_respi), col = Species)) + 
  geom_ribbon(aes(x = log(Surface_area_cm2), ymin = log(Q2.5), ymax = log(Q97.5), fill = Species), alpha = .5) + 
  theme_classic() + geom_point(alpha = .5) + scale_color_viridis_d() + scale_fill_viridis_d() + 
  facet_wrap(~Species, ncol = 6, labeller = labeller(Species = c("Acropora hyacinthus" = "A. hyacinthus", 
  "Montipora verilli" = "M. verilli", "Napopora irregularis" = "N. irregularis", "Astrea curta" = "A. curta",
  "Pocillopora verrucosa" = "P. verrucosa",  "Porites lutea" = "P. lutea"))) +
  scale_y_continuous(name = expression("log(Respiration rate std.) (mg.cm"^-2*"h"^-1*")")) + 
  scale_x_continuous(name = expression("log(Surface Area) (cm"^2*")")) +
  theme(legend.text = element_text(face = "italic")) 
Fig_4B = mixt_respi_PB %>% spread_draws(c(b_b_Intercept,sd_Species__b_Intercept), r_Species__b[Species,]) %>%
  mutate(condition_mean = b_b_Intercept + r_Species__b,
         condition_sd = sd_Species__b_Intercept + r_Species__b) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = mean(condition_sd)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
  "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
  "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the respiration allometric model std."))                                             
Fig_4E = mixt_respi_PB %>% spread_draws(c(b_a_Intercept,sd_Species__a_Intercept), r_Species__a[Species,]) %>%
  mutate(condition_mean = b_a_Intercept + r_Species__a,
         condition_sd = sd_Species__a_Intercept + r_Species__a) %>% 
  group_by(Species) %>% summarise(mean_mean = log(mean(condition_mean)), mean_sd = mean(condition_sd)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
  "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
  "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(alpha~"coefficient of the respiration allometric model std.")) 

#############################################################################################################
##############################  Launching the model for the photosynthesis rate ############################# 
#############################################################################################################

# Model
load("Models/mixt_photo_PB.RData")
mixt_photo_PB = brm(bf(productivity_photo ~ a*Surface_area_cm2^b,
                       a ~ 1 + (1|Species), 
                       b ~ 1 + (1|Species), 
                       nl = TRUE), iter = 5000,
                    data = Data_Raw_Metabo, family = gaussian(),
                    prior = c(prior(normal(.2,.2), nlpar = "a"), prior(normal(0,.3), nlpar = "b")),
                    control = list(adapt_delta = 0.999, max_treedepth = 30), chains = 3, cores = 1)
photo_PB_df = cbind(fitted(mixt_photo_PB),mixt_photo_PB$data) 
## Plot the results
Fig_3C = ggplot(photo_PB_df, aes(x = log(Surface_area_cm2), y = log(productivity_photo), col = Species)) + 
  geom_ribbon(aes(x = log(Surface_area_cm2), ymin = log(Q2.5), ymax = log(Q97.5), fill = Species), alpha = .5, show.legend = F) + 
  theme_classic() + geom_point(alpha = .5, show.legend = F) + scale_color_viridis_d() + scale_fill_viridis_d() + 
  facet_wrap(~Species, ncol = 6, labeller = labeller(Species = c("Acropora hyacinthus" = "A. hyacinthus", 
  "Montipora verilli" = "M. verilli", "Napopora irregularis" = "N. irregularis", "Astrea curta" = "A. curta",
  "Pocillopora verrucosa" = "P. verrucosa",  "Porites lutea" = "P. lutea"))) +
  scale_y_continuous(name = expression("log(Photosynthesis rate std.) (mg.cm"^-2*"h"^-1*")")) + 
  scale_x_continuous(name = expression("log(Surface Area) (cm"^2*")")) +
  theme(legend.text = element_text(face = "italic")) 
Fig_4C = mixt_photo_PB %>% spread_draws(c(b_b_Intercept,sd_Species__b_Intercept), r_Species__b[Species,]) %>%
  mutate(condition_mean = b_b_Intercept + r_Species__b,
         condition_sd = sd_Species__b_Intercept + r_Species__b) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = mean(condition_sd)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
  "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
  "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the photosynthesis allometric model std.")) 
Fig_4F = mixt_photo_PB %>% spread_draws(c(b_a_Intercept,sd_Species__a_Intercept), r_Species__a[Species,]) %>%
  mutate(condition_mean = b_a_Intercept + r_Species__a,
         condition_sd = sd_Species__a_Intercept + r_Species__a) %>% 
  group_by(Species) %>% summarise(mean_mean = log(mean(condition_mean)), mean_sd = mean(condition_sd)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
  "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
  "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(alpha~"coefficient of the photosynthesis allometric model std.")) 

#############################################################################################################
##########################################  Final Production Results ######################################## 
#############################################################################################################

Productivity = ((Fig_3A/Fig_3B/Fig_3C) | (Fig_4A/Fig_4B/Fig_4C) | (Fig_4D/Fig_4E/Fig_4F)) + plot_layout(guides = "collect", widths = c(5, 1, 1))
ggsave(Productivity, filename = paste(Results_directory,"Figure_2.eps", sep = "/"), device=cairo_ps, 
       fallback_resolution = 600, width = 60, height = 25, units = "cm")

#############################################################################################################
#############################################################################################################
######################################## IMPORTANCE OF THE VARIABLES ######################################## 
#############################################################################################################
#############################################################################################################

# Calcification
Modele_selection = step(glm(formula = log_calcif ~ log_area*Species*BEF_Darling*Morpho, data = Data_Raw_Metabo))
Modele_selection$formula ; step(glm(log_calcif ~ log_area + Species + log_area:Species, data = Data_Raw_Metabo))

data_dev = data.frame(Deviation = 1 - abs(c((024.45-342.30)/342.30, (041.06-342.30)/342.30, (210.25-342.30)/342.30)), 
                      Variables = c("C","B","A"))
data_dev = data_dev %>% mutate(Deviation = Deviation/sum(Deviation)*100) %>% 
  mutate(Variables = fct_recode(Variables, "log(Area)" = "A", "Species" = "B", "log(Area):Species" = "C"))

Fig_5A = ggplot(data_dev, aes(y = Deviation, x = "", fill = Variables)) + 
  geom_bar(stat = "identity", col = "black", alpha = 0.7, show.legend = F) +
  coord_polar("y", start = 0) + theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values=c("#C70039","#FF5733","goldenrod1")) +
  scale_y_continuous(breaks= cumsum(data_dev$Deviation) - data_dev$Deviation/2, 
                     labels= paste(round(data_dev$Deviation),"%",sep=""))

# Respiration
Modele_selection = step(glm(formula = log_respi ~ log_area*Species*BEF_Darling*Morpho, data = Data_Raw_Metabo))
Modele_selection$formula ; step(glm(log_respi ~ log_area + Species, data = Data_Raw_Metabo))

data_dev = data.frame(Deviation = 1 - abs(c((101.42-426.5)/426.5, (178.38-426.5)/426.5, (397.06-426.5)/426.5)), 
                      Variables = c("C","B","A"))
data_dev = data_dev %>% mutate(Deviation = Deviation/sum(Deviation)*100) %>% 
  mutate(Variables = fct_recode(Variables, "log(Area)" = "A", "Species" = "B", "log(Area):Species" = "C"))

Fig_5B = ggplot(data_dev, aes(y = Deviation, x = "", fill = Variables)) + 
  geom_bar(stat = "identity", col = "black", alpha = 0.7) +
  coord_polar("y", start = 0) + theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values=c("#C70039","#FF5733","goldenrod1")) +
  scale_y_continuous(breaks= cumsum(data_dev$Deviation) - data_dev$Deviation/2, 
                     labels= paste(round(data_dev$Deviation),"%",sep=""))

# Photosynthesis
Modele_selection = step(glm(formula = log_photo ~ log_area*Species*BEF_Darling*Morpho, data = Data_Raw_Metabo))
Modele_selection$formula ; step(glm(log_photo ~ log_area + Species, data = Data_Raw_Metabo))

data_dev = data.frame(Deviation = 1 - abs(c((96.78-402.6)/402.6, (151.57-402.6)/402.6, (369.44-402.6)/402.6)), 
                      Variables = c("C","B","A"))
data_dev = data_dev %>% mutate(Deviation = Deviation/sum(Deviation)*100) %>% 
  mutate(Variables = fct_recode(Variables, "log(Area)" = "A", "Species" = "B", "log(Area):Species" = "C"))

Fig_5C = ggplot(data_dev, aes(y = Deviation, x = "", fill = Variables)) + 
  geom_bar(stat = "identity", col = "black", alpha = 0.7, show.legend = F) +
  coord_polar("y", start = 0) + theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values=c("#C70039","#FF5733","goldenrod1")) +
  scale_y_continuous(breaks= cumsum(data_dev$Deviation) - data_dev$Deviation/2, 
                     labels= paste(round(data_dev$Deviation),"%",sep=""))

Importance = ((Fig_5A/Fig_5B/Fig_5C)) + plot_layout(guides = "collect", ncol = 3)
ggsave(Importance, filename = paste(Results_directory,"Figure_3.eps", sep = "/"), device=cairo_ps, 
       fallback_resolution = 600, width = 10, height = 4, units = "cm")

#############################################################################################################
#############################################################################################################
################################################## EMMEANS ################################################## 
#############################################################################################################
#############################################################################################################

# Calcification
AOV_Calcif = Data_Raw_Metabo %>% 
  dplyr::select(productivity_calcif, Surface_area_cm2, Species) %>% 
  anova_test(productivity_calcif ~ Surface_area_cm2 + Species)
Pairwise_Calcif = Data_Raw_Metabo %>% 
  dplyr::select(productivity_calcif, Surface_area_cm2, Species) %>% 
  rstatix::emmeans_test(productivity_calcif ~ Species, covariate = Surface_area_cm2) %>% 
  add_xy_position(x = "group", fun = "mean_se")
Pairwise_Calcif$y.position = Pairwise_Calcif$y.position + 0.025
Pairwise_Calcif_plot = ggline(get_emmeans(Pairwise_Calcif), x = "Species", 
                              y = "emmean",  linetype = "dashed", color = "Species", legend = 'none') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, color = Species), width = 0.2) + 
  stat_pvalue_manual(Pairwise_Calcif, hide.ns = T, tip.length = F) + scale_color_viridis_d() +
  labs(subtitle = get_test_label(AOV_Calcif, detailed = TRUE))

# Respiration
AOV_Respi = Data_Raw_Metabo %>% 
  dplyr::select(productivity_respi, Surface_area_cm2, Species) %>% 
  anova_test(productivity_respi ~ Surface_area_cm2 + Species)
Pairwise_Respi = Data_Raw_Metabo %>% 
  dplyr::select(productivity_respi, Surface_area_cm2, Species) %>% 
  rstatix::emmeans_test(productivity_respi ~ Species, covariate = Surface_area_cm2) %>% 
  add_xy_position(x = "group", fun = "mean_se")
Pairwise_Respi$y.position = Pairwise_Respi$y.position + 0.025
Pairwise_Respi_plot = ggline(get_emmeans(Pairwise_Respi), x = "Species", 
                             y = "emmean",  linetype = "dashed", color = "Species", legend = 'none') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, color = Species), width = 0.2) + 
  stat_pvalue_manual(Pairwise_Respi, hide.ns = T, tip.length = F) + scale_color_viridis_d() +
  labs(subtitle = get_test_label(AOV_Respi, detailed = TRUE))

# Photosynthesis
AOV_Photo = Data_Raw_Metabo %>% 
  dplyr::select(productivity_photo, Surface_area_cm2, Species) %>% 
  anova_test(productivity_photo ~ Surface_area_cm2 + Species)
Pairwise_Photo = Data_Raw_Metabo %>% 
  dplyr::select(productivity_photo, Surface_area_cm2, Species) %>% 
  rstatix::emmeans_test(productivity_photo ~ Species, covariate = Surface_area_cm2) %>% 
  add_xy_position(x = "group", fun = "mean_se")
Pairwise_Photo$y.position = Pairwise_Photo$y.position + 0.025
Pairwise_Photo_plot = ggline(get_emmeans(Pairwise_Photo), x = "Species", 
                             y = "emmean",  linetype = "dashed", color = "Species", legend = 'none') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, color = Species), width = 0.2) + 
  stat_pvalue_manual(Pairwise_Photo, hide.ns = T, tip.length = F) + scale_color_viridis_d() +
  labs(subtitle = get_test_label(AOV_Photo, detailed = TRUE),caption = get_pwc_label(Pairwise_Photo))

# Final_Plot Supplementary
Supplementary = Pairwise_Calcif_plot | Pairwise_Respi_plot | Pairwise_Photo_plot
ggsave(Supplementary, filename = paste(Results_directory,"Figure_S1.eps", sep = "/"), device=cairo_ps, 
       fallback_resolution = 600, width = 60, height = 20, units = "cm")

# CONFIDENCE INTERVAL
# 1)
CI_Calcif_PB = mixt_calcif_PB %>% spread_draws(c(b_b_Intercept,sd_Species__b_Intercept), r_Species__b[Species,]) %>%
  mutate(condition_mean = b_b_Intercept + r_Species__b,
         condition_sd = sd_Species__b_Intercept + r_Species__b) %>% 
  mutate(condition = condition_mean + condition_sd) %>% group_by(Species) %>% group_split(~Species)
results = vector(length = 6) ; for (i in 1:6) {
  results[i] = paste(round(length(CI_Calcif_PB[[i]]$condition[CI_Calcif_PB[[i]]$condition>0])/7500,2)*100,"%", sep = "")}
CI_Results_Calcif_PB = data.frame(Species = c("A. hyacinthus", "A. curta", "M. verilli", 
                                              "N. irregularis", "P. verrucosa", "P. lutea"),CI = results)
CI_Calcif = mixt_calcif %>% spread_draws(c(b_a_Intercept,sd_Species__a_Intercept), r_Species__a[Species,]) %>%
  mutate(condition_mean = b_a_Intercept + r_Species__a,
         condition_sd = sd_Species__a_Intercept + r_Species__a) %>% 
  mutate(condition = condition_mean) %>% group_by(Species) %>% group_split(~Species)
results = vector(length = 6) ; for (i in 1:6) {
  results[i] = paste(round(length(CI_Calcif[[i]]$condition[CI_Calcif[[i]]$condition>1])/7500,2)*100,"%", sep = "")}
CI_Results_Calcif = data.frame(Species = c("A. hyacinthus", "A. curta", "M. verilli", 
                                              "N. irregularis", "P. verrucosa", "P. lutea"),CI = results)

# 2)
CI_Respi_PB = mixt_respi_PB %>% spread_draws(c(b_b_Intercept,sd_Species__b_Intercept), r_Species__b[Species,]) %>%
  mutate(condition_mean = b_b_Intercept + r_Species__b,
         condition_sd = sd_Species__b_Intercept + r_Species__b) %>% 
  mutate(condition = condition_mean) %>% group_by(Species) %>% group_split(~Species)
results = vector(length = 6) ; for (i in 1:6) {
  results[i] = paste(round(length(CI_Respi_PB[[i]]$condition[CI_Respi_PB[[i]]$condition>0])/7500,2)*100,"%", sep = "")}
CI_Results_Respi_PB = data.frame(Species = c("A. hyacinthus", "A. curta", "M. verilli", 
                                              "N. irregularis", "P. verrucosa", "P. lutea"),CI = results)
CI_Respi = mixt_respi %>% spread_draws(c(b_a_Intercept,sd_Species__a_Intercept), r_Species__a[Species,]) %>%
  mutate(condition_mean = b_a_Intercept + r_Species__a,
         condition_sd = sd_Species__a_Intercept + r_Species__a) %>% 
  mutate(condition = condition_mean) %>% group_by(Species) %>% group_split(~Species)
results = vector(length = 6) ; for (i in 1:6) {
  results[i] = paste(round(length(CI_Respi[[i]]$condition[CI_Respi[[i]]$condition>1])/7500,2)*100,"%", sep = "")}
(CI_Results_Respi = data.frame(Species = c("A. hyacinthus", "A. curta", "M. verilli", 
                                           "N. irregularis", "P. verrucosa", "P. lutea"), CI = results)
)
# 3)
CI_Photo_PB = mixt_photo_PB %>% spread_draws(c(b_b_Intercept,sd_Species__b_Intercept), r_Species__b[Species,]) %>%
  mutate(condition_mean = b_b_Intercept + r_Species__b,
         condition_sd = sd_Species__b_Intercept + r_Species__b) %>% 
  mutate(condition = condition_mean) %>% group_by(Species) %>% group_split(~Species)
results = vector(length = 6) ; for (i in 1:6) {
  results[i] = paste(round(length(CI_Photo_PB[[i]]$condition[CI_Photo_PB[[i]]$condition<0.2])/7500,2)*100,"%", sep = "")}
CI_Results_Photo_PB = data.frame(Species = c("A. hyacinthus", "A. curta", "M. verilli", 
                                             "N. irregularis", "P. verrucosa", "P. lutea"), CI = results)
CI_Photo = mixt_photo %>% spread_draws(c(b_a_Intercept,sd_Species__a_Intercept), r_Species__a[Species,]) %>%
  mutate(condition_mean = b_a_Intercept + r_Species__a,
         condition_sd = sd_Species__a_Intercept + r_Species__a) %>% 
  mutate(condition = condition_mean) %>% group_by(Species) %>% group_split(~Species)
results = vector(length = 6) ; for (i in 1:6) {
  results[i] = paste(round(length(CI_Photo[[i]]$condition[CI_Photo[[i]]$condition<1.2])/7500,2)*100,"%", sep = "")}
(CI_Results_Photo = data.frame(Species = c("A. hyacinthus", "A. curta", "M. verilli", 
                                          "N. irregularis", "P. verrucosa", "P. lutea"), CI = results))

data_bayesien_Coef = cbind(CI_Results_Calcif, CI_Results_Calcif_PB[,2], CI_Results_Respi[,2], 
      CI_Results_Respi_PB[,2], CI_Results_Photo[,2], CI_Results_Photo_PB[,2])
colnames(data_bayesien_Coef) = c("Species","Calcif_prod","Calcif_P/B","Respi","Respi_P/B","Photo","Photo_P/B")
data_bayesien_Coef
