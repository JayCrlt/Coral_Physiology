# SET UP THE SCRIPT
rm(list=ls())
Results_directory = paste(getwd(),"/Results",sep="") ; Data_directory = paste(getwd(),"/Data",sep="")
setwd(Data_directory) ; options(mc.cores = parallel::detectCores())
library('tidyverse') ; library('parallel') ; library('brms') ; library('abind') ; library('readxl')
library('patchwork') ; library("tidybayes") ; Data_Raw_Metabo <- read_excel("Data_Metabo.xlsx")

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
load("data_Carlot.RData")
Data_Raw_Metabo = Data_Raw_Metabo %>% 
  mutate(productivity_calcif = Calcif_kg_yr*100/Surface_area_cm2,
         productivity_respi  = (NightValue)/Surface_area_cm2,
         productivity_photo  = (NightValue-Dayvalue)/Surface_area_cm2)

#############################################################################################################
##############################  Launching the model for the calcification rate ############################## 
#############################################################################################################

# Model
mixt_calcif_PB = brm(bf(productivity_calcif ~ a*Surface_area_cm2^b,
                     a ~ 1 + (1|Species), 
                     b ~ 1 + (1|Species), 
                     nl = TRUE), iter = 5000,
                  data = Data_Raw_Metabo, family = gaussian(),
                  prior = c(prior(normal(.25,.25), nlpar = "a"), prior(normal(-0.2,0.1), nlpar = "b")),
                  control = list(adapt_delta = 0.999, max_treedepth = 30), chains = 5)
CaCO3_PB_df = cbind(fitted(mixt_calcif_PB),mixt_calcif_PB$data) 
## Plot the results
Fig_3A = ggplot(CaCO3_PB_df, aes(x = log(Surface_area_cm2), y = log(productivity_calcif), col = Species)) + 
  geom_ribbon(aes(x = log(Surface_area_cm2), ymin = log(Q2.5), ymax = log(Q97.5), fill = Species), alpha = .5, show.legend = F) + 
  theme_classic() + geom_point(alpha = .5, show.legend = F) + scale_color_viridis_d() + scale_fill_viridis_d() + 
  facet_wrap(~Species, ncol = 6, labeller = labeller(Species = c("Acropora hyacinthus" = "A. hyacinthus", 
  "Montipora verilli" = "M. verilli", "Napopora irregularis" = "N. irregularis", "Astrea curta" = "A. curta",
  "Pocillopora verrucosa" = "P. verrucosa",  "Porites lutea" = "P. lutea"))) +
  scale_y_continuous(name = expression("log(Calcification std. rate) (kg.m"^-2*"yr"^-1*")")) + 
  scale_x_continuous(name = expression("log(Surface Area) (cm"^2*")")) +
  theme(legend.text = element_text(face = "italic")) 
Fig_4A = mixt_calcif_PB %>% spread_draws(c(b_b_Intercept,sd_Species__b_Intercept), r_Species__b[Species,]) %>%
  mutate(condition_mean = b_b_Intercept + r_Species__b,
         condition_sd = sd_Species__b_Intercept + r_Species__b) %>% 
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = sd(condition_mean)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
  "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
  "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the calcification std. allometric model"))          

#############################################################################################################
###############################  Launching the model for the respiration rate ############################### 
#############################################################################################################

# Model
mixt_respi_PB = brm(bf(productivity_respi ~ a*Surface_area_cm2^b,
                        a ~ 1 + (1|Species), 
                        b ~ 1 + (1|Species), 
                        nl = TRUE), iter = 5000,
                     data = Data_Raw_Metabo, family = gaussian(),
                     prior = c(prior(normal(.2,.2), nlpar = "a"), prior(normal(0,.3), nlpar = "b")),
                     control = list(adapt_delta = 0.999, max_treedepth = 30), chains = 5)
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
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = sd(condition_mean)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
  "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
  "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the respiration std. allometric model"))                                             


#############################################################################################################
##############################  Launching the model for the photosynthesis rate ############################# 
#############################################################################################################

# Model
mixt_photo_PB = brm(bf(productivity_photo ~ a*Surface_area_cm2^b,
                       a ~ 1 + (1|Species), 
                       b ~ 1 + (1|Species), 
                       nl = TRUE), iter = 5000,
                    data = Data_Raw_Metabo, family = gaussian(),
                    prior = c(prior(normal(.2,.2), nlpar = "a"), prior(normal(0,.3), nlpar = "b")),
                    control = list(adapt_delta = 0.999, max_treedepth = 30), chains = 5)
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
  group_by(Species) %>% summarise(mean_mean = mean(condition_mean), mean_sd = sd(condition_mean)) %>%
  mutate(Species = fct_recode(Species, "A. hyacinthus" = "Acropora.hyacinthus",
  "M. verilli" = "Montipora.verilli", "N. irregularis" = "Napopora.irregularis", 
  "A. curta" = "Astrea.curta", "P. verrucosa" = "Pocillopora.verrucosa", "P. lutea" = "Porites.lutea")) %>% 
  ggplot(aes(x = Species, dist = "norm", arg1 = mean_mean, arg2 = mean_sd, fill = Species)) + 
  stat_dist_eye(position = "dodge", aes(col = Species), alpha = .7, show.legend = F) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .5) +
  theme_classic() + scale_color_viridis_d() + scale_fill_viridis_d() + ylab("") + 
  scale_x_discrete(name = expression(beta~"coefficient of the photosynthesis std. allometric model"))                                             

#############################################################################################################
##########################################  Final Production Results ######################################## 
#############################################################################################################

Productivity = ((Fig_3A/Fig_3B/Fig_3C) | (Fig_4A/Fig_4B/Fig_4C)) + plot_layout(guides = "collect", widths = c(3, 2))
ggsave(Productivity, filename = paste(Results_directory,"Figure_S2_raw.eps", sep = "/"), device=cairo_ps, 
       fallback_resolution = 600, width = 40, height = 25, units = "cm")
