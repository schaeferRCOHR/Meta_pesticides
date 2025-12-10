# R code to reproduce the results of
# Schäfer et al. 2026
# How pesticide exposure and effects match with the intention of European pesticide regulation – a mini review

library(ggplot2)
library(scales)
library(dplyr)
library(forcats)
library(metafor)
library(betareg)
library(nnet)

# Set local path to working directory, where files and history is stored
setwd("Your/path")

# Load data sets
exp <- read.csv("https://raw.githubusercontent.com/schaeferRCOHR/Meta_pesticides/refs/heads/main/Exposure_overview.csv", sep = ";", dec = ".", na.string = "NA")

# remove studies based on exclusion criteria
exp2 <- exp %>%
		filter(is.na(Reason.for.exclusion))

summary(factor(exp2$Regulatory.threshold))
# do not consider time

summary(factor(exp2$Observation.unit))
# focus on pesticide residues across sites, sites and pesticides

exp3 <- exp2 %>%
			 filter(Regulatory.threshold != "DT", !Observation.unit %in% c(
      		 "Days where at least one pesticide exhibits exceedance",
      		 "Measured pesticide residues across all sites"))
			 
exp4 <- exp3 %>%
			 mutate(Observation.unit.new = recode(Observation.unit, 
			 "Detected pesticide residues across all sites" = "Concentrations", 
			 "Sites" = "Field sites"))

exp5 <- exp4 %>%
			 mutate(Observation.unit.new = factor(Observation.unit.new, 
			 levels = c("Field sites", "Pesticides", "Concentrations")),
			 Threshold_prop = Values_larger_threshold_rel/100,
			 Compartment_new = recode(Compartment, 
			 "Soil with cereals" = "Soil", 
			 "Soil with vegetables" = "Soil",
			 "Soil with grass" = "Soil",
			 "Soil with cereal-grass rotation" = "Soil",
			 "Soil of organic fields" = "Soil")
			 )
			 

head(exp5)

# test for relationship between proportions and observation unit
model1 <- betareg(Threshold_prop ~ Observation.unit.new + Compartment_new, data = exp5)
summary(model1)

# model checking
model1_res <- residuals(model1, type = "quantile")

summary(model1_res)
hist(model1_res, breaks = 20, main = "Histogram of randomized quantile residuals")
# looks roughly normal
qqnorm(model1_res)
# good fit
plot(fitted(model1), model1_res, xlab="Fitted values", ylab="Randomized quantile residuals")
# no systematic pattern

ggplot(exp5, aes(Observation.unit.new, Threshold_prop, fill = Observation.unit.new)) + 
      geom_violin(width = 0.9, alpha = 0.25, colour = NA, trim = FALSE) +
      geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6, colour = "grey20") +
      geom_point(aes(color = Observation.unit.new), position = position_jitter(width = 0.08, height = 0), 
      size = 2.4, alpha = 0.9) +
      scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1)) +
      scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
      scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
      labs(x = NULL, y = "Exceedance of threshold or model prediction") +
      theme_minimal(base_size = 18) + theme(legend.position = "none")
   
#################################   
# Analysis for biodiversity data  
#################################
    
# Preprocessing of data from Beaumelle et al 2023
# available at: https://github.com/leabeaumelle/MAPesti
beaum <- read.csv("https://raw.githubusercontent.com/leabeaumelle/MAPesti/refs/heads/main/Data/02PesticidesEffectsClean.csv", sep = ",", dec = ".", na.string = "NA")    

table(beaum$Country)
# remove studies not from Europe

beaum2 <- beaum %>%
			 filter(!Country %in% c(
      		 "Australia", "Brazil", "China", "Egypt", "Ghana", 
      		 "Israel", "Madagascar", "Mexico", "New Zealand", 
      		 "Philippines", "Sri Lanka", "USA"), !RecommendedRate %in% c(
      		 "not at recommended rate", "higher") 
      		 )

# reduce to authorised substances
auth_pest <- read.csv("https://raw.githubusercontent.com/schaeferRCOHR/Meta_pesticides/refs/heads/main/Authorisation_information.csv", sep = ";", dec = ".", na.string = "NA")
 
beaum2_filtered <- beaum2 %>%
  filter(!PollutionNameCorrected %in% auth_pest$PollutionNameCorrected[auth_pest$Authorised == "No"], PollutionNameCorrected != "1,3 - dichloropropene; chloropicrin")    
# remove related mixture

table(beaum2_filtered$TaxaGroup)
table(beaum2_filtered$Measurement)

# Need to convert SE to SD
beaum2_filtered <- beaum2_filtered %>%
  mutate(
    # Convert Control errors to SD
    Control_SD_new = case_when(
      Error == "SD" ~ Control_SD,
      Error == "SE" ~ Control_SD * sqrt(Control_N),
      TRUE ~ NA_real_
    ),
    
    # Convert Treatment errors to SD
    Treatment_SD_new = case_when(
      Error == "SD" ~ Treatment_SD,
      Error == "SE" ~ Treatment_SD * sqrt(Treatment_N),
      TRUE ~ NA_real_
    )
  )

# compute log response ratios
beaum_lrr <- escalc(
  measure = "ROM",      # Ratio of Means = exp(lnRR)
  m1i = Treatment_mean,
  sd1i = Treatment_SD_new,
  n1i = Treatment_N,
  m2i = Control_mean,
  sd2i = Control_SD_new,
  n2i = Control_N,
  data = beaum2_filtered
)

# check numbers of studies
summary_lrr1 <- beaum_lrr %>%
  group_by(BiodivMetric, PollutantClass) %>%
  summarise(
    n_studies = n_distinct(ID),   # Anzahl unterschiedlicher Studien
    n_observations = n(),         # Anzahl an Reihen / Fällen
    .groups = "drop"
  )

# Conversion for aggregation
beaum_lrr <- beaum_lrr %>%
			 mutate(PollutantClass_new = recode(PollutantClass,
			 "Pesticide (broad spectrum)" = "Insecticide",
			 "Herbicides,fungicides" = "Herbicide"))

summary_lrr2 <- beaum_lrr %>%
  group_by(BiodivMetric, PollutantClass_new) %>%
  summarise(
    n_studies = n_distinct(ID),   # Anzahl unterschiedlicher Studien
    n_observations = n(),         # Anzahl an Reihen / Fällen
    .groups = "drop"
  )

beaum_lrr <- beaum_lrr %>%
  mutate(
    BiodivMetric_new = case_when(
      BiodivMetric == "Abundance" ~ "Pop_or_Comm",
      BiodivMetric == "Diversity indices" ~ "Pop_or_Comm",
      BiodivMetric == "Biomass" ~ "Pop_or_Comm",
      BiodivMetric == "Richness" ~ "Pop_or_Comm",
      BiodivMetric == "Evenness indices" ~ "Pop_or_Comm",
      TRUE ~ BiodivMetric
    )
  )

summary_lrr3 <- beaum_lrr %>%
  group_by(BiodivMetric_new, PollutantClass_new) %>%
  summarise(
    n_studies = n_distinct(ID),   # Anzahl unterschiedlicher Studien
    n_observations = n(),         # Anzahl an Reihen / Fällen
    .groups = "drop"
  )

# aggregate by grouping variables
agg <- beaum_lrr %>%
  group_split(PollutantClass_new) %>%
  map_dfr(~{
    fit <- rma(yi = .x$yi, vi = .x$vi, method = "REML")
    tibble(
      Group = unique(.x$PollutantClass_new),
      yi = coef(fit),
      vi = fit$se^2
    )
  })

# compute standard deviation
agg <- agg %>%
  mutate(
    SE = sqrt(vi),
    CI_lower = yi - 1.96 * SE,
    CI_upper = yi + 1.96 * SE
  )

# 
summary_lrr3 <- beaum_lrr %>%
  group_by(BiodivMetric_new, PollutantClass_new) %>%
  summarise(
    n_studies = n_distinct(ID),   
    n_observations = n(),        
    .groups = "drop"
  )

#################################   
# Analysis for biodiversity data II  
#################################      
# needs to be downloaded locally from https://github.com/schaeferRCOHR/Meta_pesticides/blob/main/Biodiv_overview.csv
eff <- read.csv("Biodiv_overview.csv", sep = ";", dec = ".", na.string = "NA")

summary(factor(eff$Response.measure))

# summarize positive and negative effect sizes
table(factor(eff$Effect.direction))

# consider only significant ones
table(factor(eff$Effect.direction[eff$Stat_Sign_Based_on_CI == "Y"]))

table(factor(eff$Effect.direction[eff$Only_authorised == "Y"]))

# limit to cases that are plotted
eff2 <- eff %>%
			 filter(!Response.measure %in% c(
      		 "Log Odds ratio",
      		 "Raw data"), Case_ID != 8)
# remove study on microbial communities, only with %negative effect

summary(factor(eff$Response))

eff2 	%>% 
			count(Response.measure, Response)
			
# Plot for % positive, negative and no effect
effects_keep <- c("% positive effect", "% no effect", "% negative effect")

exp_effects <- eff2 %>%
  # normalize whitespace if needed
  mutate(Response.measure = stringr::str_squish(as.character(Response.measure))) %>%
  filter(Response.measure %in% effects_keep) %>%
  # make sure effect factor has the desired order (positive, no, negative)
  mutate(Response.measure = factor(Response.measure,
                                   levels = c("% positive effect", "% no effect", "% negative effect")),
         Mean_rel = Mean/100                          
                                   )

sample_sizes <- exp_effects %>%
  group_by(Response) %>%
  summarise(N = mean(Number.of.observations))

# use multinomial model to test for significant difference
model2 <- multinom(Mean_rel ~ Response, data = exp_effects)
summary(model2)
# compare to null model
null <- multinom(Mean_rel ~ 1, data = exp_effects)
anova(null, model2, test="Chisq")
    
ggplot(exp_effects, aes(x = Response, y = Mean_rel, fill = Response.measure)) +
  geom_col(color = "black", width = 0.7) +
  geom_text(aes(label = ifelse(Mean > 5, paste0(round(Mean,1), "%"), "")),
            position = position_stack(vjust = 0.5), size = 6) +
  geom_text(data = sample_sizes,
            aes(x = Response, y = -0.05,
                label = paste0("n = ", N)),
            size = 6, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(-0.1, 1)) +            
  labs(x = NULL, y = "Percent of cases", title = NULL,
       fill = "Effect type") +
  theme_minimal(base_size = 22)  +
  theme(axis.text.x = element_text(vjust = 1)) +
  scale_fill_manual(
    values = c(
      "% positive effect" = "green",
      "% no effect" = "lightgray", # oder "beige", "lightyellow"
      "% negative effect" = "red"
    )
  )



## summarise data
summary_eff_a <- eff2 %>%
			filter(Response.measure == "Natural log response ratio")

# remove community indices that only are used once
# omit growth as less sensitive parameter
eff3 <- eff2 %>%
			  filter(!Response %in% c(
      		 "Enrichment index",
      		 "Plant feeders", "Plant parasite index", 
      		 "Structure index"), 
      		 Response.measure == "Natural log response ratio", 
      		 Response != "Growth"
      		 ) %>%			 
			 mutate(Pesticide.type_new = recode(Pesticide.type, 
			 "Pesticide (except nematicide and herbicide)" = "Insecticide",
			 "Nematicide" = "Insecticide" 
			 ), Organism.group_new = recode(Organism.group,
			 "Nematodes" = "Soil invertebrates"		 
			 ))
# aggregate Nematicide and Pesticide with Insecticides

## summarise data
summary_eff_b <- eff3 %>%
  group_by(Pesticide.type_new, Response.level, Biome, Organism.group_new) %>%
  summarise(
    n_observations = n(),  
    .groups = "drop"
  )

## aggregate responses using random effects meta-analysis
eff3  <- eff3 %>%
  mutate(var = (Upper.value - Lower.value) / (2 * 1.96)^2)

agg_results <- eff3 %>%
  group_by(Pesticide.type_new, Response.level, Biome, Organism.group_new) %>%
  group_modify(~ {
    meta <- rma(yi = .x$Mean, vi = .x$var, method = "REML")
    tibble(
      pooled_effect = meta$b,
      pooled_lower  = meta$ci.lb,
      pooled_upper  = meta$ci.ub,
      k             = meta$k,
      Number.of.observations = sum(.x$Number.of.observations, na.rm = TRUE),
      Number.of.studies      = sum(.x$Number.of.studies, na.rm = TRUE)
    )
  })

# sort table
agg_results_a <-	agg_results %>%
			arrange(
    		Biome,
    		Organism.group_new,
    		Response.level,
    		Pesticide.type_new
    		) %>% mutate(across(Number.of.observations, as.character))


## figure for SETAC session
agg_results_terr <- agg_results_a %>%
					filter(Biome == "Terrestrial") %>%
					rename("mean" = "pooled_effect", "lower" = "pooled_lower", "upper" = 						"pooled_upper")

# create label text after removing a few columns
agg_results_terr_b  <-  agg_results_terr %>%
						ungroup %>%
						select(-Biome, -k, -Number.of.studies)
						
# add first row with names
agg_results_terr_c <- add_row(agg_results_terr_b, Pesticide.type_new = "Pesticide type", Response.level = "Response level", Organism.group_new = "Organism group", mean = NA , lower = NA, upper = NA, Number.of.observations = "n", .before = 1)

labeltext <- agg_results_terr_c[ , c("Organism.group_new", "Response.level", "Pesticide.type_new", "Number.of.observations")]

data <- agg_results_terr_c[ , c("mean", "lower", "upper")]

forestplot(
  labeltext = labeltext,
  x = data, 
  is.summary = c(TRUE, rep(FALSE, nrow(labeltext)-1)),
  zero = 0,
  xlab = "Effect size (LRR)",
  boxsize = 0.2,
  lineheight = "auto",
  colgap = unit(6, "mm"),
  col = fpColors(box = "royalblue", line = "darkblue", summary = "black"),
  lwd.ci = 2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.1,
  xticks = seq(-1, 1.5, 0.5),
  txt_gp = fpTxtGp(xlab  = gpar(cex = 1.5),  # x-axis label
    ticks = gpar(cex = 1.2),  # tick labels
    label = gpar(cex = 1)     # left-hand labels
  )
)

## create plot for aquatic and terrestrial together
agg_results_b <- agg_results_a %>%
					rename("mean" = "pooled_effect", "lower" = "pooled_lower", "upper" = 						"pooled_upper") %>%
						ungroup %>%
						select(-k, -Number.of.studies)
						
# add first row with names
agg_results_c <- add_row(agg_results_b, Pesticide.type_new = "Pesticide type", Response.level = "Response level", Biome = "Biome", Organism.group_new = "Organism group", mean = NA , lower = NA, upper = NA, Number.of.observations = "n", .before = 1)

labeltext_a <- agg_results_c[ , c("Biome", "Organism.group_new", "Response.level", "Pesticide.type_new", "Number.of.observations")]

data_a <- agg_results_c[ , c("mean", "lower", "upper")]

forestplot(
  labeltext = labeltext_a,
  x = data_a, 
  is.summary = c(TRUE, rep(FALSE, nrow(labeltext_a)-1)),
  zero = 0,
  xlab = "Effect size (LRR)",
  boxsize = 0.2,
  lineheight = "auto",
  colgap = unit(6, "mm"),
  col = fpColors(box = "royalblue", line = "darkblue", summary = "black"),
  lwd.ci = 2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.1,
  xticks = seq(-1, 1.5, 0.5),
  txt_gp = fpTxtGp(xlab  = gpar(cex = 1.5),  # x-axis label
    ticks = gpar(cex = 1.2),  # tick labels
    label = gpar(cex = 1)     # left-hand labels
  )
)

## run mixed effects meta-analysis to evaluate influence of moderators
# compute variance first
data_rma <-  agg_results_a %>%
  mutate(var = (pooled_upper - pooled_lower) / (2 * 1.96)^2, 
         across(Number.of.observations, as.numeric))

# run model
model_rma <- rma(yi = pooled_effect, vi = var, ni = Number.of.observations, mods = ~ Biome + Organism.group_new * Pesticide.type_new, data = data_rma)
summary(model_rma)
anova(model_rma)

# drop interaction and organism group
model_rma2 <- rma(yi = pooled_effect, vi = var, ni = Number.of.observations, mods = ~ Biome + Pesticide.type_new, data = data_rma)
summary(model_rma2)

# drop interaction and organism group
model_rma3 <- rma(yi = pooled_effect, vi = var, ni = Number.of.observations, mods = ~ Biome, data = data_rma)
summary(model_rma3)






  
      