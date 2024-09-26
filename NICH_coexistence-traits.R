#*********************************************************************************
# R codes used in: Functional traits predict outcomes of current and novel competition under warmer climate
# Author: Shengman Lyu (shengman.lyu@gmail.com)
# Date: 13.9.2024
# R codes used for population modelling and coexstence analyses can be found in:  Lyu, S. and J. M. Alexander (2023). "Compensatory responses of vital rates attenuate impacts of competition on population growth and promote coexistence." Ecology Letters 26(3): 437-447.
#*********************************************************************************
rm(list=ls())

library(MuMIn)
library(car)
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(vegan)
library(lme4)
library(car)
library(lmerTest)
library(ggradar)

#****************************************************
#  ---- 1. Functions  ----
#****************************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot ND and RFD (log-scale) frame indicating coexistence outcomes ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
coex.frame <- function(x1=-0.5,x2=0.9) {
  # x1, x2 lower and upper limits of ND, in which x1 < -0.001
  x <- c(seq(x1,-0.001, length.out = 100), seq(0.001,x2,length.out = 100))
  y <- log(1/(1-x))
  p <- ggplot() + 
    geom_line(data=data.frame(x=x, y=y), aes(x=x, y=y), inherit.aes = FALSE) +
    geom_line(data=data.frame(x=x, y=-y), aes(x=x, y=y), inherit.aes = FALSE) +
    geom_line(data=data.frame(x=x, y=rep(0,length(x))), aes(x=x, y=y), inherit.aes = FALSE, linetype="dashed")
  return(p)
}
coex.frame()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add frame
frame.coex = function(x1 = -0.01, x2= 0.91) {
  # x1, x2 lower and upper limits of ND, in which x1 < -0.001
  x <- c(seq(x1,-0.001, length.out = 100), seq(0.001,x2,length.out = 100))
  y <- log(1/(1-x))
  return(data.frame(x=x,y=y))
}
frame.coex()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add polygan
# polygan of coexistence area
polygan.coex = function(x1 = 0.0001, x2= 0.91) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = rep(x2,100)
  xx3 = seq(x2,x1,length.out = 100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(log(1/(1-xx1)), log(1/(1-xx2)),-log(1/(1-xx3)))
  return(data.frame(x=xx,y=yy))
}
ggplot(polygan.coex(), aes(x=x,y=y)) +
  geom_polygon()

polygan.coex.no = function(x1 = -1, x2=0) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = seq(x2,x1,length.out = 100)
  xx3 = rep(x1,100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(-xx1, xx2, seq(x1,-x1, length.out=100))
  return(data.frame(x=xx,y=yy))
}
ggplot(polygan.coex.no(), aes(x=x,y=y)) +
  geom_polygon()

polygan.coex.fd.abs = function(x1 = -1, x2=0) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = seq(x2,x1,length.out = 100)
  xx3 = rep(x1,100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(-xx1, rep(0,100), seq(0,-x1, length.out=100))
  return(data.frame(x=xx,y=yy))
}


polygan.prio = function(x1 = -3, x2= -0.001) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = rep(-0.001,100)
  xx3 = seq(x2,x1,length.out = 100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(log(1/(1-xx1)), log(1/(1-xx2)),-log(1/(1-xx3)))
  return(data.frame(x=xx,y=yy))
}
ggplot(polygan.prio(), aes(x=x,y=y)) +
  geom_polygon()

polygan.prio.no = function(x1 = 0, x2=1) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = rep(x2,100)
  xx3 = seq(x2,x1,length.out = 100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(xx1,seq(x2,0, length.out=100), -xx3)
  return(data.frame(x=xx,y=yy))
}
ggplot(polygan.prio.no(), aes(x=x,y=y)) +
  geom_polygon()

polygan.prio.fd.abs = function(x1 = 0, x2=1) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = rep(x2,100)
  xx3 = seq(x2,x1,length.out = 100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(xx1,seq(x2,0, length.out=100), rep(0, 100))
  return(data.frame(x=xx,y=yy))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to extract statics from lmer ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
coef.lmer <- function(model) {
  av <- as.data.frame(Anova(model))
  r2 <- r.squaredGLMM(model)
  factors <- rownames(as.data.frame(Anova(model)))
  pf <- NULL
  for(i in factors) {
    pf <- c(pf, as.numeric(av[i,]))
  }
  pf <- c(pf, r2[1], r2[2])
  names(pf) <- c(paste(rep(factors, each=3), c("F", "d.f.", "P")), c("R2.marginal", "R2.conditional"))
  pf
}

# to extract slopes
extract.lmer <- function(model) {
  slope <- fixef(model)[2]
  p.value <- as.numeric(summary(model)$coef[2,5])
  r2 <- r.squaredGLMM(model)
  aic <- AICc(model)
  pf <- c(slope, p.value, r2[1], r2[2], aic)
  names(pf) <- c("Slope", "P value","R2 marginal", "R2 conditional", "AICc")
  pf
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to extract statics from lm ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# to extract AIC
extract.lm <- function(model) {
  slope <- coef(model)[2]
  sm <- summary(model)
  p.value <- as.numeric(sm$coef[2,4])
  r2 <- sm$r.square
  r2.adj <- sm$adj.r.squared
  aic <- AICc(model)
  pf <- c(slope, p.value, r2, r2.adj, aic)
  names(pf) <- c("Slope", "P value","R2", "R2 sdjusted", "AICc")
  pf
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compuate confidence interval using linear model ----
# same as: Hmisc::mean_cl_normal
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci <- function(x) {
  fit <- lm(x~1)
  c(mean(x), confint(fit)[1, ])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compuate SE-based confidence interval ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci_se <- function(x) {
  u = mean(x)
  se = sqrt(var(x)/length(x))
  data.frame(y=u, ymin=u-1.96*se, ymax=u+1.96*se)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compuate SD-based confidence interval ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci_sd <- function(x, ...) {
  u = mean(x, ...)
  sd = sqrt(var(x, ...))
  data.frame(y=u, ymin=u-1.96*sd, ymax=u+1.96*sd)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compuate mean, sd, cv, min, man and range ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_cv <- function(x, ...) {
  u = mean(x, ...)
  sd = sqrt(var(x, ...))
  cv = abs(sd/u)
  min = range(x, ...)[1]
  max = range(x, ...)[2]
  rg = max - min
  data.frame(n = length(x[!is.na(x)]), mean=u, sd=sd, cv=cv, range = rg, min=min, max=max)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compuate mean, sd, cv, min, man and range ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_sd <- function(x, ...) {
  u = mean(x, ...)
  sd = sqrt(var(x, ...))
  min = u - sd
  max = u + sd
  data.frame(n = length(x[!is.na(x)]), mean=u, sd=sd, min=min, max=max)
}

mean_sd.ggplot <- function(x) {
  x <- stats::na.omit(x)
  u = mean(x)
  sd = sqrt(var(x))
  min = u - sd
  max = u + sd
  data_frame(y = u, ymin = min, ymax = max)
  #new_data_frame(list(y = u, ymin = min, ymax = max), n = 1)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compuate percentile-based confidence interval ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci_quantile <- function(x, ...) {
  u = mean(x, ...)
  x.min = quantile(x, probs = 0.025, ...)
  x.max = quantile(x, probs = 0.975, ...)
  data.frame(y=u, ymin=x.min, ymax=x.max) 
}

median_ci_quantile <- function(x, ...) {
  u = median(x, ...)
  x.min = quantile(x, probs = 0.025, ...)
  x.max = quantile(x, probs = 0.975, ...)
  data.frame(y=u, ymin=x.min, ymax=x.max) 
}

#****************************************************
# ---- 2. PCA ----
# to show species distribtuion in fucntional space
# to derive integrative traits
#****************************************************
#************************************
# ** - trait data ----
#************************************
# all traits without gaps
traits <- read_excel("Trait.csv", na= "NA")
traits

# remove NA (Armo)
trait.nona <- traits[apply(is.na(traits),1,sum) == 0,]

# trait list
#trait.list <- read_excel("/Users/slyu1/LVSM/NICH/Analysis/NICH-CH3/Trait list.xlsx")
trait.list <- read_excel("Trait list.xlsx")

# remove NA (Armo)
trait.nona <- traits[apply(is.na(traits),1,sum) == 0, ]

# 15 traits
trait.id <- trait.list$trait.ID[!is.na(trait.list$trait.ID)]

# trait data transformation
trait.transform <- trait.nona
for(i in trait.id) {
  # i = trait.id[1]
  
  # trait data
  trait.i <- trait.nona[, i]
  
  # transform trait data
  transform.i <- trait.list[trait.list$trait.ID == i & !is.na(trait.list$trait.ID), ]$Transformation
  if(transform.i == "log") {
    trait.transform.i = log(trait.i) 
  }
  if(transform.i == "sqrt") {
    trait.transform.i = sqrt(trait.i) 
  }
  if(transform.i == "none") {
    trait.transform.i = trait.i
  }
  
  # trait
  trait.transform[, i] <- trait.transform.i
}

# scale all tratis
trait.scale <- as.data.frame(sapply(trait.transform[trait.id],scale))

colnames(trait.scale) <- trait.list$Abbreviation[!is.na(trait.list$trait.ID)]

#************************************
# ** - PCA ----
#************************************

#************************************
# PCA: variable in rows or columns?
# variables in columns
# PC1: 22.78%, PC2 22.30%, PC3 15.43%
trait.pca <- prcomp(trait.scale)
trait.pca

summary(trait.pca)
trait.pca$sdev
trait.pca$rotation
trait.pca$x

#************************************
#  plot PC axis
pcs <- data.frame(species = trait.nona$species, 
                  origin = trait.nona$origin, 
                  site = trait.nona$site, 
                  elevation = trait.nona$elevation, 
                  PC1 = trait.pca$x[,"PC1"],
                  PC2 = trait.pca$x[,"PC2"],
                  PC3 = trait.pca$x[,"PC3"])
pcs
# PC 1 and 2 loadings
# loading
fig_pc12_loadings <- data.frame(trait.id = rownames(trait.pca$rotation), xx = 0, yy= 0, xx.end = trait.pca$rotation[,"PC1"], yy.end = trait.pca$rotation[,"PC2"]) %>%
  ggplot(aes(x=xx, y= yy, xend=xx.end, yend=yy.end, label=trait.id))  +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_segment(arrow = arrow(length = unit(0.1, "inches")), col="black", alpha=0.8) +
  geom_text(data = data.frame(xx.end = trait.pca$rotation[,"PC1"], 
                              yy.end = trait.pca$rotation[,"PC2"],
                              trait.id = rownames(trait.pca$rotation)),
            aes(x=xx.end, y = yy.end, ), col="black", alpha=0.8) +
  scale_x_continuous(name = "PC 1 (22.78%)") +
  scale_y_continuous(name = "PC 2 (22.30%)")
fig_pc12_loadings

# PC 1 and 2 for sites and origins
fig_pc12_site <- pcs %>%
  mutate(site = dplyr::recode(site, "Les Posses"= "Low site", "Solalex" = "Middle site", "Anzeindaz"= "High site")) %>%
  mutate(site = factor(site, levels=c("Low site", "Middle site", "High site"))) %>%
  mutate(origin = dplyr::recode(origin, "Highland"= "Highland species", "Lowland" = "Lowland species")) %>%
  mutate(origin = factor(origin, levels=c("Lowland species", "Highland species"))) %>%
  ggplot(aes(x=PC1, y = PC2, label=species, col=site)) +
  geom_point(data=filter(pcs, origin=="Highland"), aes(x=PC1, y=PC2), size=12, shape=1, color="grey70") +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_text(fontface='bold', size=4, show.legend = FALSE) +
  scale_alpha_manual (values=c("Lowland species" = 1, "Highland species" = 0.6)) +
  scale_color_manual(values=c("Low site" = "#FB9A06FF", "Middle site" = "#6DCD59FF", "High site" = "#3E4A89FF")) +
  scale_x_continuous(name = "PC 1 (22.78%)") +
  scale_y_continuous(name = "PC 2 (22.30%)")
fig_pc12_site

# PC 1 and 3 loadings
fig_pc13_loadings <- data.frame(trait.id = rownames(trait.pca$rotation), xx = 0, yy= 0, xx.end = trait.pca$rotation[,"PC1"], yy.end = trait.pca$rotation[,"PC3"]) %>%
  ggplot(aes(x=xx, y= yy, xend=xx.end, yend=yy.end, label=trait.id))  +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_segment(arrow = arrow(length = unit(0.1, "inches")), col="black", alpha=0.8) +
  geom_text(data = data.frame(xx.end = trait.pca$rotation[,"PC1"], 
                              yy.end = trait.pca$rotation[,"PC3"],
                              trait.id = rownames(trait.pca$rotation)),
            aes(x=xx.end, y = yy.end, ), col="black", alpha=0.8) +
  scale_x_continuous(name = "PC 1 (22.78%)") +
  scale_y_continuous(name = "PC 3 (15.43%)")
fig_pc13_loadings

# PC 1 and 3 site
fig_pc13_site <- pcs %>%
  mutate(site = dplyr::recode(site, "Les Posses"= "Low site", "Solalex" = "Middle site", "Anzeindaz"= "High site")) %>%
  mutate(site = factor(site, levels=c("Low site", "Middle site", "High site"))) %>%
  mutate(origin = dplyr::recode(origin, "Highland"= "Highland species", "Lowland" = "Lowland species")) %>%
  mutate(origin = factor(origin, levels=c("Lowland species", "Highland species"))) %>%
  ggplot(aes(x=PC1, y = PC3, label=species, col=site)) +
  geom_point(data=filter(pcs, origin=="Highland"), aes(x=PC1, y=PC3), size=12, shape=1, color="grey70") +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_text(fontface='bold', size=4, show.legend = FALSE) +
  scale_alpha_manual (values=c("Lowland species" = 1, "Highland species" = 0.6)) +
  scale_color_manual(values=c("Low site" = "#FB9A06FF", "Middle site" = "#6DCD59FF", "High site" = "#3E4A89FF")) +
  scale_x_continuous(name = "PC 1 (22.78%)") +
  scale_y_continuous(name = "PC 3 (15.43%)")
fig_pc13_site

#************************************
# ** - combine PC axises to trait data ----
#************************************
histogram(pcs$PC1)
histogram((pcs$PC2)) 
histogram((pcs$PC3))
pcs.scale <- as.data.frame(sapply(pcs[,c("PC1", "PC2", "PC3")],scale))
trait.scale <- cbind(trait.nona[,1:4], trait.scale, pcs.scale)

#*****************************************************************************
# ---- 3. Traits between species and across the elevation  ----
#*****************************************************************************
# scaled traits
trait.scale
trait.ID <- colnames(trait.scale)[5:22]

# species
sps <- c("Anal", "Armo", "Asal", "Plal", "Poal", "Seca", "Trba",  # 7 alpine
         "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr")  # 7 lowland

#************************************
# ** - site level traits between species ----
#************************************
# Convert to long-format
trait.scale.long <- trait.scale  %>%
  pivot_longer(cols=trait.ID, names_to="trait.id", values_to="trait.value") %>%
  mutate(species = factor(species, levels=c("Anal","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low site", "Solalex" = "Middle site", "Anzeindaz" = "High site")) %>%
  mutate(site = factor(site, levels=c("Low site", "Middle site", "High site"))) %>%
  mutate(trait.id = dplyr::recode(trait.id, "RH" = "Reproductive height (RH)", 
                                  "VH" = "Vegetative height (VH)",
                                  "LA" = "Leaf area (LA)",
                                  "LDMC" = "Leaf dry mass content (LDMC)",
                                  "SLA" = "Specific leaf area (SLA)", 
                                  "LN" = "Leaf nitrogen content (LN)",
                                  "LC" = "Leaf carbon content (LC)",
                                  "δ13C" = "Carbon isotope ratio (δ13C)",
                                  "RMF" = "Root-mass fraction (RMF)",
                                  "FRP" = "Fine root proportion (FRP)",
                                  "FRD" = "Fine root density (FRD)", 
                                  "SRL" = "Specific root length (SRL)",
                                  "SM" = "Seed mass (SM)", 
                                  "FFD" = "First flowering date (FFD)",
                                  "Light" = "Light interception (Light)",
                                  "PC1" = "PC1", 
                                  "PC2" = "PC2", 
                                  "PC3" = "PC3")) %>%
  mutate(trait.id = factor(trait.id, levels = c("Reproductive height (RH)",
                                                "Vegetative height (VH)",
                                                "Leaf area (LA)",
                                                "Leaf dry mass content (LDMC)",
                                                "Specific leaf area (SLA)", 
                                                "Leaf nitrogen content (LN)",
                                                "Leaf carbon content (LC)",
                                                "Carbon isotope ratio (δ13C)",
                                                "Root-mass fraction (RMF)",
                                                "Fine root proportion (FRP)",
                                                "Fine root density (FRD)", 
                                                "Specific root length (SRL)",
                                                "Seed mass (SM)", 
                                                "First flowering date (FFD)",
                                                "Light interception (Light)",
                                                "PC1", 
                                                "PC2", 
                                                "PC3")))

# between species
fig.trait_species <-  trait.scale.long %>% 
  ggplot(aes(x=species, y = trait.value, col=site)) +
  geom_hline(yintercept = 0, linetype="dotted", col="grey") +
  geom_point() +
  #stat_summary(fun.data = mean_sd, col="black", alpha=0.6)  +
  facet_wrap(~trait.id, ncol=3, scales="free_x") +
  scale_color_manual(values=c("Low site" = "#FB9A06FF", "Middle site" = "#6DCD59FF", "High site" = "#3E4A89FF")) +
  scale_x_discrete(name = "Species") +
  scale_y_continuous(name = "Trait value (standardised)")
fig.trait_species

# test
trait.scale.long %>%
  split(.$trait.id) %>%
  purrr::map(.x = ., .f = lmer, formula = trait.value ~ site*origin + (1|species)) %>%
  purrr::map(Anova)
  #purrr::map(rand)
  #purrr::map(summary)

# 4 traits differ significantly between sites
# 1 trait differ significantly between lowland and highland species
# all traits differ significantly between species

#************************************
# ** - correlation between traits ----
#************************************
trait.cor <- cor(trait.scale[trait.ID])
trait.cor[upper.tri(trait.cor, diag = TRUE)] <- NA
trait.cor
#write.csv(trait.cor, "/Users/slyu/LVSM/NICH/Results/coexistence_traits/trait.correlation.csv")

#*****************************************************************************
# ---- 4. Coexistence ~ sites ----
#*****************************************************************************

#************************************
# ** - coexistence data ----
#************************************
outcome <- read_csv("Coexistence.csv", na="NA")
outcome <- filter(outcome, bootstrap == 0)

outcome$sps1.site <- paste(outcome$sps1, outcome$site)
outcome$sps2.site <- paste(outcome$sps2, outcome$site)

# how many pairs
outcome %>%
  filter(!is.na(nd)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  .$site %>%
  table()

# NDFD
outcome %>%
  filter(!is.na(nd)) %>%
  ggplot(aes(x=log(1-nd), y=log(fd), col=site)) +
  geom_point()

# valid pairs across the three site
pair.table <- outcome %>%
  filter(!is.na(nd)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  .$pair %>% 
  table()

# common pairs
pair.common <- names(pair.table)[pair.table > 2]
#pair.common <- names(pair.table)[pair.table == 3]

#************************************
# ** - Outcomes of competition ~ sites ----
#************************************
# coexistence plot
fig.coexistence_site <- outcome %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel", "Current")) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low site", "Solalex" = "Middle site", "Anzeindaz" = "High site")) %>%
  mutate(site = factor(site, levels=c("Low site", "Middle site", "High site"))) %>%
  filter(!is.na(nd)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  # common pairs
  #filter(pair %in% pair.common) %>%
  mutate(NO.log = log(1-nd)) %>%
  mutate(RFD.log = log(fd)) %>% 
  ggplot(aes(x = NO.log, y = abs(RFD.log), col=is.novel)) +
  geom_polygon(data=polygan.coex.fd.abs(x1=-3), aes(x=x,y=y), col="grey30",fill="grey90", inherit.aes = FALSE) + 
  geom_polygon(data=polygan.prio.fd.abs(x2=3), aes(x=x,y=y), col="grey30", fill="grey90", inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", col="grey30") +
  geom_hline(yintercept = 0, linetype = "dashed", col="grey30") +
  geom_point() +
  coord_cartesian(xlim=c(-2.7,2.5), ylim=c(NA, 2)) +
  scale_color_manual(values=c("orange", "black")) +
  scale_x_continuous(name = "ln(Niche overlap)") +
  scale_y_continuous(name = "ln(Relative fitness difference)", expand = c(0, 0)) +
  facet_wrap(~site)
fig.coexistence_site

#************************************
# ** - summarise outcomes ----
#************************************
# all sites
outcome %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel", "Current")) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low site", "Solalex" = "Middle site", "Anzeindaz" = "High site")) %>%
  mutate(site = factor(site, levels=c("Low site", "Middle site", "High site"))) %>%
  filter(!is.na(nd)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  .$outcome.ndfd %>%
  table()
# 61/(61+32)

# current vs novel
outcome %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel", "Current")) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low site", "Solalex" = "Middle site", "Anzeindaz" = "High site")) %>%
  mutate(site = factor(site, levels=c("Low site", "Middle site", "High site"))) %>%
  filter(!is.na(nd)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  filter(is.novel == "Novel") %>%
  .$outcome.ndfd %>%
  table()
# 35/(35+12)
# 26/(26+20)


# low site
outcome %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel", "Current")) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low site", "Solalex" = "Middle site", "Anzeindaz" = "High site")) %>%
  mutate(site = factor(site, levels=c("Low site", "Middle site", "High site"))) %>%
  filter(site == "Low site") %>%
  filter(!is.na(nd)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  .$outcome.ndfd %>%
  table()
# 16/(16+13)

# middle site
outcome %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel", "Current")) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low site", "Solalex" = "Middle site", "Anzeindaz" = "High site")) %>%
  filter(site == "Middle site") %>%
  filter(!is.na(nd)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  .$outcome.ndfd %>%
  table()
# 23/(23+5)

# high site
outcome %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel", "Current")) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low site", "Solalex" = "Middle site", "Anzeindaz" = "High site")) %>%
  mutate(site = factor(site, levels=c("Low site", "Middle site", "High site"))) %>%
  filter(site == "High site") %>%
  filter(!is.na(nd)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  .$outcome.ndfd %>%
  table()
# 22/(22+14)

#************************************
# ** - ND, RFD ~ sites ----
#************************************
# ND ~ site plot
fig.ND_site.novel <- outcome %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel", "Current")) %>%
  filter(!is.na(nd)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  # common pairs
  filter(pair %in% pair.common) %>%
  mutate(NO.log = log(1-nd)) %>%
  ggplot(aes(x=site, y = NO.log, col=is.novel)) +
  geom_point(position = position_dodge(width=0.5)) +
  #stat_summary(fun.data=mean_sd.ggplot, position = position_dodge(width=0.5)) +
  geom_boxplot(position = position_dodge(width=0.5), width=0.5, alpha = 0.5, outlier.alpha = 0) +
  scale_color_manual(values=c("orange", "black")) +
  scale_x_discrete(name = "Sites") +
  scale_y_continuous(name = "ln(Niche overlap)")
fig.ND_site.novel

# ND ~ site test
# species has smaller NO at the high site
# current competitors have smaller NO
outcome %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel", "Current")) %>%
  filter(!is.na(nd)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  # common pairs
  filter(pair %in% pair.common) %>%
  mutate(NO.log = log(1-nd)) %>%
  lmer(NO.log ~ site * is.novel + (1|sps1) + (1|sps2), data = .) %>%
  #summary() 
  Anova()

# RFD ~ site plot
fig.RFD_site.novel <- outcome %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel", "Current")) %>%
  filter(!is.na(nd)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  # common pairs
  filter(pair %in% pair.common) %>%
  mutate(RFD.log = log(fd)) %>% 
  ggplot(aes(x=site, y = abs(RFD.log), col=is.novel)) +
  geom_point(position = position_dodge(width=0.5)) +
  #stat_summary(fun.data=mean_sd.ggplot, position = position_dodge(width=0.5)) +
  geom_boxplot(position = position_dodge(width=0.5), width=0.5, alpha = 0.5, outlier.alpha = 0) +
  scale_color_manual(values=c("orange", "black")) +
  scale_x_discrete(name = "Sites") +
  scale_y_continuous(name = "ln(Relative fitness difference)")
fig.RFD_site.novel

# RFD ~ site test
outcome %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel", "Current")) %>%
  filter(!is.na(nd)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  # common pairs
  filter(pair %in% pair.common) %>%
  mutate(RFD.log = log(fd)) %>% 
  lmer(abs(RFD.log) ~ site * is.novel + (1|sps1) + (1|sps2), data = .) %>%
  #summary()
  Anova()

#*****************************************************************************
# ---- 5. Invasion growth rates ~ trait differences ----
#*****************************************************************************
# data
trait.scale 
trait.ID <- colnames(trait.scale)[5:22]

# species
sps <- c("Anal", "Armo", "Asal", "Plal", "Poal", "Seca", "Trba",  # 7 alpine
         "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr")  # 7 lowland

#************************************
# ** - population growth rates ----
#************************************
pgr <- read_csv("Population growth.csv", na="NA")
pgr <- filter(pgr, bootstrap == 0)

# how many pairs
pgr %>%
  filter(pair == "yes") %>%
  #filter(focal.species != "Armo") %>%
  filter(!(background.species %in% c("site", "site_none", "species", "species_none"))) %>%
  filter(!is.na(pgr)) %>%
  filter(background.species != "none") %>%
  mutate(p = paste(focal.species, background.species)) %>%
  .$p %>%
  unique()

#************************************
# ** - calculate response ratio ----
#************************************
# match intrinsic lambdas
pgr.intrinsic <- filter(pgr, background.species == "none")
pgr.intrinsic$ID_focal.species <- paste(pgr.intrinsic$focal.species, pgr.intrinsic$site)
pgr$ID_focal.species <- paste(pgr$focal.species, pgr$site)
pgr$pgr.intrinsic.focal <- pgr.intrinsic[match(pgr$ID_focal.species, pgr.intrinsic$ID_focal.species), ]$pgr
pgr

# log response ratio
pgr$response.ratio <- pgr$pgr/pgr$pgr.intrinsic.focal
hist(pgr$response.ratio)
hist(log(pgr$response.ratio))

# how many pairs?
pgr %>%
  filter(focal.species != "Armo") %>%
  filter(focal.species != background.species) %>%
  filter(background.species != "none") %>%
  filter(background.species != "site") %>%
  filter(background.species != "site_none") %>%
  filter(background.species != "species") %>%
  filter(background.species != "species_none") %>%
  filter(!is.na(response.ratio)) %>%
  .$site %>%
  table()

# Intrinsic growth rate of background species
pgr.intrinsic$ID_focal.species <- paste(pgr.intrinsic$focal.species, pgr.intrinsic$site)
pgr$ID_background.species <- paste(pgr$background.species, pgr$site)
pgr$pgr.intrinsic.background <- pgr.intrinsic[match(pgr$ID_background.species, pgr.intrinsic$ID_focal.species), ]$pgr
pgr

#************************************
# ** - Response ratio across sites ----
# competition (after excluding facilitation) decrease with increasing elevation overall
# competition is the same between current vs novel competitors
#************************************

#************************************
# mean and SD
fig.rr.distribution <- pgr %>%
  filter(!(background.species %in% c("none", "site", "site_none", "species", "species_none"))) %>%
  filter(!is.na(pgr)) %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(focal.species != background.species) %>% # remove intra pairs
  filter(response.ratio != 0) %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland_Highland", "Novel", "Current")) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low", "Solalex" = "Middle", "Anzeindaz" = "High")) %>%
  mutate(site = factor(site, levels=c("Low", "Middle", "High"))) %>%
  ggplot(aes(x = log(response.ratio), fill=site, col=site)) +
  #facet_wrap(~site, ncol=1) + 
  scale_color_manual(values=c("Low" = "#FB9A06FF", "Middle" = "#6DCD59FF", "High" = "#3E4A89FF")) +
  scale_fill_manual(values=c("Low" = "#FB9A06FF", "Middle" = "#6DCD59FF", "High" = "#3E4A89FF")) +
  geom_histogram(alpha=0.6)
fig.rr.distribution

# response ratio
fig.rr_site.novel <- pgr %>%
  filter(!(background.species %in% c("none", "site", "site_none", "species", "species_none"))) %>%
  filter(!is.na(pgr)) %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(focal.species != background.species) %>% # remove intra pairs
  mutate(is.novel = ifelse(origin.pair == "Lowland_Highland", "Novel interactions", "Current interactions")) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low", "Solalex" = "Middle", "Anzeindaz" = "High")) %>%
  mutate(site = factor(site, levels=c("Low", "Middle", "High"))) %>% 
  ggplot(aes(x=site, y=log(response.ratio), col=is.novel)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_point(position = position_dodge(width=0.4))+
  geom_boxplot(alpha=0.4, width=0.4, outlier.alpha = 0) +
  scale_color_manual(values = c("black", "#FB9A06FF")) +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "ln(Response ratio)")
fig.rr_site.novel

# test
pgr %>%
  filter(!(background.species %in% c("none", "site", "site_none", "species", "species_none"))) %>%
  filter(!is.na(pgr)) %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(focal.species != background.species) %>% # remove intra pairs
  mutate(is.novel = ifelse(origin.pair == "Lowland_Highland", "Novel", "Current")) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low", "Solalex" = "Middle", "Anzeindaz" = "High")) %>%
  mutate(site = factor(site, levels=c("Low", "Middle", "High"))) %>% 
  filter(log(response.ratio) < 0) %>%
  lmer(log(response.ratio) ~ site *is.novel + (1|focal.species) + (1|background.species), data=.) %>%
  Anova()

#************************************
# ** - calculate trait differences ----
#************************************
trait.pgr.focal <-  trait.pgr.bg <- td.pgr.hie <-  td.pgr.abs <- data.frame(x=pgr$pgr)
for(i in trait.ID) {
  # i = trait.ID[1]
  
  # trait data
  fc.site.match <- match(paste(pgr$focal.species, pgr$site), paste(trait.scale$species,trait.scale$site))
  bg.site.match <- match(paste(pgr$background.species, pgr$site), paste(trait.scale$species,trait.scale$site))
  
  # hierarchical and absolute trait dissimilarity
  td.hie.focal.i <- trait.scale[fc.site.match, i]
  td.hie.bg.i <- trait.scale[bg.site.match, i]
  td.hie.i <- trait.scale[fc.site.match, i] - trait.scale[bg.site.match, i]
  td.abs.i <- abs(td.hie.i)
  
  # cbind
  trait.pgr.focal <- cbind(trait.pgr.focal, td.hie.focal.i)
  trait.pgr.bg <- cbind(trait.pgr.bg, td.hie.bg.i)
  td.pgr.hie <- cbind(td.pgr.hie, td.hie.i)
  td.pgr.abs <- cbind(td.pgr.abs, td.abs.i)
}

# 
trait.pgr.focal <- trait.pgr.focal[-1]
trait.pgr.bg <- trait.pgr.bg[-1]
td.pgr.hie <- td.pgr.hie[-1]
colnames(td.pgr.hie) <- trait.ID

# absolute differences
td.pgr.abs <- abs(td.pgr.hie)

#************************************
# ** - Combine PGR with TD ----
#************************************
# data: trait ID acronym
data_td.pgr.hie_long <- bind_cols(pgr, as_tibble(td.pgr.hie)) %>%
  filter(!(background.species %in% c("none", "site", "site_none", "species", "species_none"))) %>%
  filter(!is.na(pgr)) %>%
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  filter(focal.species != background.species) %>% # remove intra pairs
  mutate(pair.id = paste(focal.species, background.species, sep="_")) %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland_Highland", "Novel pairs", "Current pairs")) %>%
  mutate(is.novel = factor(is.novel, levels = c("Current pairs", "Novel pairs")) ) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low", "Solalex" = "Middle", "Anzeindaz" = "High")) %>%
  mutate(site = factor(site, levels=c("Low", "Middle", "High"))) %>%
  mutate(origin.pair = factor(origin.pair)) %>%
  pivot_longer(cols=colnames(td.pgr.hie), names_to="trait.id", values_to = "TD") %>%
  mutate(trait.id = factor(trait.id, levels=c("RH", "VH", "RMF",
                                              "LA", "LDMC", "SLA", "LN", "LC", 
                                              "FRP", "FRD", "SRL",
                                              "SM", 
                                              "δ13C",
                                              "FFD",
                                              "Light",
                                              "PC1", "PC2", "PC3")))

# data: trait ID
data_td.pgr.hie_long <- bind_cols(pgr, as_tibble(td.pgr.hie)) %>%
  filter(!(background.species %in% c("none", "site", "site_none", "species", "species_none"))) %>%
  filter(!is.na(pgr)) %>%
  filter(focal.species != background.species) %>% # remove intra pairs
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland_Highland", "Novel interactions", "Current interactions")) %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland_Highland", "Novel interactions", "Current interactions")) %>%
  mutate(is.novel = factor(is.novel, levels = c("Current competition", "Novel competition")) ) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low site", "Solalex" = "Middle site", "Anzeindaz" = "High site")) %>%
  mutate(site = factor(site, levels=c("Low site", "Middle site", "High site"))) %>%
  pivot_longer(cols=colnames(td.pgr.hie), names_to="trait.id", values_to = "TD") %>%
  mutate(trait.id = dplyr::recode(trait.id, 
                                  "RH" = "Reproductive height", 
                                  "VH" = "Vegetative height",
                                  "LA" = "Leaf area",
                                  "LDMC" = "Leaf dry matter content",
                                  "SLA" = "Specific leaf area", 
                                  "LN" = "Leaf nitrogen content",
                                  "LC" = "Leaf carbon content",
                                  "δ13C" = "Carbon isotope",
                                  "RMF" = "Root mass fraction",
                                  "FRP" = "Fine root proportion",
                                  "FRD" = "Fine root density", 
                                  "SRL" = "Specific root length",
                                  "SM" = "Seed mass", 
                                  "FFD" = "First flowering date",
                                  #"FO" = "Flowering overlap"
                                  "Light" = "Light interception",
                                  "PC1" = "PC1", 
                                  "PC2" = "PC2", 
                                  "PC3" = "PC3")) %>%
  mutate(trait.id = factor(trait.id, levels = c("PC3",
                                                "PC2", 
                                                "PC1", 
                                                "Light interception",
                                                "First flowering date",
                                                "Carbon isotope",
                                                "Seed mass",
                                                "Specific root length",
                                                "Fine root density", 
                                                "Fine root proportion",
                                                "Leaf carbon content",
                                                "Leaf nitrogen content",
                                                "Specific leaf area", 
                                                "Leaf dry matter content",
                                                "Leaf area",
                                                "Root mass fraction",
                                                "Vegetative height",
                                                "Reproductive height")))
data_td.pgr.hie_long

#************************************
# ** - Variation in TD  ----
#************************************
# species become more similar at high elevation
# TD decreases in mean but not in SD with increasing elevation (with all pairs)
# TD decreases in mean AND SD with increasing elevation (with only pairs with competition data)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mean and SD
fig.td.hie_site <- data_td.pgr.hie_long %>%
  filter(TD != 0) %>% # remove intra pairs
  #filter(!is.na(response.ratio)) %>% # remove pairs without competition data
  #filter(!(trait.id %in% c("PC1", "PC2", "PC3"))) %>% # exclude PC axes
  ggplot(aes(x=site, y= abs(TD), col=is.novel)) +
  #geom_point(position = position_dodge(width=0.75)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  facet_wrap(~trait.id, ncol=4, scales="free_y") +
  scale_color_manual(values=c("Current pairs" = "black",  "Novel pairs" = "#FB9A06FF")) +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Interspecific trait difference")
fig.td.hie_site

# test: all traits together
data_td.pgr.hie_long %>% 
  filter(TD != 0) %>% # remove intra pairs
  filter(!(trait.id %in% c("PC1", "PC2", "PC3"))) %>% # exclude PC axes
  lmer(abs(TD) ~ site *is.novel + (1|trait.id) + (1|focal.species) + (1|background.species), data=.) %>%
  Anova()
#summary()

# test: each trait separately
# TD differed significantly between sites
td.lmer <- data_td.pgr.hie_long %>% 
  filter(TD != 0) %>% # remove intra pairs
  split(.$trait.id) %>%
  purrr::map(.x = ., .f = lmer, formula = abs(TD) ~ site*is.novel + (1|focal.species) + (1|background.species))

# ANOVA
td.lmer %>% purrr::map(Anova, type=2)

#************************************
# ** - IGR ~ TD: lmer ----
#************************************
# lmer
igr.lmer_td.hie <- data_td.pgr.hie_long %>%
  split(.$trait.id) %>%
  #purrr::map(.x = ., .f = lmer, formula = log(pgr) ~ TD*site*is.novel + (1|focal.species) + (1|background.species))
  purrr::map(.x = ., .f = lmer, formula = log(pgr) ~ TD*site + is.novel + TD:is.novel + (1|focal.species) + (1|background.species))

igr.lmer_td.hie %>% purrr::map(Anova, type=2)
igr.lmer_td.hie %>% purrr::map(summary)
#lmer.pgr_td.hie %>% purrr::map(rand)

# summary of the lmer results
igr.lmer.summary <- 
  igr.lmer_td.hie %>% purrr::map_dfr(coef.lmer, .id="Trait")
igr.lmer.summary
#write.csv(igr.lmer.summary, "/Users/slyu1/LVSM/NICH/Results/coexistence_traits/igr.summary.lmer_novel.csv")

# lmer
igr.lmer_td.hie_low <- data_td.pgr.hie_long %>%
  filter(site == "Low") %>%
  split(.$trait.id) %>%
  purrr::map(.x = ., .f = lmer, formula = log(pgr) ~ TD + (1|focal.species) + (1|background.species))

# ANOVA
igr.lmer_td.hie_low %>%  purrr::map(Anova, type=2)

igr.lmer_td.hie_low %>% purrr::map_dfr(coef.lmer, .id="Trait")

#************************************
# ** - IGR ~ TD: slopes  ----
#************************************
# function to extract slope and confidence interval
lmer.slope <- function(site = NULL, is.novel = NULL, formula = NULL, data = NULL)  {
  
  # re-level the data
  data$site <- relevel(data$site, ref=site)
  data$is.novel <- relevel(data$is.novel, ref=is.novel)
  
  # fit the model
  mi <- lmer(formula = as.formula(formula), data=data)
  av <- Anova(mi)
  ci <- confint(mi)
  
  # slopes and CI
  intercept = as.numeric(fixef(mi)["(Intercept)"])
  slope <- as.numeric(fixef(mi)["TD"])
  slope.min <- as.numeric(ci["TD", 1])
  slope.max <- as.numeric(ci["TD", 2])
  slope.p <- as.numeric(summary(mi)$coef[2,5])
  
  # if significant
  if(slope.p < 0.05) significant.slope <- "*"
  else significant.slope <- " "
  
  # output
  out <- data.frame(intercept = intercept, slope = slope, slope.min = slope.min, slope.max = slope.max, slope.p = slope.p, significant.slope = significant.slope)
  return(out)
}

# make the table
igr.slope_site.novel <- 
  data_td.pgr.hie_long %>%
  group_by(trait.id, site, is.novel) %>%
  summarise(mean_cv(TD, na.rm=TRUE))
igr.slope_site.novel

# get slopes for each trait at each site
for(i in 1:nrow(igr.slope_site.novel)) {
  # i = 1
  di <- igr.slope_site.novel[i,]
  tr <- as.character(di$trait.id)
  st <- as.character(di$site)
  no <- as.character(di$is.novel)
  print(i)
  
  # get slope
  slope.i <- lmer.slope(site = st, is.novel = no, 
                        formula = "log(pgr) ~ TD*site + is.novel + TD:is.novel + (1|focal.species) + (1|background.species)",
                        data = filter(data_td.pgr.hie_long, trait.id==tr))
  
  # slopes and CI
  igr.slope_site.novel[i, "slope"] = slope.i$slope
  igr.slope_site.novel[i, "slope.min"] = slope.i$slope.min
  igr.slope_site.novel[i, "slope.max"] = slope.i$slope.max
  igr.slope_site.novel[i, "slope.p"] = slope.i$slope.p
  igr.slope_site.novel[i, "significant.slope"] = slope.i$significant.slope
}
igr.slope_site.novel
# write.csv(igr.slope_site.novel, "/Users/slyu1/LVSM/NICH/writing/CH3.0_traits-coexistence/Figures/IGR.slope.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot: slope ~ trait
igr.slope_site.novel[igr.slope_site.novel$trait.id == "First flowering date" & 
  igr.slope_site.novel$site == "High site" & 
  igr.slope_site.novel$is.novel == "Novel interactions", "significant.slope"] <- " "

fig_igr.slope_trait <- igr.slope_site.novel %>%
  mutate(site = dplyr::recode(site, "Low" = "Low site", "Middle" = "Middle site", "High" = "High site")) %>%
  mutate(site = factor(site, levels=c("Low site", "Middle site", "High site"))) %>%
  mutate(is.novel = factor(is.novel, levels=c("Current interactions", "Novel interactions"))) %>%
  ggplot(aes(x=slope, y = trait.id, col = is.novel, xmin=slope.min, xmax=slope.max, alpha=significant.slope)) +
  geom_vline(xintercept = 0, col="grey") +
  geom_pointrange(position = position_dodge2(width = 0.5), size=0.1) +
  coord_cartesian(xlim=c(-0.70,0.5)) +
  facet_wrap(~site, ncol=3) +
  scale_alpha_manual(values=c(0.2,1)) +
  scale_color_manual(values=c("black", "#FB9A06FF")) + 
  scale_x_continuous(name = "Effect size (slope)") +
  scale_y_discrete(name = "")
fig_igr.slope_trait

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot
fig.igr_td.hie <-  data_td.pgr.hie_long %>%
  ggplot(aes(x=TD, y= log(pgr), col=site)) +
  geom_point(alpha=0.4) +
  facet_wrap(~trait.id, scale="free_x", ncol=4) +
  scale_x_continuous(name = "Hierarchical trait difference (focal - competitor)") +
  scale_y_continuous(name = "ln(Invasion growth rate)") +
  scale_linetype_manual(values=c("dotted", "solid")) +
  scale_color_manual(values=c("Low" = "#FB9A06FF", "Middle" = "#6DCD59FF", "High" = "#3E4A89FF"))
fig.igr_td.hie

#*****************************************************************************
# ---- 6. RFD ~ trait differences ----
#*****************************************************************************

#************************************
# ** - read data ----
#************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ** - competitive outcomes ----
outcome <- read_csv("Coexistence.csv", na="NA")
outcome <- filter(outcome, bootstrap == 0)
outcome$sps1.site <- paste(outcome$sps1, outcome$site)
outcome$sps2.site <- paste(outcome$sps2, outcome$site)

# how many pairs
outcome %>%
  filter(!is.na(nd)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  .$site %>%
  table()

#************************************
# ** - calculate trait difference for each trait ----
#************************************
td.outcome.hie <- data.frame(x=outcome$igr11)
for(i in trait.ID) {
  # i = trait.ID[1]
  
  # trait data
  sps1.site.match <- match(paste(outcome$sps1, outcome$site), paste(trait.scale$species,trait.scale$site))
  sps2.site.match <- match(paste(outcome$sps2, outcome$site), paste(trait.scale$species,trait.scale$site))
  
  # hierarchical trait differences
  td.outcome.hie.i <- trait.scale[sps1.site.match, i] - trait.scale[sps2.site.match, i]
  
  # cbind
  td.outcome.hie <- cbind(td.outcome.hie, td.outcome.hie.i)
}
td.outcome.hie <- td.outcome.hie[-1]
colnames(td.outcome.hie) <- trait.ID

# absolute differences
td.outcome.abs <- abs(td.outcome.hie)

#************************************
# ** - combine outcomes with TD ----
#************************************
# TD hierarchical
# Trait abbreviation
data_outcome_td.hie <- bind_cols(outcome, as_tibble(td.outcome.hie)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low site", "Solalex" = "Middle site", "Anzeindaz"= "High site")) %>%
  mutate(site = factor(site, levels = c("Low site", "Middle site", "High site"))) %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel interactions", "Current interactions")) %>%
  mutate(is.novel = factor(is.novel, levels = c("Novel interactions", "Current interactions"))) %>%
  pivot_longer(cols=colnames(td.outcome.hie), names_to="trait.id", values_to = "TD") %>%
  mutate(trait.id = factor(trait.id, levels=c("RH", "VH", "RMF",
                                            "LA", "LDMC", "SLA", "LN", "LC", 
                                            "FRP", "FRD", "SRL",
                                            "SM", 
                                            "δ13C",
                                            "FFD",
                                            "Light",
                                            "PC1", "PC2", "PC3")))

data_outcome_td.hie

# TD absolute
# Trait abbreviation
data_outcome_td.abs <- bind_cols(outcome, as_tibble(td.outcome.abs)) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low site", "Solalex" = "Middle site", "Anzeindaz"= "High site")) %>%
  mutate(site = factor(site, levels = c("Low site", "Middle site", "High site"))) %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel interactions", "Current interactions")) %>%
  mutate(is.novel = factor(is.novel, levels = c("Novel interactions", "Current interactions"))) %>%
  pivot_longer(cols=colnames(td.outcome.abs), names_to="trait.id", values_to = "TD") %>%
  mutate(trait.id = factor(trait.id, levels=c("RH", "VH", "RMF",
                                              "LA", "LDMC", "SLA", "LN", "LC", 
                                              "FRP", "FRD", "SRL",
                                              "SM", 
                                              "δ13C",
                                              "FFD",
                                              "Light",
                                              "PC1", "PC2", "PC3")))
data_outcome_td.abs

#************************************
# ** - RFD ~ TD: plot ----
#************************************
# hierarchical difference
fig_fd_td.hie <- data_outcome_td.hie %>% 
  #ggplot(aes(x=TD, y= log(fd), col=site, label=pair)) +
  ggplot(aes(x=abs(TD), y= abs(log(fd)), col=site, label=pair)) +
  geom_point() +
  #geom_text() +
  #geom_smooth(method="lm", se=FALSE) +
  scale_color_manual(values=c("Low site" = "#FB9A06FF", "Middle site" = "#6DCD59FF", "High site" = "#3E4A89FF")) +
  facet_wrap(~trait.id, scales = "free")
fig_fd_td.hie

#************************************
# ** - RFD ~ TD: lmer ----
#************************************
fd.lmer <- data_outcome_td.hie %>%
  split(.$trait.id) %>%
  purrr::map(.x = ., .f = lmer, formula = log(fd) ~ TD*site + is.novel + TD:is.novel  + (1|sps1) + (1|sps2))
  #purrr::map(.x = ., .f = lmer, formula = log(fd) ~ TD*site*is.novel  + (1|sps1) + (1|sps2))
  #purrr::map(.x = ., .f = lmer, formula = abs(log(fd)) ~ abs(TD)*site + is.novel + abs(TD):is.novel  + (1|sps1) + (1|sps2))

fd.lmer %>% purrr::map(Anova, type=2)
fd.lmer %>% purrr::map(summary)
#fd.lmer %>% purrr::map(rand)

# summary of the lmer results
fd.lmer.summary <- 
  fd.lmer %>% purrr::map_dfr(coef.lmer, .id="Trait")
fd.lmer.summary
#write.csv(fd.lmer.summary, "/Users/slyu1/Downloads/RFD.summary.lmer.csv")

#************************************
# ** - RFD ~ TD: slopes ----
#************************************
# make the table
fd.slope_site.novel <- data_outcome_td.hie %>% 
  group_by(trait.id, site, is.novel) %>%
  summarise(mean_cv(TD, na.rm=TRUE))
fd.slope_site.novel

# data
for(i in 1:nrow(fd.slope_site.novel)) {
  # i = 25
  di <- fd.slope_site.novel[i,]
  tr <- as.character(di$trait.id)
  st <- as.character(di$site)
  no <- as.character(di$is.novel)
  print(i)
  
  # get slope (hierarchical TD)
  slope.i <- lmer.slope(site = st, is.novel = no, 
                        #formula = "log(fd) ~ TD*site + is.novel + TD:is.novel + (1|sps1) + (1|sps2)",
                        formula = "log(fd) ~ TD*site + is.novel + TD:is.novel + (1|sps1) + (1|sps2)",
                        data = filter(data_outcome_td.hie, trait.id==tr))
  
  # get slope (absolute TD)
  #slope.i <- lmer.slope(site = st, is.novel = no, 
   #                     formula = "abs(log(fd)) ~ TD*site + is.novel + TD:is.novel + (1|sps1) + (1|sps2)",
    #                    data = filter(data_outcome_td.abs, trait.id==tr))
  
  # slopes and CI
  fd.slope_site.novel[i, "slope"] = slope.i$slope
  fd.slope_site.novel[i, "slope.min"] = slope.i$slope.min
  fd.slope_site.novel[i, "slope.max"] = slope.i$slope.max
  fd.slope_site.novel[i, "slope.p"] = slope.i$slope.p
  fd.slope_site.novel[i, "significant.slope"] = slope.i$significant.slope
}

fd.slope_site.novel
# write.csv(fd.slope_site.novel, "/Users/slyu1/LVSM/NICH/writing/CH3.0_traits-coexistence/Figures/RFD.slope.csv")

#************************************
# slopes ~ traits radar plot
# Organize the data to the wide format and plot
spec <- fd.slope_site.novel %>% 
  select(c(site, is.novel, slope, trait.id)) %>%
  build_wider_spec(names_from = trait.id, values_from = slope)

# radar plot
fig_fd.slope_radar <-
  fd.slope_site.novel %>%
  select(c(site, is.novel, slope, trait.id)) %>%
  pivot_wider_spec(spec) %>%
  mutate(group = paste(site, is.novel, sep="_")) %>%
  #filter(is.novel == "Current interactions") %>% group_by(site) %>%
  filter(is.novel == "Novel interactions") %>% group_by(site) %>%
  #filter(site == "Low site") %>% group_by(is.novel) %>%
  #filter(site == "Middle site") %>% group_by(is.novel) %>%
  #filter(site == "High site") %>% group_by(is.novel) %>%
  group_by(group) %>%
  select(3:20)  %>%
  ggradar(grid.min = -0.5, grid.mid = 0, grid.max = 0.5, values.radar = c("-0.5", "0", "0.5"),
          font.radar = "Arial",  background.circle.colour = "white",
          base.size=60, grid.label.size = 5, grid.line.width = 0.3, axis.label.size = 4.5, group.line.width = 1, group.point.size = 3.5,
          # group.colours = c("#FB9A06FF", "black"),
          group.colours = c("#3E4A89FF", "#FB9A06FF", "#6DCD59FF"),
          plot.legend=FALSE
          )
fig_fd.slope_radar

# across each site
# fig_fd.slope_radar_low <- fig_fd.slope_radar 
# fig_fd.slope_radar_middle <- fig_fd.slope_radar 
# fig_fd.slope_radar_high <- fig_fd.slope_radar 

# between current vs novel
# fig_fd.slope_radar_current <- fig_fd.slope_radar 
# fig_fd.slope_radar_novel <- fig_fd.slope_radar 

#*****************************************************************************
# 7. ND ~ trait difference ----
#*****************************************************************************

#************************************
# ** - ND ~ TD: plot ----
#************************************
# hierarchical difference
fig_nd_td.hie <- data_outcome_td.hie %>% 
  ggplot(aes(x=TD, y= log(1 - nd), col=site, label=pair)) +
  geom_point() +
  #geom_text() +
  #geom_smooth(method="lm", se=FALSE) +
  scale_color_manual(values=c("Low site" = "#FB9A06FF", "Middle site" = "#6DCD59FF", "High site" = "#3E4A89FF")) +
  facet_wrap(~trait.id, scales = "free")
fig_nd_td.hie

# absolute difference
fig_nd_td.abs <- data_outcome_td.abs %>% 
  ggplot(aes(x=TD, y= log(1 - nd), col=site, label=pair)) +
  geom_point() +
  #geom_text() +
  #geom_smooth(method="lm", se=FALSE) +
  scale_color_manual(values=c("Low site" = "#FB9A06FF", "Middle site" = "#6DCD59FF", "High site" = "#3E4A89FF")) +
  facet_wrap(~trait.id, scales = "free")
fig_nd_td.abs

#************************************
# ** - ND ~ TD: lmer ----
#************************************
nd.lmer <- 
  # hierarchical TD
  #data_outcome_td.hie %>%
  # absolute TD
  data_outcome_td.abs %>%
  split(.$trait.id) %>%
  purrr::map(.x = ., .f = lmer, formula = log(1-nd) ~ TD*site + is.novel + TD:is.novel  + (1|sps1) + (1|sps2))
  #purrr::map(.x = ., .f = lmer, formula = log(1-nd) ~ TD*site*is.novel  + (1|sps1) + (1|sps2))

nd.lmer %>% purrr::map(Anova)
nd.lmer %>% purrr::map(summary)
#fd.lmer %>% purrr::map(rand)

# summary of the lmer results
nd.lmer.summary <- 
  nd.lmer %>% purrr::map_dfr(coef.lmer, .id="Trait")
nd.lmer.summary
#write.csv(nd.lmer.summary, "/Users/slyu1/Downloads/NO.summary.lmer.csv")


#************************************
# ** - ND ~ TD: slopes ----
#************************************
# make the table
nd.slope_site.novel <- data_outcome_td.abs %>% 
  group_by(trait.id, site, is.novel) %>%
  summarise(mean_cv(abs(TD), na.rm=TRUE))
nd.slope_site.novel

# data
for(i in 1:nrow(nd.slope_site.novel)) {
  # i = 1
  di <- nd.slope_site.novel[i,]
  tr <- as.character(di$trait.id)
  st <- as.character(di$site)
  no <- as.character(di$is.novel)
  print(i)
  
  # get slope (hierarchical TD)
  slope.i <- lmer.slope(site = st, is.novel = no, 
                        formula = "log(1-nd) ~ TD*site + is.novel + TD:is.novel + (1|sps1) + (1|sps2)",
                        data = filter(data_outcome_td.abs, trait.id==tr))
  
  # slopes and CI
  nd.slope_site.novel[i, "slope"] = slope.i$slope
  nd.slope_site.novel[i, "slope.min"] = slope.i$slope.min
  nd.slope_site.novel[i, "slope.max"] = slope.i$slope.max
  nd.slope_site.novel[i, "slope.p"] = slope.i$slope.p
  nd.slope_site.novel[i, "significant.slope"] = slope.i$significant.slope
}
nd.slope_site.novel
# write.csv(nd.slope_site.novel, "/Users/slyu1/LVSM/NICH/writing/CH3.0_traits-coexistence/Figures/ND.slope.csv")

#************************************
# slopes ~ traits radar plot
# Organize the data to the wide format and plot
spec <- nd.slope_site.novel %>% 
  select(c(site, is.novel, slope, trait.id)) %>%
  build_wider_spec(names_from = trait.id, values_from = slope)

# radar plot
fig_nd.slope_radar <-
  nd.slope_site.novel %>%
  select(c(site, is.novel, slope, trait.id)) %>%
  mutate(slope = ifelse(slope > 0.5, 0.5, slope)) %>%
  mutate(slope = ifelse(slope < -0.5, -0.5, slope)) %>%
  pivot_wider_spec(spec) %>%
  mutate(group = paste(site, is.novel, sep="_")) %>%
  #filter(is.novel == "Current interactions") %>% group_by(site) %>%
  filter(is.novel == "Novel interactions") %>% group_by(site) %>%
  group_by(group) %>%
  select(3:20)  %>%
  ggradar(grid.min = -0.5, grid.mid = 0, grid.max = 0.5, values.radar = c("-0.5", "0", "0.5"),
          font.radar = "Arial", background.circle.colour = "white",
          base.size=60, grid.label.size = 5, grid.line.width = 0.3, axis.label.size = 4.5, group.line.width = 1, group.point.size = 3.5,
          #group.colours = c("#FB9A06FF", "black"), 
          group.colours = c("#3E4A89FF", "#FB9A06FF", "#6DCD59FF"),
          plot.legend=FALSE)
fig_nd.slope_radar

# current vs novel
# fig_nd.slope_radar_current <- fig_nd.slope_radar 
# fig_nd.slope_radar_novel <- fig_nd.slope_radar 

#*****************************************************************************
# 9. Coexistence ~ Multiple trait differences ----
#*****************************************************************************

#************************************
# ** - read data ----
#************************************
trait.scale
trait.ID <- colnames(trait.scale)[5:22]

# species
sps <- c("Anal", "Armo", "Asal", "Plal", "Poal", "Seca", "Trba",  # 7 alpine
         "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr")  # 7 lowland

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# population growth rates
pgr


#************************************
# ** - functions to generate all combinations ----
#************************************
combn.trait <- function(t, nn) {
  # t: trait pool
  # nn: number of traits for combination
  
  # number of traits
  n.trait <- length(t)
  
  # make a matrix
  trait.combn <- matrix(NA, ncol=n.trait)
  
  # combination for each number of traits
  for(i in nn) {
    # i = 5
    comb.i <- combn(1:n.trait, i, simplify = FALSE)
    trait.combn.i <- matrix(0, ncol=n.trait, nrow=length(comb.i))
    
    for(j in 1:nrow(trait.combn.i)) {
      # j = 1
      trait.combn.i[j, comb.i[[j]]] <- 1
    }
    # combine 
    trait.combn <- rbind(trait.combn, trait.combn.i)
  }
  
  trait.combn <- as.data.frame(trait.combn)
  colnames(trait.combn) <- t
  return(trait.combn)
}

# 
trait.combn <- combn.trait(trait.ID[1:15], 1:15)

# 32768 combinations
nrow(trait.combn)

# number of traits each combination
trait.number <- apply(trait.combn == 1 & !is.na(trait.combn), 1, sum)

#************************************
# ** - Calculate MTD for IGR ----
#************************************
# including pairs only used in previous analyses
pgr.nona <- 
  pgr %>% 
  filter(!(background.species %in% c("none", "site", "site_none", "species", "species_none"))) %>%
  filter(!is.na(pgr)) %>%
  filter(focal.species != background.species) %>% # remove intra pairs
  filter(focal.species != "Armo") %>%
  filter(!(focal.species == "Trba" & site == "Les Posses")) %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland_Highland", "Novel interactions", "Current interactions")) %>%
  mutate(is.novel = factor(is.novel, levels = c("Current interactions", "Novel interactions")) ) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low site", "Solalex" = "Middle site", "Anzeindaz" = "High site")) %>%
  mutate(site = factor(site, levels=c("Low site", "Middle site", "High site")))
  
# IGR at each site
pgr.low <- filter(pgr.nona, site == "Low site") # n = 75
pgr.middle <- filter(pgr.nona, site == "Middle site") # n = 83
pgr.high <- filter(pgr.nona, site == "High site") # n = n = 92

# traits at each site
trait.low <- data.frame(trait.scale[trait.scale$site == "Les Posses", ])
trait.middle <- data.frame(trait.scale[trait.scale$site == "Solalex", ])
trait.high <- data.frame(trait.scale[trait.scale$site == "Anzeindaz", ])

# species at each site
sps.low <- trait.low$species
sps.middle <- trait.middle$species
sps.high <- trait.high$species

# Trait combinations
set.seed(1)
id.trait.combn <- c(2:16, sample(17:nrow(trait.combn), 1623))

# new matrix to hold MTD
mtd.pgr <- matrix(NA, ncol=length(id.trait.combn), nrow=(nrow(pgr.low) + nrow(pgr.middle) + nrow(pgr.high)))

# calculate MTD for each combination
for(i in 1:length(id.trait.combn)) {
  # i = 1
  # trait
  trait.combn.i <- trait.combn[id.trait.combn[i], ]
  trait.i <- colnames(trait.combn.i)[trait.combn.i ==1] 
  
  # trait distance at low site
  dist.low <- as.matrix(dist(data.frame(trait.low[, trait.i])))
  colnames(dist.low) <- sps.low
  rownames(dist.low) <- sps.low
  td.low <- diag(dist.low[pgr.low$focal.species, pgr.low$background.species])
  
  # trait distance at middle site
  dist.middle <- as.matrix(dist(data.frame(trait.middle[, trait.i])))
  colnames(dist.middle) <- sps.middle
  rownames(dist.middle) <- sps.middle
  td.middle <- diag(dist.middle[pgr.middle$focal.species, pgr.middle$background.species])
  
  # trait distance at high site
  dist.high <- as.matrix(dist(data.frame(trait.high[, trait.i])))
  colnames(dist.high) <- sps.high
  rownames(dist.high) <- sps.high
  td.high <- diag(dist.high[pgr.high$focal.species, pgr.high$background.species])
  
  # calculate TD for each species pairs
  mtd.pgr[,i] <- c(td.low, td.middle, td.high)
}
mtd.pgr <- as.data.frame(mtd.pgr)
colnames(mtd.pgr) <- paste("Trait.combn", id.trait.combn, sep="_")
mtd.pgr$Trait.combn_3

mean_cv(mtd.pgr$Trait.combn_325)
mean_cv(mtd.pgr$Trait.combn_30430)

#************************************
# ** - IGR ~ MTD  ----
#************************************
# combine outcome and MTD
data_pgr_mtd <-  
  bind_cols(bind_rows(pgr.low, pgr.middle, pgr.high), mtd.pgr) %>%
  pivot_longer(cols = colnames(mtd.pgr), values_to = "MTD", names_to = "Trait.combination")

#************************************
# low site
models_pgr.MTD_low <- data_pgr_mtd %>%
  filter(site == "Low site") %>%
  split(.$Trait.combination) %>%
  purrr::map(.x=., .f=lmer, formula = log(pgr) ~ MTD + (1|focal.species) + (1|background.species))

# model performance
model.performance_pgr.MTD_low <- models_pgr.MTD_low %>%
  purrr::map_dfr(.f = extract.lmer, .id = "Trait.combination") %>%
  separate(col = Trait.combination, sep = "_", into = c("Trait.combination", "ID.trait.combn"))

# Number of traits for each combination
model.performance_pgr.MTD_low$trait.number <- trait.number[as.numeric(model.performance_pgr.MTD_low$ID.trait.combn)]

# delta AICc
model.performance_pgr.MTD_low <- model.performance_pgr.MTD_low %>%
  mutate(delta.AICc = AICc - min(AICc)) %>%
  mutate(site = "Low site")

# AIC ~ trait number plot
fig_pgr.MTD_AICc_low <- model.performance_pgr.MTD_low %>%
  ggplot(aes(x=factor(trait.number), y = delta.AICc)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)
fig_pgr.MTD_AICc_low

# AIC ~ trait number test
model.performance_pgr.MTD_low %>%
  lm(delta.AICc ~ trait.number, data=.) %>%
  summary()
#Anova()

# trait combination with low AICc
model.performance_pgr.MTD_low %>%
  arrange(AICc)

# R2 marginal ~ trait number
model.performance_pgr.MTD_low %>%
  ggplot(aes(x=factor(trait.number), y = `R2 marginal`)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# trait combination with highest R2
model.performance_pgr.MTD_low %>%
  arrange(desc(`R2 marginal`))

# slope ~ trait number
model.performance_pgr.MTD_low %>%
  ggplot(aes(x=factor(trait.number), y = Slope)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# P value ~ trait number
model.performance_pgr.MTD_low %>%
  ggplot(aes(x=factor(trait.number), y = log(`P value`))) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# correlations between AICc, R2, slope and P value
model.performance_pgr.MTD_low %>%
  select(AICc, Slope, `P value`, `R2 marginal`) %>%
  mutate(`P value` = log(`P value`)) %>%
  pairs()

# plot
data_pgr_mtd %>%
  filter(site == "Low site") %>%
  filter(Trait.combination == "Trait.combn_9198") %>%
  ggplot(aes(x=MTD, y = log(pgr))) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE)

#************************************
# Middle site
models_pgr.MTD_middle <- data_pgr_mtd %>%
  filter(site == "Middle site") %>%
  split(.$Trait.combination) %>%
  purrr::map(.x=., .f=lmer, formula = log(pgr) ~ MTD + (1|focal.species) + (1|background.species))
#purrr::map(.x=., .f=lmer, formula = log(pgr) ~ MTD*is.novel + (1|focal.species) + (1|background.species))

# model performance
model.performance_pgr.MTD_middle <- models_pgr.MTD_middle %>%
  purrr::map_dfr(.f = extract.lmer, .id = "Trait.combination") %>%
  separate(col = Trait.combination, sep = "_", into = c("Trait.combination", "ID.trait.combn"))

# Number of traits for each combination
model.performance_pgr.MTD_middle$trait.number <- trait.number[as.numeric(model.performance_pgr.MTD_middle$ID.trait.combn)]

# delta AICc
model.performance_pgr.MTD_middle <- 
  model.performance_pgr.MTD_middle %>%
  mutate(delta.AICc = AICc - min(AICc)) %>%
  mutate(site = "Middle site")

# AIC ~ trait number
fig_pgr.MTD_AICc_middle <- 
  model.performance_pgr.MTD_middle %>%
  ggplot(aes(x=factor(trait.number), y = delta.AICc)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)
fig_pgr.MTD_AICc_middle

# AIC ~ trait number
model.performance_pgr.MTD_middle %>%
  lm(delta.AICc ~ trait.number, data=.) %>%
  summary()
  #Anova()

# trait combination with low AICc
model.performance_pgr.MTD_middle %>%
  arrange(AICc)

# slope ~ trait number
model.performance_pgr.MTD_middle %>%
  ggplot(aes(x=factor(trait.number), y = Slope)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

#************************************
# High site
models_pgr.MTD_high <- data_pgr_mtd %>%
  filter(site == "High site") %>%
  split(.$Trait.combination) %>%
  purrr::map(.x=., .f=lmer, formula = log(pgr) ~ MTD + (1|focal.species) + (1|background.species))
#purrr::map(.x=., .f=lmer, formula = log(pgr) ~ MTD*is.novel + (1|focal.species) + (1|background.species))

# model performance
model.performance_pgr.MTD_high <- models_pgr.MTD_high %>%
  purrr::map_dfr(.f = extract.lmer, .id = "Trait.combination") %>%
  separate(col = Trait.combination, sep = "_", into = c("Trait.combination", "ID.trait.combn"))

# Number of traits for each combination
model.performance_pgr.MTD_high$trait.number <- trait.number[as.numeric(model.performance_pgr.MTD_high$ID.trait.combn)]

# delta AICc
model.performance_pgr.MTD_high <- 
  model.performance_pgr.MTD_high %>%
  mutate(delta.AICc = AICc - min(AICc)) %>%
  mutate(site = "High site")
  
# AIC ~ trait number plot
fig_pgr.MTD_AICc_high <- 
  model.performance_pgr.MTD_high %>%
  ggplot(aes(x=factor(trait.number), y = delta.AICc)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)
fig_pgr.MTD_AICc_high

# AIC ~ trait number test
model.performance_pgr.MTD_high %>%
  lm(delta.AICc ~ trait.number, data=.) %>%
  summary()
  #Anova()

# trait combination with low AICc
model.performance_pgr.MTD_high %>%
  arrange(AICc)

# slope ~ trait number
model.performance_pgr.MTD_high %>%
  ggplot(aes(x=factor(trait.number), y = Slope)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

#************************************
# combine the three sites
# AICc
fig_pgr.MTD_AICc <- bind_rows(model.performance_pgr.MTD_low, 
          model.performance_pgr.MTD_middle,
          model.performance_pgr.MTD_high) %>%
  mutate(site = factor(site, levels=c("Low site", "Middle site", "High site"))) %>%
  ggplot(aes(x=trait.number, y = delta.AICc, group = factor(trait.number)),) +
  geom_jitter(width=0.15) +
  stat_summary(fun.data = mean_sd.ggplot, aes(col=site)) + 
  facet_wrap(~ site) +
  scale_color_manual(values=c("Low site" = "#FB9A06FF", "Middle site" = "#6DCD59FF", "High site" = "#3E4A89FF")) +
  scale_x_continuous(name = "Number of traits", n.breaks = 10) +
  scale_y_continuous(name = "delta AICc")
fig_pgr.MTD_AICc

# slope
fig_pgr.MTD_slope <- 
  bind_rows(model.performance_pgr.MTD_low, 
          model.performance_pgr.MTD_middle,
          model.performance_pgr.MTD_high) %>%
  mutate(is.p.significant = ifelse(`P value` < 0.05, "*", "")) %>%
  mutate(site = factor(site, levels=c("Low site", "Middle site", "High site"))) %>%
  ggplot(aes(x=trait.number, y = Slope, group = factor(trait.number)),) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(width=0.15, aes(alpha=is.p.significant)) +
  stat_summary(fun.data = mean_sd.ggplot, aes(col=site)) + 
  facet_wrap(~ site) +
  scale_color_manual(values=c("Low site" = "#FB9A06FF", "Middle site" = "#6DCD59FF", "High site" = "#3E4A89FF")) +
  scale_x_continuous(name = "Number of traits", n.breaks = 10) +
  scale_y_continuous(name = "Slope")
fig_pgr.MTD_slope
fig_pgr.MTD_AICc


#************************************
# ** - Calculate MTD for outcome  ----
#************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# competitive outcomes
outcome

# exclude pairs without predicted outcomes
outcome.nona <- outcome[!is.na(outcome$nd) & !is.na(outcome$fd), ]

# outcomes at each sites
outcome.low <- outcome.nona[outcome.nona$site == "Les Posses",] # n = 32
outcome.middle <- outcome.nona[outcome.nona$site == "Solalex",] # n = 34
outcome.high <- outcome.nona[outcome.nona$site == "Anzeindaz",] # n = 41

# traits at each site
trait.low <- data.frame(trait.scale[trait.scale$site == "Les Posses", ])
trait.middle <- data.frame(trait.scale[trait.scale$site == "Solalex", ])
trait.high <- data.frame(trait.scale[trait.scale$site == "Anzeindaz", ])

# species at each site
sps.low <- trait.low$species
sps.middle <- trait.middle$species
sps.high <- trait.high$species

# Trait combinations
set.seed(1)
id.trait.combn <- c(2:16, sample(17:nrow(trait.combn), 1623))

# new matrix to hold MTD
mtd.outcome <- matrix(NA, ncol=length(id.trait.combn), nrow=(nrow(outcome.low) + nrow(outcome.middle) + nrow(outcome.high)))

# calculate MTD for each combination
for(i in 1:length(id.trait.combn)) {
  # i = 1612
  # trait
  trait.combn.i <- trait.combn[id.trait.combn[i], ]
  trait.i <- colnames(trait.combn.i)[trait.combn.i ==1] 
  
  # trait distance at low site
  dist.low <- as.matrix(dist(data.frame(trait.low[, trait.i])))
  colnames(dist.low) <- sps.low
  rownames(dist.low) <- sps.low
  td.low <- diag(dist.low[outcome.low$sps1, outcome.low$sps2])
  
  # trait distance at middle site
  dist.middle <- as.matrix(dist(data.frame(trait.middle[, trait.i])))
  colnames(dist.middle) <- sps.middle
  rownames(dist.middle) <- sps.middle
  td.middle <- diag(dist.middle[outcome.middle$sps1, outcome.middle$sps2])
  
  # trait distance at high site
  dist.high <- as.matrix(dist(data.frame(trait.high[, trait.i])))
  colnames(dist.high) <- sps.high
  rownames(dist.high) <- sps.high
  td.high <- diag(dist.high[outcome.high$sps1, outcome.high$sps2])
  
  # calculate TD for each species pairs
  mtd.outcome[,i] <- c(td.low, td.middle, td.high)
}
mtd.outcome <- as.data.frame(mtd.outcome)
colnames(mtd.outcome) <- paste("Trait.combn", id.trait.combn, sep="_")
mtd.outcome[,500]

#************************************
# ** - ND ~ MTD  ----
#************************************
# combine outcome and MTD
data_outcome_mtd <-  
  bind_cols(bind_rows(outcome.low, outcome.middle, outcome.high), mtd.outcome) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low site", "Solalex" = "Middle site", "Anzeindaz"= "High site")) %>%
  mutate(site = factor(site, levels = c("Low site", "Middle site", "High site"))) %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel interactions", "Current interactions")) %>%
  mutate(is.novel = factor(is.novel, levels = c("Novel interactions", "Current interactions"))) %>%
  pivot_longer(cols = colnames(mtd.outcome), values_to = "MTD", names_to = "Trait.combination")
  
#************************************
# low site
models_ND.MTD_low <- data_outcome_mtd %>%
  filter(site == "Low site") %>%
  split(.$Trait.combination) %>%
  purrr::map(.x=., .f=lm, formula = log(1-nd) ~ MTD)
  #purrr::map(.x=., .f=lmer, formula = log(1-nd) ~ MTD + (1|sps1) + (1|sps2))

# how to deal with singularity?!
summary(models_ND.MTD_low$Trait.combn_10)

# model performance
model.performance_ND.MTD_low <- models_ND.MTD_low %>%
  purrr::map_dfr(.f = extract.lm, .id = "Trait.combn") %>%
  #purrr::map_dfr(.f = extract.lmer, .id = "Trait.combn") %>%
  separate(col = Trait.combn, sep = "_", into = c("Trait.combn", "ID.trait.combn"))

# trait number
model.performance_ND.MTD_low$trait.number <- trait.number[as.numeric(model.performance_ND.MTD_low$ID.trait.combn)]

# delta AICc
model.performance_ND.MTD_low <- 
  model.performance_ND.MTD_low %>%
  mutate(delta.AICc = AICc - min(AICc)) %>%
  mutate(site = "Low site")

# AIC ~ trait number plot
model.performance_ND.MTD_low %>%
  ggplot(aes(x=factor(trait.number), y = delta.AICc)) +
  geom_jitter(width=0.15) +
  stat_summary(fun.data = mean_sd.ggplot, col="#FB9A06FF")

# AIC ~ trait number test
model.performance_ND.MTD_low %>%
  lm(delta.AICc ~ trait.number, data=.) %>%
  summary()
#Anova()

# Slope ~ trait number plot
model.performance_ND.MTD_low %>%
  ggplot(aes(x=factor(trait.number), y = Slope)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# P value ~ trait number
model.performance_ND.MTD_low %>%
  ggplot(aes(x=factor(trait.number), y = `P value`)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# R2 marginal ~ trait number
model.performance_ND.MTD_low %>%
  ggplot(aes(x=factor(trait.number), y = `R2`)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)


#************************************
# middle site
models_ND.MTD_middle <- data_outcome_mtd %>%
  filter(site == "Middle site") %>%
  split(.$Trait.combination) %>%
  purrr::map(.x=., .f=lm, formula = log(1-nd) ~ MTD)
  #purrr::map(.x=., .f=lmer, formula = log(1-nd) ~ MTD + (1|sps1) + (1|sps2))

summary(models_ND.MTD_middle$Trait.combn_10)
# how to deal with singularity?!

# model performance
model.performance_ND.MTD_middle <- models_ND.MTD_middle %>%
  purrr::map_dfr(.f = extract.lm, .id = "Trait.combn") %>%
  separate(col = Trait.combn, sep = "_", into = c("Trait.combn", "ID.trait.combn"))

# trait number
model.performance_ND.MTD_middle$trait.number <- trait.number[as.numeric(model.performance_ND.MTD_middle$ID.trait.combn)]

# delta AICc
model.performance_ND.MTD_middle <- 
  model.performance_ND.MTD_middle %>%
  mutate(delta.AICc = AICc - min(AICc)) %>%
  mutate(site = "Middle site")

# AIC ~ trait number plot
model.performance_ND.MTD_middle %>%
  #filter(delta.AICc > 2) %>%
  ggplot(aes(x=factor(trait.number), y = delta.AICc)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# AIC ~ trait number test
model.performance_ND.MTD_middle %>%
  filter(delta.AICc > 2) %>%
  lm(delta.AICc ~ trait.number, data=.) %>%
  summary()
#Anova()

# Slope ~ trait number plot
model.performance_ND.MTD_middle %>%
  ggplot(aes(x=factor(trait.number), y = Slope)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# R2 marginal ~ trait number
model.performance_ND.MTD_middle %>%
  ggplot(aes(x=factor(trait.number), y = `R2`)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

#************************************
# High site
models_ND.MTD_high <- data_outcome_mtd %>%
  filter(site == "High site") %>%
  split(.$Trait.combination) %>%
  purrr::map(.x=., .f=lm, formula = log(1-nd) ~ MTD)
  #purrr::map(.x=., .f=lmer, formula = log(1-nd) ~ MTD + (1|sps1) + (1|sps2))

summary(models_ND.MTD_high$Trait.combn_10)
# how to deal with singularity?!

# model performance
model.performance_ND.MTD_high <- models_ND.MTD_high %>%
  purrr::map_dfr(.f = extract.lm, .id = "Trait.combn") %>%
  separate(col = Trait.combn, sep = "_", into = c("Trait.combn", "ID.trait.combn"))

# trait number
model.performance_ND.MTD_high$trait.number <- trait.number[as.numeric(model.performance_ND.MTD_high$ID.trait.combn)]

# delta AICc
model.performance_ND.MTD_high <- 
  model.performance_ND.MTD_high %>%
  mutate(delta.AICc = AICc - min(AICc)) %>%
  mutate(site = "High site")

# AIC ~ trait number plot
model.performance_ND.MTD_high %>%
  ggplot(aes(x=factor(trait.number), y = delta.AICc)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# AIC ~ trait number test
model.performance_ND.MTD_high %>%
  lm(delta.AICc ~ trait.number, data=.) %>%
  summary()
#Anova()

# Slope ~ trait number plot
model.performance_ND.MTD_high %>%
  ggplot(aes(x=factor(trait.number), y = Slope)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# P value ~ trait number
model.performance_ND.MTD_high %>%
  ggplot(aes(x=factor(trait.number), y = `P value`)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# R2 marginal ~ trait number
model.performance_ND.MTD_high %>%
  ggplot(aes(x=factor(trait.number), y = `R2`)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

#************************************
# ** - RFD ~ MTD  ----
#************************************
# combine outcome and MTD
data_outcome_mtd <-  
  bind_cols(bind_rows(outcome.low, outcome.middle, outcome.high), mtd.outcome) %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low site", "Solalex" = "Middle site", "Anzeindaz"= "High site")) %>%
  mutate(site = factor(site, levels = c("Low site", "Middle site", "High site"))) %>%
  mutate(is.novel = ifelse(origin.pair == "Lowland-highland", "Novel interactions", "Current interactions")) %>%
  mutate(is.novel = factor(is.novel, levels = c("Novel interactions", "Current interactions"))) %>%
  pivot_longer(cols = colnames(mtd.outcome), values_to = "MTD", names_to = "Trait.combination")

#************************************
# low site
models_RFD.MTD_low <- data_outcome_mtd %>%
  filter(site == "Low site") %>%
  split(.$Trait.combination) %>%
  purrr::map(.x=., .f=lmer, formula = abs(log(fd)) ~ MTD +  + (1|sps1) + (1|sps2))
  #purrr::map(.x=., .f=lmer, formula = abs(log(fd)) ~ MTD*is.novel +  + (1|sps1) + (1|sps2))

# how to deal with singularity?!
summary(models_RFD.MTD_low$Trait.combn_10)

# model performance
model.performance_RFD.MTD_low <- models_RFD.MTD_low %>%
  purrr::map_dfr(.f = extract.lmer, .id = "Trait.combn") %>%
  separate(col = Trait.combn, sep = "_", into = c("Trait.combn", "ID.trait.combn"))

# trait number
model.performance_RFD.MTD_low$trait.number <- trait.number[as.numeric(model.performance_RFD.MTD_low$ID.trait.combn)]

# delta AICc
model.performance_RFD.MTD_low <- 
  model.performance_RFD.MTD_low %>%
  mutate(delta.AICc = AICc - min(AICc)) %>%
  mutate(site = "Low site")

# AIC ~ trait number plot
model.performance_RFD.MTD_low %>%
  ggplot(aes(x=factor(trait.number), y = delta.AICc)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# AIC ~ trait number test
model.performance_RFD.MTD_low %>%
  lm(delta.AICc ~ trait.number, data=.) %>%
  summary()
#Anova()

# Slope ~ trait number plot
model.performance_RFD.MTD_low %>%
  ggplot(aes(x=factor(trait.number), y = Slope)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# P value ~ trait number
model.performance_RFD.MTD_low %>%
  ggplot(aes(x=factor(trait.number), y = `P value`)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# R2 marginal ~ trait number
model.performance_RFD.MTD_low %>%
  ggplot(aes(x=factor(trait.number), y = `R2 marginal`)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# correlation between AICc, slope, P value, 
model.performance_RFD.MTD_low %>%
  select(AICc, Slope, `P value`, `R2 marginal`) %>%
  pairs()

#************************************
# middle site
models_RFD.MTD_middle <- data_outcome_mtd %>%
  filter(site == "Middle site") %>%
  split(.$Trait.combination) %>%
  purrr::map(.x=., .f=lmer, formula = abs(log(fd)) ~ MTD +  + (1|sps1) + (1|sps2))

summary(models_RFD.MTD_middle$Trait.combn_10)
# how to deal with singularity?!

# model performance
model.performance_RFD.MTD_middle <- models_RFD.MTD_middle %>%
  purrr::map_dfr(.f = extract.lmer, .id = "Trait.combn") %>%
  separate(col = Trait.combn, sep = "_", into = c("Trait.combn", "ID.trait.combn"))

# trait number
model.performance_RFD.MTD_middle$trait.number <- trait.number[as.numeric(model.performance_RFD.MTD_middle$ID.trait.combn)]

# delta AICc
model.performance_RFD.MTD_middle <- 
  model.performance_RFD.MTD_middle %>%
  mutate(delta.AICc = AICc - min(AICc)) %>%
  mutate(site = "Middle site")

# AIC ~ trait number plot
model.performance_RFD.MTD_middle %>%
  ggplot(aes(x=factor(trait.number), y = delta.AICc)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# AIC ~ trait number test
model.performance_RFD.MTD_middle %>%
  lm(delta.AICc ~ trait.number, data=.) %>%
  summary()
#Anova()

# Slope ~ trait number plot
model.performance_RFD.MTD_middle %>%
  ggplot(aes(x=factor(trait.number), y = Slope)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# R2 marginal ~ trait number
model.performance_RFD.MTD_middle %>%
  ggplot(aes(x=factor(trait.number), y = `R2 marginal`)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

#************************************
# High site
models_RFD.MTD_high <- data_outcome_mtd %>%
  filter(site == "High site") %>%
  split(.$Trait.combination) %>%
  purrr::map(.x=., .f=lmer, formula = abs(log(fd)) ~ MTD +  + (1|sps1) + (1|sps2))

summary(models_RFD.MTD_high$Trait.combn_10)
# how to deal with singularity?!

# model performance
model.performance_RFD.MTD_high <- models_RFD.MTD_high %>%
  purrr::map_dfr(.f = extract.lmer, .id = "Trait.combn") %>%
  separate(col = Trait.combn, sep = "_", into = c("Trait.combn", "ID.trait.combn"))

# trait number
model.performance_RFD.MTD_high$trait.number <- trait.number[as.numeric(model.performance_RFD.MTD_high$ID.trait.combn)]

# delta AICc
model.performance_RFD.MTD_high <- 
  model.performance_RFD.MTD_high %>%
  mutate(delta.AICc = AICc - min(AICc)) %>%
  mutate(site = "High site")

# AIC ~ trait number plot
model.performance_RFD.MTD_high %>%
  ggplot(aes(x=factor(trait.number), y = delta.AICc)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# AIC ~ trait number test
model.performance_RFD.MTD_high %>%
  lm(delta.AICc ~ trait.number, data=.) %>%
  summary()
#Anova()

# Slope ~ trait number plot
model.performance_RFD.MTD_high %>%
  ggplot(aes(x=factor(trait.number), y = Slope)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# P value ~ trait number
model.performance_RFD.MTD_high %>%
  ggplot(aes(x=factor(trait.number), y = `P value`)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

# R2 marginal ~ trait number
model.performance_RFD.MTD_high %>%
  ggplot(aes(x=factor(trait.number), y = `R2 marginal`)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.6)

#*****************************************************************************
# Figures ----
#*****************************************************************************

#************************************
# ** - Figure 1 ----
#************************************
p <- fig.coexistence_site +
  theme_bw(base_family = "Arial") + 
  theme(
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.title=element_text(size=12),
    axis.line = element_line(colour = "black"),
    panel.grid = element_blank(),
    
    legend.title = element_blank(),
    legend.position = c(0.88, 0.82),
    legend.text = element_text(size=8),
    
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size=12)) 
p
ggsave(p, filename = "/Users/slyu1/Downloads/fig_coexistence_site.pdf", device = cairo_pdf, 
       width = 7.6, height = 3, units = "in")

#************************************
# ** - Figure 2 ----
#************************************

# traits
p <- fig_igr.slope_trait +
  theme_bw(base_family = "Arial") + 
  theme(
    axis.text.x=element_text(size=8),
    axis.text.y=element_text(size=11, colour="black"),
    axis.title=element_text(size=11),
    axis.line = element_line(colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "right", 
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size=10))
p
ggsave(p, filename = "/Users/slyu1/Downloads/fig.1.1.pdf", device = cairo_pdf, 
       width = 8.2, height = 5.8, units = "in")

#************************************
# ** - Figure 3 ----
#************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RFD slopes of current vs novel pairs
# Current pairs
p <- fig_fd.slope_radar_current +
  theme_void(base_family = "Arial") + 
  theme(legend.position = "null") 
p
ggsave(p, filename = "/Users/slyu1/Downloads/figure_FD.slope_current.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")

# Novel pairs
p <- fig_fd.slope_radar_novel +
  theme_void(base_family = "Arial") + 
  theme(legend.position = "null") 
p
ggsave(p, filename = "/Users/slyu1/Downloads/figure_FD.slope_novel.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ND slopes between current vs novel pairs
# current pairs
p <- fig_nd.slope_radar_current +
  theme_void(base_family = "Arial") + 
  theme(legend.position = "null") 
p
ggsave(p, filename = "/Users/slyu1/Downloads/figure_ND.slope_current_scaled.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")

# current pairs
p <- fig_nd.slope_radar_novel +
  theme_void(base_family = "Arial") + 
  theme(legend.position = "null") 
p
ggsave(p, filename = "/Users/slyu1/Downloads/figure_ND.slope_novel_scaled.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")
