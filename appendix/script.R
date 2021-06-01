#####################################################################################
## R code for "The Effect of Legislature Size on Public Spending: A Meta-Analysis" ##
## Huzeyfe Alptekin, Danilo Freire, Umberto Mignozzetti, and Catarina Roman        ##
## June 2021                                                                       ##
#####################################################################################

## Install and load required packages
set.seed(732578) # From random.org

# Required packages
pkgs <- c("tidyverse", "meta", "metafor",
          "readxl", "devtools", "data.table",
          "knitr", "gridGraphics", "gridExtra",
          "ggpubr", "kableExtra", "magick",
          'stargazer', 'pander', 'broom')

# Install if not already installed
install <- function(x) {
  if (x %in% rownames(installed.packages()) == FALSE)
    install.packages(x, dependencies = T,
                     repos = "http://cran.us.r-project.org")
}
lapply(pkgs, install)
devtools::install_github("isubirana/compareGroups")
devtools::install_github("MathiasHarrer/dmetar")

# Load packages
lapply(pkgs, require, character.only = T)
library("compareGroups"); library(dmetar)


## Auxiliary functions to estimate meta-analysis models
## and create plots

# Broom the multilevel model
broom_mod <- function(mod, subgroup = FALSE) {
  aux <- data.frame(predict(mod))
  if (subgroup) {
    teeff = aux$pred
    seteeff = aux$se
    labeff = c("Subgroup Effect")
    loeff = c(aux$ci.lb)
    upeff = c(aux$ci.ub)
  } else {
    teeff = c(aux$pred, NA)
    seteeff = c(aux$se, NA)
    labeff = c("Overall Effect",
               "Prediction Interval")
    loeff = c(aux$ci.lb, aux$cr.lb)
    upeff = c(aux$ci.ub, aux$cr.ub)
  }
  mod2 <- tibble(
    TE = as.numeric(mod$yi),
    seTE = sqrt(mod$vi),
    studlab = mod$slab,
    lower = as.numeric(mod$yi)-1.96*sqrt(mod$vi),
    upper = as.numeric(mod$yi)+1.96*sqrt(mod$vi),
    group = 'A') %>%
    bind_rows(.,
              aux = tibble(
                TE = teeff,
                seTE = seteeff,
                studlab = labeff,
                lower = loeff,
                upper = upeff,
                group = "B")) %>%
    group_by(studlab) %>%
    mutate(studlab2 = paste0(studlab, "_", 1:n())) %>%
    ungroup()
  return(mod2)
}

# Estimation of heterogeneous effects 
estim_het <- function(dat, yi, v, random, slab, hetvar = NULL) {
  if (is.null(hetvar)) {
    stop('Hetvar should be different than null.')
  }
  levshetvar <- unique(hetvar)
  fullmod <- rma.mv(yi = yi, V = v, random = random,
                    data = dat, slab = slab, method = 'REML',
                    test = 't', tdist = T)
  res <- tibble()
  for (i in levshetvar) {
    partmod <- rma.mv(yi = yi, V = v, 
                      random = random,
                      data = dat, slab = slab, 
                      method = 'REML',
                      test = 't', tdist = T, 
                      subset = hetvar == i)
    partmod <- broom_mod(partmod, subgroup = T)
    partmod$byvar = i
    partmod <- bind_rows(
      tibble(TE = NA, seTE = NA,
             studlab = toupper(i), lower = NA,
             upper = NA, group = "B", byvar = i),
      partmod
    )
    res <- bind_rows(res, partmod)
  }
  aux <- data.frame(predict(fullmod))
  res <- bind_rows(res,
                   tibble(
                     byvar = NA,
                     TE = c(aux$pred, NA),
                     seTE = c(aux$se, NA),
                     studlab = c("Overall Effect", "Prediction Interval"),
                     lower = c(aux$ci.lb, aux$cr.lb),
                     upper = c(aux$ci.ub, aux$cr.ub),
                     group = "B"))
  res <- data.frame(res)
  res$byvar <- toupper(res$byvar)
  return(res)
}

# Build plot function for forest plots
build_forest <- function(mod, capt, lsize = 22, ttl = NULL) {
  if(class(mod)[1] == "rma.mv") {
    mod2 <- broom_mod(mod)
  } else {
  # Build dataset for plot
  mod2 <- tibble(
    TE = mod$TE,
    seTE = mod$seTE,
    studlab = mod$studlab,
    lower = mod$lower,
    upper = mod$upper,
    group = "A") %>%
    bind_rows(.,
              aux = tibble(
                TE = c(mod$TE.random, NA),
                seTE = c(mod$seTE.random, NA),
                studlab = c("Overall Effect",
                            "Prediction Interval"),
                lower = c(mod$lower.random,
                          mod$lower.predict),
                upper = c(mod$upper.random,
                          mod$upper.predict),
                group = "B")) %>%
    group_by(studlab) %>%
    mutate(studlab2 = paste0(studlab, "_", 1:n())) %>%
    ungroup()
  }
  # Graph limits
  limg <- max(abs(c(mod2$lower, mod2$upper)))
  # Build plot
  p <- mod2 %>%
    ggplot(aes(y = reorder(studlab2, TE),
               x = TE, xmin = lower, xmax = upper)) +
    geom_point(aes(color = group)) +
    geom_errorbarh(aes(color = group),
                   height = 0.1) +
    scale_color_manual(values = c("#000000", "#8b0000")) +
    scale_x_continuous(limits = c(-1.1 * limg, 1.1 * limg)) +
    scale_y_discrete(
      labels = function(x)
        str_replace(x, "_[0-9]*$", "")) +
    geom_vline(xintercept = 0,
               color = "#000000", linetype = "dashed") +
    labs(x = "",
         y = "") +
    facet_grid(group~., scales = "free", space = "free") +
    labs(caption = capt,
         title = ttl) +
    theme_minimal() %+replace%
    theme(strip.text.y = element_blank(),
          legend.position = "none",
          axis.text.y = element_text(size = .8 * lsize,
                                     hjust = 1),
          axis.text.x = element_text(size = .6 * lsize,
                                     hjust = 1.1),
          plot.caption = element_text(size = lsize),
          plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    margin = margin(0, 0, 10, 0)),
          panel.grid.major = element_blank())
  return(p)
}

# Build forest plot for heterogeneous analysis
build_forest_het <- function(dat, coef, v, slab, capt, lsize = 22, ttl = NULL, hetvar = NULL) {
  mod2 <- estim_het(dat = dat, 
                    yi = coef, 
                    v = v, 
                    random = ~ 1 | id_level1/id_level2, 
                    slab = slab, 
                    hetvar = hetvar)
  mod2 <- data.frame(mod2)
  mod2$byvar <- toupper(mod2$byvar)
  TEaux <- mod2$TE
  TEaux[mod2$studlab== 'Subgroup Effect'] = TEaux[mod2$studlab== 'Subgroup Effect'] - 100
  # Graph limits
  limg <- max(abs(c(mod2$lower, mod2$upper)))
  # Build plot
  p <- mod2 %>%
    ggplot(aes(y = reorder(studlab, TEaux),
               x = TE, xmin = lower, xmax = upper)) +
    geom_point(aes(color = group)) +
    geom_errorbarh(aes(color = group),
                   height = 0.1) +
    scale_color_manual(values = c("#000000", "#8b0000")) +
    scale_x_continuous(limits = c(-1.1 * limg, 1.1 * limg)) +
    scale_y_discrete(
      labels = function(x)
        str_replace(x, "_[0-9]*$", "")) +
    geom_vline(xintercept = 0,
               color = "#000000", linetype = "dashed") +
    labs(x = "",
         y = "") +
    facet_grid(byvar~., scales = "free", space = "free") +
    labs(caption = capt,
         title = ttl) +
    theme_minimal() %+replace%
    theme(strip.text.y = element_blank(),
          legend.position = "none",
          axis.text.y = element_text(size = .8 * lsize,
                                     hjust = 1),
          axis.text.x = element_text(size = .6 * lsize,
                                     hjust = 1.1),
          plot.caption = element_text(size = lsize),
          plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    margin = margin(0, 0, 10, 0)),
          panel.grid.major = element_blank())
  return(p)
}

## Load datasets
load("../dataset/dataCoefs.RData")

## Code for figures included in Section A of the Supplementary Material
xstar <- function(x, beta, theta, alpha = 0.7, s = 0.5, m = 10) {
  n = x
  a = (alpha/theta)^(1/(theta-alpha))
  b = (n*m)^((beta+theta-1)/(theta-alpha))
  c = ((n^(2-beta-theta))/(n-n*s+s))^(1/(theta-alpha))
  return (a*b*c)
}

curve(xstar(x, beta = 0.35, theta = 0.8), from = 10, to = 20,
      main = 'Increasing project size when varying legislature size \n (high deadweight losses)', xlab = 'Project Size', ylab = 'Number of Districts (legislature size)')
curve(xstar(x, beta = 0.35, theta = 0.5), from = 10, to = 20,
      main = 'Decreasing project size when varying legislature size \n (low deadweight losses)', xlab = 'Project Size', ylab = 'Number of Districts (legislature size)')

## Descriptive Statistics

# Study year
dat %>%
  select(id, year) %>%
  unique() %>%
  ggplot(aes(x = as.factor(year))) +
    geom_bar(color = "black") +
  labs(x = "",
       y = "") +
  theme_bw()

# Frequency of published papers
dat %>%
  select(id, published) %>%
  unique() %>%
  ggplot(aes(x = as.factor(published))) +
    geom_bar(color = "black") +
  labs(x = "Published Study?",
       y = "") +
  theme_bw()

# Electoral System
dat %>%
  select(id, elecsys2) %>%
  unique() %>%
  ggplot(aes(x = as.factor(elecsys2))) +
    geom_bar(color = "black") +
  labs(x = "Electoral Systems",
       y = "") +
  theme_bw()

# Dependent Variables
dat %>%
  select(id, depvar2) %>%
  unique() %>%
  mutate(depvar2 = factor(depvar2,
                          labels = c("Expenditure Per Capita",
                                     "Expenditure as Percentage GDP",
                                     "Log Expenditure Per Capita"))) %>%
  ggplot(aes(x = depvar2)) +
    geom_bar(color = "black") +
  labs(x = "",
       y = "") +
  coord_flip() +
  theme_bw()

# Independent Variables
dat %>%
  select(id, indepvar2) %>%
  unique() %>%
  mutate(indepvar2 = factor(indepvar2, 
                            labels = c("Upper Chamber Size",
                                       "Log of Lower Chamber Size",
                                       "Lower Chamber Size"
                                                  ))) %>%
  ggplot(aes(x = indepvar2)) +
    geom_bar(color = "black") +
  labs(x = "",
       y = "") +
  coord_flip() +
  theme_bw()

# Histogram of the Coefficients and the Standard Errors
# Coefficients:
dat %>%
  ggplot(aes(x = coef)) +
  geom_histogram(bins = 15, color = "black") +
  labs(x = "Coefficients", y = '') +
  theme_bw()

# Standard Errors
dat %>%
  ggplot(aes(x = SE)) +
  geom_histogram(bins = 10, color = "black") +
  labs(x = "Standard Errors", y = '') +
  theme_bw()

# Sign Coefficients
dat %>%
  ggplot(aes(x=as.factor(scoef))) +
  geom_bar(color = "black") +
  labs(x = "Coefficient Sign",
       y = "") +
  theme_bw()

## Descriptive Statistics of Moderators
fulldat$usemeta2 <- factor(fulldat$usemeta)
levels(fulldat$usemeta2) <- c("Extended Sample", "Main Sample")
aux <- select(fulldat, usemeta2, indepvar2, year, published, method,
              instdesign, elecsys2) %>%
  mutate(indepvar2 = recode(indepvar2,
                            `N` = "Lower Chamber Size",
                            `K` = "Upper Chamber Size",
                            `logN` = "Log of Lower Chamber Size"),
         elecsys2 = recode(elecsys2,
                           `Non-Maj` = "Non-Majoritarian",
                           `Maj` = "Majoritarian")) %>%
  rename(`Independent Variables` = indepvar2,
         `Year`                  = year,
         `Published work`        = published,
         `Estimation method`     = method,
         `Institutional Design`  = instdesign,
         `Electoral system`      = elecsys2)
export2md(descrTable(~.-usemeta2, aux, y = aux$usemeta2,
                        show.p.overall = F, show.all = T),
          caption = "Descriptive Statistics of Moderators",
          format  = "latex")

## Binomial Tests for Coefficient Signs

# For the number of legislators in the lower house, the results follow below.
aux <- filter(dat, indepvar2 == "N")
aux2 <- binom.test(table(aux$scoef)[2], sum(table(aux$scoef)), p = 0.5)
pander(tidy(aux2)[,-c(2,4,5,6)][,c(3,4,1,2)])

# For the log of the number of legislators in the lower house, the results are as follows:
aux <- filter(dat, indepvar2 == "logN")
aux2 <- binom.test(table(aux$scoef)[2], sum(table(aux$scoef)), p = 0.5)
pander(tidy(aux2)[,-c(2,4,5,6)][,c(3,4,1,2)])

# For the number of legislators in the upper house, the results are:
aux <- filter(dat, indepvar2=='K')
aux2 <- binom.test(table(aux$scoef)[2], sum(table(aux$scoef)), p=0.5)
pander(tidy(aux2)[,-c(2,4,5,6)][,c(3,4,1,2)])

# Unicameral vesus non-unicameral systems:
aux <- list()
aux[[1]] <- filter(dat, indepvar2 == "N", instdesign2 == 'Unicameral')
aux[[2]] <- filter(dat, indepvar2 == "N", instdesign2 != 'Unicameral')
aux[[3]] <- filter(dat, indepvar2 == "logN", instdesign2 == 'Unicameral')
aux[[4]] <- filter(dat, indepvar2 == "logN", instdesign2 != 'Unicameral')
aux[[1]] <- binom.test(table(aux[[1]]$scoef)[2],
                       sum(table(aux[[1]]$scoef)), p = 0.5)
aux[[2]] <- binom.test(table(aux[[2]]$scoef)[2],
                       sum(table(aux[[2]]$scoef)), p = 0.5)
aux[[3]] <- binom.test(table(aux[[3]]$scoef)[2],
                       sum(table(aux[[3]]$scoef)), p = 0.5)
aux[[4]] <- binom.test(table(aux[[4]]$scoef)[2],
                       sum(table(aux[[4]]$scoef)), p = 0.5)
aux2 <- tidy(aux[[1]])[,-c(2,4,5,6)][,c(3,4,1,2)]
aux2 <- bind_rows(aux2, tidy(aux[[2]])[,-c(2,4,5,6)][,c(3,4,1,2)])
aux2 <- bind_rows(aux2, tidy(aux[[3]])[,-c(2,4,5,6)][,c(3,4,1,2)])
aux2 <- bind_rows(aux2, tidy(aux[[4]])[,-c(2,4,5,6)][,c(3,4,1,2)])
aux2 <- aux2 %>%
  mutate(`Indep. Variable` = c('Lower House Size', 'Lower House Size',
                               'Log of Lower House Size', 'Log of Lower House Size'),
         `Legislative Inst.` = c('Unicameral', 'Non-unicameral',
                                 'Unicameral', 'Non-unicameral')) %>%
  relocate(`Indep. Variable`, `Legislative Inst.`, .before = method)
aux2$method <- NULL
aux2$alternative <- NULL
pander(aux2)

## Meta-Analysis

# Lower House Size and Expenditure per Capita
aux <- dat %>%
  filter(indepvar2 == 'N',
         depvar2 == 'ExpPC')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
f1 <- build_forest(mod, NULL, lsize = 15, ttl = '1.1 - Lower Chamber Size\nand Expenditure Per Capita')
funnel(mod)

# Electoral System Subgroup Analysis
build_forest_het(aux, aux$coef, aux$VAR, slab = aux$authoryear,
                 capt = NULL, hetvar = aux$elecsys2)
                 capt = NULL, hetvar = aux$elecsys2)

# Institutional Design Subgroup Analysis
build_forest_het(aux, aux$coef, aux$VAR, slab = aux$authoryear,
                 capt = NULL, hetvar = aux$instdesign2)

# Upper House Size and Expenditure per Capita
aux <- dat %>%
  filter(indepvar2 == 'K',
         depvar2 == 'ExpPC')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)
f2 <- build_forest(mod, NULL, 15, ttl = '1.6 - Upper Chamber Size\nand Expenditure Per Capita')

# Lower House Size and Log Expenditure Per Capita
aux <- dat %>%
  filter(indepvar2 == 'N',
         depvar2 == 'logExpPC')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)
f3 <- build_forest(mod, NULL, 15, ttl = '1.2 - Lower Chamber Size\nand Log Expenditure Per Capita')

# Log of Lower House Size and Log of Expenditure Per Capita
aux <- dat %>%
  filter(indepvar2 == 'logN',
         depvar2 == 'logExpPC')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)
f4  <- build_forest(mod, NULL, 15, ttl = '1.4 - Log Lower Chamber Size\nand Log Expenditure Per Capita')

# Lower House Size and Expenditure as Percentage of GDP
aux <- dat %>%
  filter(indepvar2 == 'N',
         depvar2 == 'PCTGDP')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)
f5 <- build_forest(mod, NULL, 15, ttl = '1.3 - Lower Chamber Size\nand Expenditure as Percentage of GDP')

# Log Lower House Size and Expenditure as Percentage of GDP
aux <- dat %>%
  filter(indepvar2 == 'logN',
         depvar2 == 'PCTGDP')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")

mod
build_forest(mod, NULL)
funnel(mod)
f6 <- build_forest(mod, NULL, 15, ttl = '1.5 - Log Lower Chamber Size\nand Expenditure as Percentage of GDP')

# Upper House Size and Expenditure as Percentage of GDP
aux <- dat %>%
  filter(indepvar2 == 'K',
         depvar2 == 'PCTGDP')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)
f7 <- build_forest(mod, NULL, 15, ttl = '1.7 - Upper Chamber Size\nand Expenditure as Percentage of GDP')

# Lower House Size and Expenditure per Capita (IV)
aux <- dat %>%
  filter(indepvar2 == 'N',
         depvar2 == 'ExpPC',
         method %in% c('IV'))

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)

aux <- dat %>%
  filter(indepvar2 == 'N',
         depvar2 == 'ExpPC')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")

f8 <- build_forest_het(aux, aux$coef, 
                       aux$VAR, slab = aux$authoryear,
                       capt = NULL, lsize = 15, 
                       ttl = 'Lower House Size and Expenditure per Capita\n(subgrouping by the estimation technique)', 
                       hetvar = aux$method)

# Regression Method Subgroup Analysis
aux <- dat %>%
  filter(indepvar2 == 'N',
         depvar2 == 'ExpPC')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")

build_forest_het(aux, aux$coef, aux$VAR, slab = aux$authoryear,
                 capt = NULL, hetvar = aux$method)

# Lower House Size and Log of Expenditure per Capita (RDD)
aux <- dat %>%
  filter(indepvar2 == 'N',
         depvar2 == 'logExpPC',
         method == 'RDD')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)

aux <- dat %>%
  filter(indepvar2 == 'N',
         depvar2 == 'logExpPC',
         method == 'RDD')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")

f9 <- build_forest(mod, capt = NULL, 
                       lsize = 15, ttl = 'Lower House Size and Log of Expenditure per Capita\n(subgrouping by RDD)')

# Plot 1
pdf('../graphs/graph1.pdf', width = 16, height = 11)
ggarrange(f1,f3,f5,f4,f6,f2,f7, align = 'hv')
dev.off()
## Plot 2
pdf('../graphs/graph2.pdf', width = 12, height = 6)
ggarrange(f8, f9, align = 'hv')
dev.off()

## Meta-Analysis (All Coefficients)

# Lower House Size and Expenditure Per Capita
aux <- fulldat %>%
  filter(indepvar2 == 'N',
         depvar2 == 'ExpPC')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)

# Electoral System Subgroup Analysis
build_forest_het(aux, aux$coef, aux$VAR, slab = aux$authoryear,
                 capt = NULL, hetvar = aux$elecsys2)

# Upper House Size and Expenditure Per Capita
aux <- fulldat %>%
  filter(indepvar2 == 'K',
         depvar2 == 'ExpPC')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)

# Lower House Size and Log of Expenditure Per Capita
aux <- fulldat %>%
  filter(indepvar2 == 'N',
         depvar2 == 'logExpPC')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)

# Log of Lower House Size and Log of Expenditure Per Capita
aux <- fulldat %>%
  filter(indepvar2 == 'logN',
         depvar2 == 'logExpPC')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)

# Lower House Size and Expenditure as Percentage of GDP
aux <- fulldat %>%
  filter(indepvar2 == 'N',
         depvar2 == 'PCTGDP')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)

# Log of Lower House Size and Expenditure as Percentage of GDP
aux <- fulldat %>%
  filter(indepvar2 == 'logN',
         depvar2 == 'PCTGDP')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)

# Upper House Size and Expenditure as Percentage of GDP
aux <- fulldat %>%
  filter(indepvar2 == 'K',
         depvar2 == 'PCTGDP')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)

# Lower House Size and Expenditure per Capita (IV)
aux <- fulldat %>%
  filter(indepvar2 == 'N',
         depvar2 == 'ExpPC',
         method %in% c('IV'))

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)

# Lower House Size and Log of Expenditure per Capita (RDD)
aux <- fulldat %>%
  filter(indepvar2 == 'N',
         depvar2 == 'logExpPC',
         method == 'RDD')

mod <- rma.mv(coef, VAR, data=aux,
          slab=paste(authoryear),
          test = 't',
          random = ~ 1 | id_level1/id_level2, 
          tdist = TRUE,
          method = "REML")
mod
build_forest(mod, NULL)
funnel(mod)

## Meta-Regressions

# Meta-Regressions for Expenditure Per Capita
# Restricted Sample
mod <- rma.mv(yi = coef,
              V = VAR,
              data = dat,
              method = "REML",
              random = ~ 1 | id_level1/id_level2, 
              mods = ~indepvar2+year+published+elecsys2+method+instdesign,
              test = "knha",
              sparse = TRUE,
              tdist = TRUE,
              subset = dat$depvar2=='ExpPC',
              slab = dat$authoryear)

summary(mod)

# Full Sample
mod <- rma.mv(yi = coef,
              V = VAR,
              data = fulldat,
              method = "REML",
              random = ~ 1 | id_level1/id_level2, 
              mods = ~indepvar2+year+published+elecsys2+method+instdesign,
              test = "knha",
              tdist = TRUE,
              sparse = TRUE,
              subset = fulldat$depvar2=='ExpPC',
              slab = fulldat$authoryear)

summary(mod)

# Meta-Regressions for Log of Expenditure Per Capita
# Restricted Sample
mod <- rma.mv(yi = coef,
              V = VAR,
              data = dat,
              method = "REML",
              random = ~ 1 | id_level1/id_level2, 
              mods = ~indepvar2+year+published+elecsys2+method+instdesign,
              test = "knha",
              tdist = TRUE,
              sparse = TRUE,
              subset = dat$depvar2=='logExpPC',
              slab = dat$authoryear)

summary(mod)

# Full Sample
mod <- rma.mv(yi = coef,
              V = VAR,
              data = fulldat,
              method = "REML",
              random = ~ 1 | id_level1/id_level2, 
              mods = ~indepvar2+year+published+elecsys2+method+instdesign,
              test = "knha",
              tdist = TRUE,
              sparse = TRUE,
              subset = fulldat$depvar2=='logExpPC',
              slab = fulldat$authoryear)

summary(mod)

# Meta-Regressions for Expenditure as a Percentage of the GDP
# Restricted Sample
mod <- rma.mv(yi = coef,
              V = VAR,
              data = dat,
              method = "REML",
              random = ~ 1 | id_level1/id_level2, 
              mods = ~indepvar2+year+published+elecsys2+method+instdesign,
              test = "knha",
              tdist = TRUE,
              sparse = TRUE,
              subset = dat$depvar2=='PCTGDP',
              slab = dat$authoryear)

summary(mod)

# Full Sample
mod <- rma.mv(yi = coef,
              V = VAR,
              data = fulldat,
              method = "REML",
              random = ~ 1 | id_level1/id_level2, 
              mods = ~indepvar2+year+published+elecsys2+method+instdesign,
              test = "knha",
              tdist = TRUE,
              sparse = TRUE,
            subset = fulldat$depvar2=='PCTGDP',
            slab = fulldat$authoryear)
summary(mod)

## Meta-Regressions (All Coefficients)
# Restricted Sample
mod <- rma.mv(yi = coef,
              V = VAR,
              data = dat,
              method = "REML",
              random = ~ 1 | id_level1/id_level2, 
              mods = ~depvar2+indepvar2+year+published+elecsys2+method+instdesign,
              test = "knha",
              tdist = TRUE,
              sparse = TRUE,
              slab = dat$authoryear)

summary(mod)

# Full Sample
mod <- rma.mv(yi = coef,
              V = VAR,
              data = fulldat,
              method = "REML",
              random = ~ 1 | id_level1/id_level2, 
              mods = ~depvar2+indepvar2+year+published+elecsys2+method+instdesign,
              test = "knha",
              tdist = TRUE,
              sparse = TRUE,
              slab = fulldat$authoryear)
summary(mod)
