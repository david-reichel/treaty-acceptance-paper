########################################
######   Replication code for  #########
########################################
### Reichel, David, Vink, Maarten and Grimheden, Jonas
### Regional Cooperation and EU conditionality: Council of Europe Treaty Ratification, 1949-2016
### In: Journal of European Public Policy, 2019
########################################

### Packages used (see below session info on package versions)
library(tidyverse)
library(broom)
library(survival)
library(coxme)
library(arm)

###################
# function on median hazard function
# from Austin 2017
mhrfun <- function(md) {
  dfv <- as.data.frame(VarCorr(md))
  mhrs <- dfv %>%
    mutate_all(function(x) exp(sqrt(x) * qnorm(0.75) * sqrt(2)))
  mhrs
}

# function for running models and extracting information needed
runmodel <- function(form = NULL, d = NULL, m = NULL){
  cmm1 <- coxme(f, data = d, ties = c("efron"))
  pred <- ifelse(cmm1$linear.predictor > 0, 1, 0)
  obs <- d$rtf
  mhrs <- mhrfun(cmm1)
  df1 <- as.data.frame(coef(cmm1))
  df1$term <- row.names(df1)
  df1$se <- sqrt(diag(vcov(cmm1)))
  df1$se95 <- 1.96*sqrt(diag(vcov(cmm1)))
  df1$se99 = 2.58*sqrt(diag(vcov(cmm1)))
  df1$model = m
  df1$n <- cmm1$n[2]
  df1$var_ISO2 <- as.data.frame(VarCorr(cmm1))[1, 1]
  df1$var_no <- as.data.frame(VarCorr(cmm1))[1, 2]
  df1$mhr_ISO2 <- mhrs$ISO2
  df1$mhr_no <- mhrs$no
  df1
}

runmodel2 <- function(form = NULL, d = NULL, m = NULL){
  coxme(f, data = d, ties = c("efron"))
}


dat2 <- read.csv("coe_treaty_ratifications_final_dataset1.csv", stringsAsFactors = FALSE)
dim(dat2)


#############################


# descriptives
t1 <- dat2 %>%
  dplyr::select(rtf, upcomingEU, propneighrtf1, polity2, ln_hudoc_violations, rats_needed, 
                recent_independence, n_signed_upon_opening_treaty_std) %>%
  summarise_all(funs(mean,min,max)) %>%
  gather(variable, value) %>%
  mutate(stt = str_extract(variable, "_[^_]+$"),
         stt = str_replace(stt, "_", ""),
         variable = str_replace(variable, "_[^_]+$", ""),
         value = round(value, 3)) %>%
  spread(key = stt, value = value) %>%
  dplyr::select(variable, min, mean, max)
t1$n <- nrow(dat2)
t1$no <- length(unique(dat2$no))
t1$nc <- length(unique(dat2$ISO2))

t1


 #----------------------------------------
##### how often do countries ratify together

dat2 %>%
  group_by(propneighrtf1) %>%
  summarise(sm_rtf = sum(rtf)) %>%
  mutate(prp = sm_rtf / sum(sm_rtf))

dat2 %>%
  group_by(propneighrtf1, ISO2) %>%
  summarise(sm_rtf = sum(rtf)) %>%
  spread(propneighrtf1, sm_rtf) %>%
  mutate(sm = `1`+`0`,
         prp = `1`/sm) %>%
  arrange(desc(prp))

## how many happen in the first years after opening
dat2 %>%
  filter(rtf == TRUE) %>%
  count(end) %>%
  mutate(prp = n / sum(n),
         cprp = cumsum(prp))



## Run regression models
## ------------------------------------------------------------------------------
# --------------------------------------------

d <- dat2

f <- as.formula(Surv(start, end, rtf) ~ 1 + (1 | ISO2) + (1 | no))

m <- "Model 0"

m0 <- runmodel2(form = f, d = d, m = m)

df0 <- mhrfun(m0) %>%
  rename(mhr_ISO2 = ISO2, mhr_no = no)
temp <- as.data.frame(VarCorr(m0)) %>% rename(var_ISO2 = ISO2, var_no = no)
df0 <- cbind(df0, temp) %>% mutate(model = m, n = nrow(d))
rm(temp)

# -

m <- "Model 0 HR"

d <- filter(dat2, type2 == "Human rights (core)")

m0hr <- runmodel2(form = f, d = d, m = m)

df0hr <- mhrfun(m0hr) %>%
  rename(mhr_ISO2 = ISO2, mhr_no = no)
temp <- as.data.frame(VarCorr(m0hr)) %>% rename(var_ISO2 = ISO2, var_no = no)
df0hr <- cbind(df0hr, temp) %>% mutate(model = m, n = nrow(d))
rm(temp)

# -

m <- "Model 0 HA"

d <- filter(dat2, type2 == "Harmonisation")

m0ha <- runmodel2(form = f, d = d, m = m)

m0ha <- mhrfun(m0ha) %>%
  rename(mhr_ISO2 = ISO2, mhr_no = no)
temp <- as.data.frame(VarCorr(m0ha)) %>% rename(var_ISO2 = ISO2, var_no = no)
m0ha <- cbind(m0ha, temp) %>% mutate(model = m, n = nrow(d))
rm(temp)


# -------------------
d <- dat2

f <- as.formula(Surv(start, end, rtf, type = "counting") ~ upcomingEU + 
                  propneighrtf1 + polity2 + recent_independence + ln_hudoc_violations +
                  type2 + 
                  (1 | ISO2) + (1 | no))
m <- "Model 1"

df1 <- runmodel(form = f, d = d, m = m)
#m1 <- runmodel2(form = f, d = d, m = m)

# -------------------

f <- as.formula(Surv(start, end, rtf, type = "counting") ~ upcomingEU + 
                  propneighrtf1 + polity2 + recent_independence + ln_hudoc_violations +
                  (1 | ISO2) + (1 | no))
m <- "Model 1 HR"
d <- filter(dat2, type2 == "Human rights (core)")

df1hr <- runmodel(form = f, d = d, m = m)
#m1hr <- runmodel2(form = f, d = d, m = m)

# -------------------

f <- as.formula(Surv(start, end, rtf, type = "counting") ~ upcomingEU + 
                  propneighrtf1 + polity2 + recent_independence + ln_hudoc_violations +
                  (1 | ISO2) + (1 | no))
m <- "Model 1 HA"
d <- filter(dat2, type2 == "Harmonisation")

df1ha <- runmodel(form = f, d = d, m = m)
#m1ha <- runmodel2(form = f, d = d, m = m)

# -------------------


f <- as.formula(Surv(start, end, rtf, type = "counting") ~ upcomingEU + 
                  propneighrtf1*polity2 + recent_independence + ln_hudoc_violations +
                  type2 + 
                  (1 | ISO2) + (1 | no))
m <- "Model 2"
d <- dat2

df2 <- runmodel(form = f, d = d, m = m)
#m2 <- runmodel2(form = f, d = d, m = m)

# -------------------

f <- as.formula(Surv(start, end, rtf, type = "counting") ~ upcomingEU + 
                  propneighrtf1*polity2 + recent_independence + ln_hudoc_violations + 
                  (1 | ISO2) + (1 | no))
m <- "Model 2 HR"
d <- filter(dat2, type2 == "Human rights (core)")

df2hr <- runmodel(form = f, d = d, m = m)


# -------------------

f <- as.formula(Surv(start, end, rtf, type = "counting") ~ upcomingEU + 
                  propneighrtf1*polity2 + recent_independence + ln_hudoc_violations +
                  (1 | ISO2) + (1 | no))
m <- "Model 2 HA"
d <- filter(dat2, type2 == "Harmonisation")

df2ha <- runmodel(form = f, d = d, m = m)


# -----------------------------------------------------------
# model 3
# -------------------


f <- as.formula(Surv(start, end, rtf, type = "counting") ~ upcomingEU + 
                  propneighrtf1*polity2 + recent_independence + ln_hudoc_violations +
                  type2 + rats_needed + n_signed_upon_opening_treaty + 
                  (1 | ISO2) + (1 | no))
m <- "Model 3"
d <- dat2

df3 <- runmodel(form = f, d = d, m = m)


# -------------------

f <- as.formula(Surv(start, end, rtf, type = "counting") ~ upcomingEU + 
                  propneighrtf1*polity2 + recent_independence + ln_hudoc_violations +
                  rats_needed + n_signed_upon_opening_treaty + 
                  (1 | ISO2) + (1 | no))
m <- "Model 3 HR"
d <- filter(dat2, type2 == "Human rights (core)")

df3hr <- runmodel(form = f, d = d, m = m)


# -------------------

f <- as.formula(Surv(start, end, rtf, type = "counting") ~ upcomingEU + 
                  propneighrtf1*polity2 + recent_independence + ln_hudoc_violations +
                  rats_needed + n_signed_upon_opening_treaty + 
                  (1 | ISO2) + (1 | no))
m <- "Model 3 HA"
d <- filter(dat2, type2 == "Harmonisation")

df3ha <- runmodel(form = f, d = d, m = m)
#m3ha <- runmodel2(form = f, d = d, m = m)


# -----------------------------------------------------------------------
# -----------------------------------------------------------------------


## results

df <- rbind(df1, df1hr, df1ha, df2, df2hr, df2ha, df3, df3hr, df3ha)

df <- df %>%
  rename(estimate = `coef(cmm1)`) %>%
  mutate(estimate.exp = exp(estimate),
         lo95 = exp(estimate - se95),
         hi95 = exp(estimate + se95),
         lo99 = exp(estimate - se99),
         hi99 = exp(estimate + se99),
         modeln = paste0(model, " (n = ", n, ")")) %>%
  mutate(term = str_replace(term, "upcomingEU", "Upcoming EU membership"),
         term = str_replace_all(term, "_", " "),
         term = str_replace(term, "propneighrtf1:polity2", "Interaction: rat. neighbour * polity"),         
         term = str_replace(term, "propneighrtf1", "Ratification in neighbouring country"),
         term = str_replace(term, "type2Foundational", "Type: Foundational"),
         term = str_replace(term, "type2Human rights", "Type: Human Rights"),
         term = str_replace(term, "polity2", "Democracy (polity)"),
         term = str_replace(term, "recent independence", "Recent independence"),
         term = str_replace(term, "ln hudoc violations", "Number of ECHR violations (log)"),
         term = str_replace(term, "n signed upon opening treaty", "Number of signatures on opening day"),
         term = str_replace(term, "rats needed", "Number of ratifications needed")) %>%
  mutate(term = as.factor(term)) %>%
  mutate(term = factor(term, levels = rev(levels(term))),
         term = fct_relevel(term, "Upcoming EU membership", "Ratification in neighbouring country",
                            "Democracy (polity)", "Interaction: rat. neighbour * polity"))

tdf2 <- df %>%
  dplyr::select(model, term, estimate, se, estimate.exp, lo95, hi95) %>%
  mutate(estimate = round(estimate, 2),
         se = round(se, 2),
         estimate.exp = round(estimate.exp, 2),
         lo95 = round(lo95, 2),
         hi95 = round(hi95, 2)) %>% 
  mutate(term = as.factor(term)) %>%
  mutate(term = factor(term, levels = rev(levels(term))),
         term = fct_relevel(term, "Upcoming EU membership", "Ratification in neighbouring country",
                            "Democracy (polity)", "Interaction: rat. neighbour * polity",
                            "Recent independence")) %>%
  mutate(result = paste0(estimate.exp, " (", lo95, "-", hi95, ")")) %>%
  dplyr::select(model, term, result) %>%
  spread(model, result)
tdf2

# write.csv2(tdf2, file = paste0("final_results_models_cox_", Sys.Date(), "b.csv"), row.names = FALSE, na = "n/a")
# I use write.csv2 because I am located in Europe and can then check in spreadsheet format

df_sub <- df[!duplicated(df[ , c("model", "n")]), ]


## visualise coefficients (not presented in the paper)
p_cox <- ggplot(df, aes(estimate, term)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_segment(aes(x = estimate - se99, xend = estimate + se99, y = term, yend = term),
               size = 1, alpha = 0.6, colour = "darkgrey") +
  geom_segment(aes(x = estimate - se95, xend = estimate + se95, y = term, yend = term),
               size = 1.5, alpha = 1, colour = "darkgrey") +
  geom_point(size = 2, colour = "darkgrey") +
  geom_text(aes(label = round(estimate, 2), vjust = 1.3),
             size = 3) +
  facet_wrap(~ modeln, ncol = 3) +
  labs(title = "",
       subtitle = "",
       x = "Coefficient",
       y = "") +
  theme_bw()
p_cox

# ggsave(p_cox, file = "graphs/cox_results_final.jpg", width = 8, height = 6)


df_modinfo <- rbind(df1, df1ha, df1hr, df2, df2ha, df2hr, df3, df3ha, df3hr) %>%
  dplyr::select(model, var_ISO2, mhr_ISO2, var_no, mhr_no, n) %>%
  rbind(df0, df0hr) %>%
  filter(!duplicated(.)) %>%
  t(.) %>% 
  as.data.frame(header = TRUE)

df_modinfo

# write.csv2(df_modinfo, file = "model_information.csv")


length(unique(dat2$no))
length(unique(dat2$ISO2))


## END

###########################################################################

# R version 3.5.2 (2018-12-20)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows >= 8 x64 (build 9200)
# 
# Matrix products: default
# 
# locale:
#   LC_COLLATE=German_Austria.1252  LC_CTYPE=German_Austria.1252   
#   LC_MONETARY=German_Austria.1252 LC_NUMERIC=C                   
#   LC_TIME=German_Austria.1252    
# 
# attached base packages:
#   stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] bindrcpp_0.2.2  arm_1.10-1      lme4_1.1-19     Matrix_1.2-15  
# [5] MASS_7.3-51.1   coxme_2.2-10    bdsmatrix_1.3-3 survival_2.43-3
# [9] broom_0.5.1     forcats_0.3.0   stringr_1.3.1   dplyr_0.7.8    
# [13] purrr_0.2.5     readr_1.3.1     tidyr_0.8.2     tibble_2.0.1   
# [17] ggplot2_3.1.0   tidyverse_1.2.1
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_0.2.5 splines_3.5.2    haven_2.0.0      lattice_0.20-38 
# [5] colorspace_1.3-2 generics_0.0.2   yaml_2.2.0       utf8_1.1.4      
# [9] rlang_0.3.1      pillar_1.3.1     nloptr_1.2.1     glue_1.3.0      
# [13] withr_2.1.2      modelr_0.1.2     readxl_1.2.0     bindr_0.1.1     
# [17] plyr_1.8.4       munsell_0.5.0    gtable_0.2.0     cellranger_1.1.0
# [21] rvest_0.3.2      coda_0.19-2      labeling_0.3     fansi_0.4.0     
# [25] Rcpp_1.0.0       scales_1.0.0     backports_1.1.3  jsonlite_1.6    
# [29] abind_1.4-5      hms_0.4.2        stringi_1.2.4    grid_3.5.2      
# [33] cli_1.0.1        tools_3.5.2      magrittr_1.5     lazyeval_0.2.1  
# [37] crayon_1.3.4     pkgconfig_2.0.2  xml2_1.2.0       lubridate_1.7.4 
# [41] assertthat_0.2.0 minqa_1.2.4      httr_1.4.0       rstudioapi_0.8  
# [45] R6_2.3.0         nlme_3.1-137     compiler_3.5.2  


########################################################################

