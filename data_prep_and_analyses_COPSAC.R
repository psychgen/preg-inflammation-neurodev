################################################################################
######### Pregnancy diet neurodevelopment - Trio models ########################
################################################################################

#Load data
load(file = "//.../food_preference.Rdata")

##Test of correlation between child and mother diet 
hist(pref_data$Diet_PC2)
cor(pref_data$Diet_PC2, pref_data$Diet_PC2_child, use = "complete.obs") 

#Test for correlation between ADHD PRS and maternal dietary pattern PRS
hist(pref_data$PC1_hirshorn_mor)
cor.test(pref_data$Diet_PC2, pref_data$PC1_hirshorn_mor, use = "complete.obs") 

#Test for correlation between DietPC2 and dietary pattern PRS
hist(pref_data$adhd_prs_mor)
cor.test(pref_data$adhd_prs_mor, pref_data$PC1_hirshorn_mor, use = "complete.obs") 


##Confounder table

## to include N
# glm
getIRR <- function(fit){
  data.frame(IRR = exp(coef(summary(fit)))[2,1] %>% round(2),
             CI.lower = exp(confint(fit))[2,1] %>% round(2),
             CI.upper = exp(confint(fit))[2,2] %>% round(2),
             p = coef(summary(fit))[2,4] %>% round(15),
             n = nrow(fit$model),
             beta=coef(summary(fit))[2,1] %>% round(3),
             se=coef(summary(fit))[2,2] %>% round(3),
             outcome = names(fit$model)[1])
}

# linear
getIRR_lm <- function(fit){
  data.frame(beta=coef(summary(fit))[2,1]%>% round(2),
             CI.lower = confint(fit)[2,1] %>% round(3),
             CI.upper = confint(fit)[2,2] %>% round(3),
             p = coef(summary(fit))[2,4] %>% as.numeric(as.character(formatC(., format = "e", digits = 15))),
             n = nrow(fit$model),
             se=coef(summary(fit))[2,2]%>% round(3),
             outcome = names(fit$model)[1])
}


#install.packages('xlsx')
library("xlsx")
#install.packages('rJava')
library(rJava)
detach(package:dplyr)
library(dplyr)

##Hirschhorn PC1

# Function
sigfunction <- function(data, col_of_interest, ID){
  
  # Determine ID columns and effect column
  interest_idx <- which(colnames(data) == col_of_interest) #interest_idx is the exposure
  colnames(data)[interest_idx] <- "effect"
  ID_idx <- which(colnames(data) == ID)
  colnames(data)[ID_idx] <- "ABCNO"
  
  # sex 
  sex <- glm(effect ~ SEX,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Sex",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)

  # Birth weight
  birthw <- glm(effect ~ bw,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Child birth weight",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Gestational age
  gest_age <- glm(effect ~ AGE_DAYS,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Child gestational age",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Maternal smoking in pregnancy
  mat_smoke <- glm(effect ~ smoking_preg,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal smoking in pregnancy",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Maternal alcohol intake in pregnancy 
  mat_alc <- glm(effect ~ alc_preg_YN_2010,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal alcohol intake in pregnancy",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Parity
  mat_parity <- glm(effect ~ as.numeric(parity),data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal parity",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)  
  
  # Maternal antibiotics usage in pregnancy
  mater_ab <- glm(effect ~ AB_preg,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal use of antibiotics in pregnancy",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Maternal preeclampsia in pregnancy 
  mater_pe <- glm(effect ~ PE,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal pre-eclampsia in pregnancy",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Paternal age at child birth
  pater_age <- glm(effect ~ Fatherage,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Paternal age at birth",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # maternal week 24 CRP
  crpmother <- glm(effect ~ log2(CRPm),data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal CRP at pregnancy week 24",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # maternal week 24 glyca
  glycamother <- glm(effect ~ glyca24,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal glycA at pregnancy week 24",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Maternal BMI
  mat_weight <- glm(effect ~ Mother_BMI,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal pre-pregnancy BMI",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Maternal age
  mater_age <- glm(effect ~ Motherage,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal age at child birth",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Maternal education
  mater_edu <- glm(effect ~ mat_edu,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal education at child birth",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Paternal education
  pater_edu <- glm(effect ~ pat_edu,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Paternal education at child birth",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Household income
  house_income <- glm(effect ~ income,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Household income at child birth",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)

  
  results <- bind_rows(sex, birthw, gest_age, mat_smoke, mat_alc, mat_parity, mater_ab, mater_pe, pater_age, crpmother, glycamother, mat_weight, mater_age, mater_edu, pater_edu, house_income)
  
  
  return(results)
}


conf_pc1_food <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_child", ID="ABCNO")


conf_pc1_food_mor <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_mor", ID="ABCNO")


conf_pc1_food_far <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_far", ID="ABCNO")


##Including maternal education, paternal education and income as categorical variables

#Maternal PRS

pref_data_inc <- pref_data %>% filter(!is.na(income) & !is.na(PC1_hirshorn_mor))
pref_data_inc<- pref_data_inc %>% mutate(income = relevel(factor(income), ref="2"))

glm(PC1_hirshorn_mor ~ income, data = pref_data_inc) %>% tidy(conf.int=TRUE)
#Likehood ratio test
fit0 <- lm(PC1_hirshorn_mor ~ 1, data=pref_data_inc %>% filter(!is.na (income)))
fit1 <- lm (PC1_hirshorn_mor ~ income, data=pref_data_inc)
anova(fit0, fit1, test="Chisq")

pref_data_edu <- pref_data %>% filter(!is.na(mat_edu) & !is.na(PC1_hirshorn_mor))

glm(PC1_hirshorn_mor ~ mat_edu, data = pref_data_edu) %>% tidy(conf.int=TRUE)
#Likehood ratio test
fit0 <- lm(PC1_hirshorn_mor ~ 1, data=pref_data_edu)
fit1 <- lm (PC1_hirshorn_mor ~ mat_edu, data=pref_data_edu)
anova(fit0, fit1, test="Chisq")

pref_data_edup <- pref_data %>% filter(!is.na(pat_edu) & !is.na(PC1_hirshorn_mor))
pref_data_edup <- pref_data_edup %>% mutate(pat_edu = relevel(factor(pat_edu), ref="1"))

glm(PC1_hirshorn_mor ~ pat_edu, data = pref_data_edup) %>% tidy(conf.int=TRUE)
#Likehood ratio test
fit0 <- lm(PC1_hirshorn_mor ~ 1, data=pref_data_edup) 
fit1 <- lm (PC1_hirshorn_mor ~ pat_edu, data=pref_data_edup)
anova(fit0, fit1, test="Chisq")

#Paternal PRS

pref_data_inc <- pref_data %>% filter(!is.na(income) & !is.na(PC1_hirshorn_far))
pref_data_inc<- pref_data_inc %>% mutate(income = relevel(factor(income), ref="2"))

glm(PC1_hirshorn_far ~ income, data = pref_data_inc) %>% tidy(conf.int=TRUE)
#Likehood ratio test
fit0 <- lm(PC1_hirshorn_far ~ 1, data=pref_data_inc %>% filter(!is.na (income)))
fit1 <- lm (PC1_hirshorn_far ~ income, data=pref_data_inc)
anova(fit0, fit1, test="Chisq")

pref_data_edu <- pref_data %>% filter(!is.na(mat_edu) & !is.na(PC1_hirshorn_far))

glm(PC1_hirshorn_far ~ mat_edu, data = pref_data_edu) %>% tidy(conf.int=TRUE)
#Likehood ratio test
fit0 <- lm(PC1_hirshorn_far ~ 1, data=pref_data_edu)
fit1 <- lm (PC1_hirshorn_far ~ mat_edu, data=pref_data_edu)
anova(fit0, fit1, test="Chisq")

pref_data_edup <- pref_data %>% filter(!is.na(pat_edu) & !is.na(PC1_hirshorn_far))
pref_data_edup <- pref_data_edup %>% mutate(pat_edu = relevel(factor(pat_edu), ref="1"))

glm(PC1_hirshorn_far ~ pat_edu, data = pref_data_edup) %>% tidy(conf.int=TRUE)
#Likehood ratio test
fit0 <- lm(PC1_hirshorn_far ~ 1, data=pref_data_edup) 
fit1 <- lm (PC1_hirshorn_far ~ pat_edu, data=pref_data_edup)
anova(fit0, fit1, test="Chisq")


#Child PRS

pref_data_inc <- pref_data %>% filter(!is.na(income) & !is.na(PC1_hirshorn_child))
pref_data_inc<- pref_data_inc %>% mutate(income = relevel(factor(income), ref="2"))

glm(PC1_hirshorn_child ~ income, data = pref_data_inc) %>% tidy(conf.int=TRUE)
#Likehood ratio test
fit0 <- lm(PC1_hirshorn_child ~ 1, data=pref_data_inc %>% filter(!is.na (income)))
fit1 <- lm (PC1_hirshorn_child ~ income, data=pref_data_inc)
anova(fit0, fit1, test="Chisq")

pref_data_edu <- pref_data %>% filter(!is.na(mat_edu) & !is.na(PC1_hirshorn_child))

glm(PC1_hirshorn_child ~ mat_edu, data = pref_data_edu) %>% tidy(conf.int=TRUE)
#Likehood ratio test
fit0 <- lm(PC1_hirshorn_child ~ 1, data=pref_data_edu)
fit1 <- lm (PC1_hirshorn_child ~ mat_edu, data=pref_data_edu)
anova(fit0, fit1, test="Chisq")

pref_data_edup <- pref_data %>% filter(!is.na(pat_edu) & !is.na(PC1_hirshorn_child))
pref_data_edup <- pref_data_edup %>% mutate(pat_edu = relevel(factor(pat_edu), ref="1"))

glm(PC1_hirshorn_child ~ pat_edu, data = pref_data_edup) %>% tidy(conf.int=TRUE)
#Likehood ratio test
fit0 <- lm(PC1_hirshorn_child ~ 1, data=pref_data_edup) 
fit1 <- lm (PC1_hirshorn_child ~ pat_edu, data=pref_data_edup)
anova(fit0, fit1, test="Chisq")



##ADHD

#install.packages('xlsx')
library("xlsx")
#install.packages('rJava')
library(rJava)
detach(package:dplyr)
library(dplyr)


# Function
sigfunction <- function(data, col_of_interest, ID){
  
  # Determine ID columns and effect column
  interest_idx <- which(colnames(data) == col_of_interest) #interest_idx is the exposure
  colnames(data)[interest_idx] <- "effect"
  ID_idx <- which(colnames(data) == ID)
  colnames(data)[ID_idx] <- "ABCNO"
  
  
  # sex 
  sex <- glm(effect ~ SEX,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Sex",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Birth weight
  birthw <- glm(effect ~ bw,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Child birth weight",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Gestational age
  gest_age <- glm(effect ~ AGE_DAYS,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Child gestational age",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Maternal smoking in pregnancy
  mat_smoke <- glm(effect ~ smoking_preg,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal smoking in pregnancy",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Maternal alcohol intake in pregnancy 
  mat_alc <- glm(effect ~ alc_preg_YN_2010,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal alcohol intake in pregnancy",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Parity
  mat_parity <- glm(effect ~ as.numeric(parity),family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal parity",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)  
  
  # Maternal antibiotics usage in pregnancy
  mater_ab <- glm(effect ~ AB_preg,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal use of antibiotics in pregnancy",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Maternal preeclampsia in pregnancy 
  mater_pe <- glm(effect ~ PE,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal pre-eclampsia in pregnancy",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Paternal age at child birth
  pater_age <- glm(effect ~ Fatherage,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Paternal age at birth",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Social circumstances
  soc_pc <- glm(effect ~ Socialcircumstances,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Social circumstances at birth",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Maternal BMI
  mat_weight <- glm(effect ~ Mother_BMI,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal pre-pregnancy BMI",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Maternal age
  mater_age <- glm(effect ~ Motherage,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal age at child birth",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Maternal education
  mater_edu <- glm(effect ~ mat_edu,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Maternal education at child birth",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Paternal education
  pater_edu <- glm(effect ~ pat_edu,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Paternal education at child birth",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # Household income
  house_income <- glm(effect ~ income,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "Household income at child birth",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  
  results <- bind_rows(sex, birthw, gest_age, mat_smoke, mat_alc, mat_parity, mater_ab, mater_pe, pater_age, soc_pc, mat_weight, mater_age, mater_edu, pater_edu, house_income)
  
  return(results)
}


conf_adhd <- sigfunction(data=pref_data ,col_of_interest = "adhd", ID="ABCNO")


##Including maternal education, paternal education and income as categorical variables


pref_data_inc <- pref_data %>% filter(!is.na(income) & !is.na(adhd))
pref_data_inc<- pref_data_inc %>% mutate(income = relevel(factor(income), ref="2"))

glm(adhd ~ income, family ="binomial", data = pref_data_inc) %>% tidy(exp=TRUE, conf.int=TRUE)
#Likehood ratio test
fit0 <- glm(adhd ~ 1, data=pref_data_inc, family ="binomial")
fit1 <- glm (adhd ~ income, data=pref_data_inc, family ="binomial")
anova(fit0, fit1, test="Chisq")

pref_data_edu <- pref_data %>% filter(!is.na(mat_edu) & !is.na(adhd))

glm(adhd ~ mat_edu, family ="binomial", data = pref_data_edu) %>% tidy(exp=TRUE, conf.int=TRUE)
#Likehood ratio test
fit0 <- glm(adhd ~ 1, data=pref_data_edu, family ="binomial")
fit1 <- glm (adhd ~ mat_edu, data=pref_data_edu, family ="binomial")
anova(fit0, fit1, test="Chisq")

pref_data_edup <- pref_data %>% filter(!is.na(pat_edu) & !is.na(adhd))
pref_data_edup <- pref_data_edup %>% mutate(pat_edu = relevel(factor(pat_edu), ref="1"))

glm(adhd ~ pat_edu, family ="binomial", data = pref_data_edup) %>% tidy(exp=TRUE, conf.int=TRUE)
#Likehood ratio test
fit0 <- glm(adhd ~ 1, data=pref_data_edup, family ="binomial") 
fit1 <- glm (adhd ~ pat_edu, data=pref_data_edup, family ="binomial")
anova(fit0, fit1, test="Chisq")


######################### Test of models - assortative mating #########################

##Checking for correlation between maternal and paternal PRS

hist(pref_data$PC1_hirshorn_far)
hist(pref_data$PC3_hirshorn_mor)

cor.test(pref_data$PC1_hirshorn_mor, pref_data$PC1_hirshorn_far, method = "pearson")


######################### Analyses - food preference PGS Cole et al ##################


#--------------------------------------------------------------------------------


## to include N
# glm
getIRR <- function(fit){
  data.frame(IRR = exp(coef(summary(fit)))[2,1] %>% round(2),
             CI.lower = exp(confint(fit))[2,1] %>% round(2),
             CI.upper = exp(confint(fit))[2,2] %>% round(2),
             p = coef(summary(fit))[2,4] %>% round(15),
             n = nrow(fit$model),
             beta=coef(summary(fit))[2,1] %>% round(3),
             se=coef(summary(fit))[2,2] %>% round(3),
             outcome = names(fit$model)[1])
}

# linear
getIRR_lm <- function(fit){
  data.frame(beta=coef(summary(fit))[2,1]%>% round(2),
             CI.lower = confint(fit)[2,1] %>% round(3),
             CI.upper = confint(fit)[2,2] %>% round(3),
             p = coef(summary(fit))[2,4] %>% as.numeric(as.character(formatC(., format = "e", digits = 15))),
             n = nrow(fit$model),
             se=coef(summary(fit))[2,2]%>% round(3),
             outcome = names(fit$model)[1])
}


#install.packages('xlsx')
library("xlsx")
#install.packages('rJava')
library(rJava)
detach(package:dplyr)
library(dplyr)


# Function
sigfunction <- function(data, col_of_interest, ID){
  
  # Determine ID columns and effect column
  interest_idx <- which(colnames(data) == col_of_interest) #interest_idx is the exposure
  colnames(data)[interest_idx] <- "effect"
  ID_idx <- which(colnames(data) == ID)
  colnames(data)[ID_idx] <- "ABCNO"
  
  
  # ADHD 
  adhd <- glm(adhd ~ effect + SEX,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD 10 years, yes/no",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD - inattentive
  adhd_in <- glm(adhd_inattention ~ effect + SEX,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD inattentive presentation, 10 years, yes/no",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)

  # ADHD - combined type
  adhd_imh <- glm(adhd_ih ~ effect + SEX,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD combined presentation, 10 years, yes/no",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)  
  
  # ADHD-RS total
  adhdrs_tot <- lm(adhd_rs_copsych_1_18 ~ effect + SEX,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS total score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS inattention
  adhdrs_att <- lm(adhd_rs_copsych_1_9 ~ effect + SEX,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS inattention score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS impulsivity/hyperactivity 
  adhdrs_ih <- lm(adhd_rs_copsych_10_18 ~ effect + SEX,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS impulsivity/hyperactivity score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
 
  
  results <- bind_rows(adhd, adhd_in, adhd_imh, adhdrs_tot, adhdrs_att, adhdrs_ih)
  
  
  return(results)
}

results_pc1_food <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_child", ID="ABCNO")


results_pc1_food_mor <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_mor", ID="ABCNO")


results_pc1_food_far <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_far", ID="ABCNO")



######################### Analyses - adjusted for child PRS and paternal PRS #############################


#install.packages('xlsx')
library("xlsx")
#install.packages('rJava')
library(rJava)
detach(package:dplyr)
library(dplyr)


# Function
sigfunction <- function(data, col_of_interest, ID){
  
  # Determine ID columns and effect column
  interest_idx <- which(colnames(data) == col_of_interest) #interest_idx is the exposure
  colnames(data)[interest_idx] <- "effect"
  ID_idx <- which(colnames(data) == ID)
  colnames(data)[ID_idx] <- "ABCNO"
  
  
  # ADHD 
  adhd <- glm(adhd ~ effect + PC1_hirshorn_child + PC1_hirshorn_far + SEX,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD 10 years, yes/no",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD inattentive type
  adhd_in <- glm(adhd_inattention ~ effect + PC1_hirshorn_child + PC1_hirshorn_far + SEX,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD inattentive presentation, 10 years, yes/no",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD combined type
  adhd_imh <- glm(adhd_ih ~ effect + PC1_hirshorn_child + PC1_hirshorn_far + SEX,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD combined presentation, 10 years, yes/no",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  
  # ADHD-RS total
  adhdrs_tot <- lm(adhd_rs_copsych_1_18 ~ effect + PC1_hirshorn_child + PC1_hirshorn_far + SEX,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS total score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS inattention
  adhdrs_att <- lm(adhd_rs_copsych_1_9 ~ effect + PC1_hirshorn_child + PC1_hirshorn_far + SEX,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS inattention score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS impulsivity/hyperactivity 
  adhdrs_ih <- lm(adhd_rs_copsych_10_18 ~ effect + PC1_hirshorn_child + PC1_hirshorn_far + SEX,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS impulsivity/hyperactivity score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  

  
  results <- bind_rows(adhd, adhd_in, adhd_imh, adhdrs_tot, adhdrs_att, adhdrs_ih)
  
  
  return(results)
}

results_pc1_food_mor <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_mor", ID="ABCNO")


#Test of full confounderadjustment of significant finding in trio model
glm(adhd_rs_copsych_1_18~PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + CRPm + Mother_BMI + pufa + mat_edu + pat_edu + Motherage + income + Fatherage + SEX + AGE_DAYS + bw, data = pref_data) %>% tidy(conf.int=T)
##Only adjusting for variables significantly associated with the PGS 
glm(adhd_rs_copsych_1_18~PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + pat_edu + mat_edu, data = pref_data) %>% tidy(conf.int=T)
glm(adhd_rs_copsych_1_18~PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX, data = pref_data) %>% tidy(conf.int=T)


######################### Analyses - adjusted for child PRS and maternal PRS #############################


#install.packages('xlsx')
library("xlsx")
#install.packages('rJava')
library(rJava)
detach(package:dplyr)
library(dplyr)


# Function
sigfunction <- function(data, col_of_interest, ID){
  
  # Determine ID columns and effect column
  interest_idx <- which(colnames(data) == col_of_interest) #interest_idx is the exposure
  colnames(data)[interest_idx] <- "effect"
  ID_idx <- which(colnames(data) == ID)
  colnames(data)[ID_idx] <- "ABCNO"

  
  
  # ADHD 
  adhd <- glm(adhd ~ effect + PC1_hirshorn_child + PC1_hirshorn_mor + SEX,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD 10 years, yes/no",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD inattentive type
  adhd_in <- glm(adhd_inattention ~ effect + PC1_hirshorn_child + PC1_hirshorn_mor + SEX,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD inattentive presentation, 10 years, yes/no",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD combined type
  adhd_imh <- glm(adhd_ih ~ effect + PC1_hirshorn_child + PC1_hirshorn_mor + SEX,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD combined presentation, 10 years, yes/no",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  
  # ADHD-RS total
  adhdrs_tot <- lm(adhd_rs_copsych_1_18 ~ effect + PC1_hirshorn_child + PC1_hirshorn_mor + SEX,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS total score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS inattention
  adhdrs_att <- lm(adhd_rs_copsych_1_9 ~ effect + PC1_hirshorn_child + PC1_hirshorn_mor + SEX,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS inattention score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS impulsivity/hyperactivity 
  adhdrs_ih <- lm(adhd_rs_copsych_10_18 ~ effect + PC1_hirshorn_child + PC1_hirshorn_mor + SEX,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS impulsivity/hyperactivity score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  
  results <- bind_rows(adhd, adhd_in, adhd_imh, adhdrs_tot, adhdrs_att, adhdrs_ih)
  
  
  return(results)
}



results_pc1_food_far <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_far", ID="ABCNO")




###################### Analyses - adjusted for parental PRS, trio model #######################


#install.packages('xlsx')
library("xlsx")
#install.packages('rJava')
library(rJava)
detach(package:dplyr)
library(dplyr)

# Function
sigfunction <- function(data, col_of_interest, ID){
  
  # Determine ID columns and effect column
  interest_idx <- which(colnames(data) == col_of_interest) #interest_idx is the exposure
  colnames(data)[interest_idx] <- "effect"
  ID_idx <- which(colnames(data) == ID)
  colnames(data)[ID_idx] <- "ABCNO"
  

  
  # ADHD 
  adhd <- glm(adhd ~ effect + PC1_hirshorn_mor + PC1_hirshorn_far + SEX,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD 10 years, yes/no",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD inattentive presentation
  adhd_in <- glm(adhd_inattention ~ effect + PC1_hirshorn_mor + PC1_hirshorn_far + SEX,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD inattentive presentation, 10 years, yes/no",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD combined presentation
  adhd_imh <- glm(adhd_ih ~ effect + PC1_hirshorn_mor + PC1_hirshorn_far + SEX,family = "binomial",data=data) %>%
    getIRR %>%
    mutate(CI =paste0(IRR," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD combined presentation, 10 years, yes/no",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  
  # ADHD-RS total
  adhdrs_tot <- lm(adhd_rs_copsych_1_18 ~ effect + PC1_hirshorn_mor + PC1_hirshorn_far + SEX,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS total score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS inattention
  adhdrs_att <- lm(adhd_rs_copsych_1_9 ~ effect + PC1_hirshorn_mor + PC1_hirshorn_far + SEX,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS inattention score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS impulsivity/hyperactivity 
  adhdrs_ih <- lm(adhd_rs_copsych_10_18 ~ effect + PC1_hirshorn_mor + PC1_hirshorn_far + SEX,data=data) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS impulsivity/hyperactivity score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
 
  results <- bind_rows(adhd, adhd_in, adhd_imh, adhdrs_tot, adhdrs_att, adhdrs_ih)
  
  
  return(results)
}


results_pc1_food_child_adjust_parent <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_child", ID="ABCNO")




######################## Performing analyses excluding ADHD cases #####################


#install.packages('xlsx')
library("xlsx")
#install.packages('rJava')
library(rJava)
detach(package:dplyr)
library(dplyr)


# Function
sigfunction <- function(data, col_of_interest, ID){
  
  # Determine ID columns and effect column
  interest_idx <- which(colnames(data) == col_of_interest) #interest_idx is the exposure
  colnames(data)[interest_idx] <- "effect"
  ID_idx <- which(colnames(data) == ID)
  colnames(data)[ID_idx] <- "ABCNO"
  
 
  # ADHD-RS total
  adhdrs_tot <- lm(adhd_rs_copsych_1_18 ~ effect + SEX,data=subset(data, adhd==0)) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS total score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS inattention
  adhdrs_att <- lm(adhd_rs_copsych_1_9 ~ effect + SEX,data=subset(data, adhd==0)) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS inattention score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS impulsivity/hyperactivity 
  adhdrs_ih <- lm(adhd_rs_copsych_10_18 ~ effect + SEX,data=subset(data, adhd==0)) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS impulsivity/hyperactivity score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  

  
  results <- bind_rows(adhdrs_tot, adhdrs_att, adhdrs_ih)
  
  
  return(results)
}


results_pc1_food_child_diet_noadhd <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_child", ID="ABCNO")


results_pc1_food_mom_diet_noadhd <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_mor", ID="ABCNO")


results_pc1_food_dad_diet_noadhd <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_far", ID="ABCNO")


## Trio child

#install.packages('xlsx')
library("xlsx")
#install.packages('rJava')
library(rJava)
detach(package:dplyr)
library(dplyr)


# Function
sigfunction <- function(data, col_of_interest, ID){
  
  # Determine ID columns and effect column
  interest_idx <- which(colnames(data) == col_of_interest) #interest_idx is the exposure
  colnames(data)[interest_idx] <- "effect"
  ID_idx <- which(colnames(data) == ID)
  colnames(data)[ID_idx] <- "ABCNO"
  
  
  # ADHD-RS total
  adhdrs_tot <- lm(adhd_rs_copsych_1_18 ~ effect + SEX + PC1_hirshorn_mor + PC1_hirshorn_far,data=subset(data, adhd==0)) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS total score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS inattention
  adhdrs_att <- lm(adhd_rs_copsych_1_9 ~ effect + SEX + PC1_hirshorn_mor + PC1_hirshorn_far,data=subset(data, adhd==0)) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS inattention score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS impulsivity/hyperactivity 
  adhdrs_ih <- lm(adhd_rs_copsych_10_18 ~ effect + SEX + PC1_hirshorn_mor + PC1_hirshorn_far,data=subset(data, adhd==0)) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS impulsivity/hyperactivity score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  
  
  results <- bind_rows(adhdrs_tot, adhdrs_att, adhdrs_ih)
  
  
  return(results)
}


results_pc1_food_child_diet_noadhd_trio <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_child", ID="ABCNO")



##Trio maternal

#install.packages('xlsx')
library("xlsx")
#install.packages('rJava')
library(rJava)
detach(package:dplyr)
library(dplyr)


# Function
sigfunction <- function(data, col_of_interest, ID){
  
  # Determine ID columns and effect column
  interest_idx <- which(colnames(data) == col_of_interest) #interest_idx is the exposure
  colnames(data)[interest_idx] <- "effect"
  ID_idx <- which(colnames(data) == ID)
  colnames(data)[ID_idx] <- "ABCNO"
  
  
  # ADHD-RS total
  adhdrs_tot <- lm(adhd_rs_copsych_1_18 ~ effect + SEX + PC1_hirshorn_child + PC1_hirshorn_far,data=subset(data, adhd==0)) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS total score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS inattention
  adhdrs_att <- lm(adhd_rs_copsych_1_9 ~ effect + SEX + PC1_hirshorn_child + PC1_hirshorn_far,data=subset(data, adhd==0)) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS inattention score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS impulsivity/hyperactivity 
  adhdrs_ih <- lm(adhd_rs_copsych_10_18 ~ effect + SEX + PC1_hirshorn_child + PC1_hirshorn_far,data=subset(data, adhd==0)) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS impulsivity/hyperactivity score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  
  
  results <- bind_rows(adhdrs_tot, adhdrs_att, adhdrs_ih)
  
  
  return(results)
}



results_pc1_food_mom_diet_noadhd_trio <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_mor", ID="ABCNO")


##Trio paternal

#install.packages('xlsx')
library("xlsx")
#install.packages('rJava')
library(rJava)
detach(package:dplyr)
library(dplyr)


# Function
sigfunction <- function(data, col_of_interest, ID){
  
  # Determine ID columns and effect column
  interest_idx <- which(colnames(data) == col_of_interest) #interest_idx is the exposure
  colnames(data)[interest_idx] <- "effect"
  ID_idx <- which(colnames(data) == ID)
  colnames(data)[ID_idx] <- "ABCNO"
  
  
  # ADHD-RS total
  adhdrs_tot <- lm(adhd_rs_copsych_1_18 ~ effect + SEX + PC1_hirshorn_child + PC1_hirshorn_mor,data=subset(data, adhd==0)) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS total score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS inattention
  adhdrs_att <- lm(adhd_rs_copsych_1_9 ~ effect + SEX + PC1_hirshorn_child + PC1_hirshorn_mor,data=subset(data, adhd==0)) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS inattention score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  # ADHD-RS impulsivity/hyperactivity 
  adhdrs_ih <- lm(adhd_rs_copsych_10_18 ~ effect + SEX + PC1_hirshorn_child + PC1_hirshorn_mor,data=subset(data, adhd==0)) %>%
    getIRR_lm %>%
    mutate(CI =paste0(beta," [", sprintf("%.2f",round(CI.lower,2)), ";",sprintf("%.2f",round(CI.upper,2)),"]"),
           PVAL = round(p,5),
           outcome = "ADHD-RS impulsivity/hyperactivity score",
           N=n) %>%
    dplyr::select(outcome,N,CI,PVAL)
  
  
  
  results <- bind_rows(adhdrs_tot, adhdrs_att, adhdrs_ih)
  
  
  return(results)
}




results_pc1_food_dad_diet_noadhd_trio <- sigfunction(data=pref_data ,col_of_interest = "PC1_hirshorn_far", ID="ABCNO")



######################## Plotting selected PRS against nutrients and food groups #####################

#### Hirschhorn PC1

## Maternal PCI Hirschhorn against food groups

#dataset = diet

# Create histograms for all variables in one plot

diet_long_hist <- diet %>% select(`Low fat dairy`:`Spices`) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
ggplot(diet_long_hist, aes(x = value)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free") +  # Separate histogram per variable
  labs(
    title = "Histograms for Multiple Variables",
    x = "Value",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14)

scale2 <- function(x, na.rm = FALSE) (x %>% scale())
AllMetabolic_CM <- diet %>% dplyr::mutate(across(c(`Low fat dairy`:`Spices`), scale2)) 
glimpse(AllMetabolic_CM)

AllMetabolic_CM <- left_join(AllMetabolic_CM, pref_data %>% select(ABCNO, PC1_hirshorn_mor), by ="ABCNO")
dataset_hth <- left_join(AllMetabolic_CM, pref_data %>% select(ABCNO, Diet_PC2), by ="ABCNO") ##Data set for a head to head comparison of the diet PC1 vs WDP

getmixedSTAT2 <- function(x){
  
  mdl <- glm(data = x, Metabolic_result ~ Screen_time_value, family= "gaussian")
  res <- mdl %>% tidy(., exp=F, conf.int=T) %>% filter(term != "(Intercept)")
}

summary(AllMetabolic_CM)
output3<- AllMetabolic_CM   %>% 
  gather(Metabolic_factor,Metabolic_result, `Low fat dairy`:Spices) %>% 
  gather(Screen_time,Screen_time_value,PC1_hirshorn_mor) %>% 
  filter(!is.na(Screen_time_value)) %>% 
  filter(!is.infinite(Screen_time_value)) %>% 
  group_by(Metabolic_factor, Screen_time) %>% do(getmixedSTAT2(x = .)) %>% ungroup()

output3$p.valuefdr <- p.adjust(output3$p.value,method="fdr" )

output3<- output3 %>% mutate(FILL=ifelse(estimate >= 0.09, "High", ifelse(estimate < 0.09 & estimate > -0.093, "mid", ifelse(estimate <-0.093, "low", NA)))) %>% mutate(FILL=as.factor(FILL))
output3<- output3 %>% mutate(p.value=ifelse(p.value<0.001,"<0.001", round(p.value,3)))
output3<- output3 %>% mutate(p.valuedfrbi=ifelse(p.valuefdr<0.05,1,0))
output3<-output3[order(output3$estimate),]  

order <- output3$Metabolic_factor


output3$Metabolic_factor<- factor(output3$Metabolic_factor, levels= output3$Metabolic_factor)

saved<-output3$Metabolic_factor


a <- ggplot(output3 %>% filter(FILL != "mid"), aes(y = Metabolic_factor, 
                                                   x = as.numeric(as.character((estimate))), 
                                                   xmin = as.numeric(as.character((conf.low))), 
                                                   xmax = as.numeric(as.character((conf.high))), 
                                                   fill = FILL)) + 
  geom_bar(stat = "identity", position = "dodge", alpha = 1) + 
  geom_errorbarh(height = 0.3) +
  geom_hline(yintercept = seq(-0.3, 0.45, by = 0.1), color = "gray80", linetype = "dotted") +
  xlim(c(-0.3, 0.45)) + 
  scale_fill_manual(values = c("#8BC34A", "#E57373")) + 
  guides(fill = "none") +
  xlab("Scaled Estimate (SD)") +
  ylab(NULL) +
  geom_vline(xintercept = 0, color = "black", size = 0.1) +
  theme_minimal(base_size = 15) +
  ggtitle("Maternal dietary pattern PGS") +
  theme(
    plot.title = element_text(hjust = 0.42, size = 16),
    axis.text.y = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  )
   
a <- a + geom_text(aes(
  label = paste0("p = ", formatC(p.value, format = "e", digits = 2), 
                 ifelse(p.valuefdr < 0.05, " *", "")), 
  x = 0.36), 
  hjust = 0.03, 
  size = 3.5) 

a



## PCA of food groups ##

#Importing total energy intake
tot_energy <- nutrients_df %>% select(ABCNO, `ENERGI: Energy (kJ/day)`)

#Merge tot_energy to diet data
pca_fg <- left_join(AllMetabolic_CM, tot_energy, by = "ABCNO")
pca_fg <- na.omit(pca_fg) %>% select(-PC1_hirshorn_mor)

#Calibrate for total energy before PCA
vars_to_adjust <- names(pca_fg)[2:44]

for (v in vars_to_adjust) {
  pca_fg[[paste0(v, "_adj")]] <- residuals(lm(pca_fg [[v]] ~ pca_fg$`ENERGI: Energy (kJ/day)`)) + mean(pca_fg[[v]], na.rm = TRUE)
}

#Check this works - it does
cor.test(pca_fg$`ENERGI: Energy (kJ/day)`, pca_fg$`Whole grains`, na.rm = TRUE)
cor.test(pca_fg$`ENERGI: Energy (kJ/day)`, pca_fg$`Whole grains_adj`, na.rm = TRUE)

#Remove non adjusted values
pca_fg <- pca_fg[c("ABCNO", grep("_adj$", names(pca_fg), value = TRUE))]

# Step 1: Separate the ABCNO column
ABCNO <- pca_fg$ABCNO

# Step 2: Remove the ABCNO column for PCA analysis
pca_data_fg <- pca_fg[, !names(pca_fg) %in% 'ABCNO']

# Step 3: Perform the PCA
data.pca_fg <- princomp(pca_data_fg)

# Step 4: Summarize the PCA
summary(data.pca_fg)

# Step 5: Keep track of the PCA scores with the corresponding ABCNO
pca_scores_with_IDs <- data.frame(ABCNO = ABCNO, data.pca_fg$scores)
pca_scores_total_dataset <- left_join(pca_scores_with_IDs, pref_data, by = "ABCNO")

# View the PCA scores with the corresponding IDs
head(pca_scores_with_IDs)

data.pca_fg$loadings[, 1:10]
fviz_eig(data.pca_fg, addlabels = TRUE)
fviz_cos2(data.pca_fg, choice = "var", axes = 1)
fviz_cos2(data.pca_fg, choice = "var", axes = 2)

fviz_pca_var(data.pca_fg, col.var = "cos2", axes = c(1,2),
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

cor.test(pca_scores_total_dataset$Comp.1, pca_scores_total_dataset$PC1_hirshorn_mor)
cor.test(pca_scores_total_dataset$Comp.1, pca_scores_total_dataset$PC1_hirshorn_mor)$estimate^2 #R squared
cor.test(pca_scores_total_dataset$Comp.2, pca_scores_total_dataset$PC1_hirshorn_mor)

fviz_pca_var(data.pca_fg, col.var = "cos2", axes = c(2,3),
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

fviz_pca_var(data.pca_fg, col.var = "cos2", axes = c(3,4),
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)



#######################################################################
### Plotting dietPC2 and Diet PGS against food groups side by side ####
#######################################################################

results <- data.frame(
  food_group = character(),
  exposure = character(),
  estimate = numeric(),
  conf_low = numeric(),
  conf_high = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

food_groups <- names(dataset_hth2)[!names(dataset_hth2) %in% c("ABCNO", "PC1_hirshorn_mor", "Diet_PC2")]

for (fg in food_groups) {
  for (exp in c("PC1_hirshorn_mor", "Diet_PC2")) {
    formula <- as.formula(paste0("`", fg, "` ~ `", exp, "`"))
    model <- lm(formula, data = dataset_hth2)
    
    coef_summary <- summary(model)$coefficients
    conf_int <- confint(model)
    
    estimate <- coef_summary[2, "Estimate"]
    p_value <- coef_summary[2, "Pr(>|t|)"]
    ci_low <- conf_int[2, 1]
    ci_high <- conf_int[2, 2]
    
    results <- rbind(results, data.frame(
      food_group = fg,
      exposure = exp,
      estimate = estimate,
      conf_low = ci_low,
      conf_high = ci_high,
      p_value = p_value
    ))
  }
}

head(results)

#Test of results
lm(`Low fat dairy`~PC1_hirshorn_mor, data = dataset_hth2) %>% tidy(conf.int=T)
lm(`Low fat dairy`~Diet_PC2, data = dataset_hth2) %>% tidy(conf.int=T)


library(dplyr)
library(ggplot2)
library(patchwork)

top20_diet <- results %>%
  filter(exposure == "Diet_PC2") %>%
  mutate(abs_estimate = abs(estimate)) %>%
  arrange(desc(abs_estimate)) %>%
  slice_head(n = 20)

diet_order <- top20_diet %>%
  arrange(estimate) %>%
  pull(food_group)

diet_plot_data <- top20_diet %>%
  mutate(
    food_group = factor(food_group, levels = diet_order),
    exposure = "Western dietary pattern PC",
    sign = ifelse(estimate >= 0, "positive", "negative")
  )

pc1_plot_data <- results %>%
  filter(exposure == "PC1_hirshorn_mor", food_group %in% diet_order) %>%
  mutate(
    estimate = estimate, 
    conf_low = conf_low, 
    conf_high = conf_high, 
    food_group = factor(food_group, levels = diet_order),
    exposure = "Maternal dietary pattern PGS",
    sign = ifelse(estimate >= 0, "positive", "negative")
  )

plot_diet <- ggplot(diet_plot_data, aes(x = food_group, y = estimate, fill = sign)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.25) +
  scale_fill_manual(values = c("positive" = "#FF7F0E", "negative" = "#1F77B4")) +
  coord_flip() +
  labs(
    title = "Maternal Western dietary pattern PC",
    x = "Food group",
    y = "Beta estimate (95% CI)",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 13) + theme(legend.position = "none")

plot_pc1 <- ggplot(pc1_plot_data, aes(x = food_group, y = estimate, fill = sign)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.25) +
  scale_fill_manual(values = c("positive" = "#FF7F0E", "negative" = "#1F77B4")) +
  coord_flip() +
  labs(
    title = "Maternal dietary pattern PGS",
    x = "Food group",
    y = "Beta estimate (95% CI)",
    fill = "Direction"
  ) + theme_minimal(base_size = 13) + theme(axis.text.y = element_blank(),axis.title.y = element_blank())

p_hth <- plot_diet + plot_pc1
p_hth


## Top 20 strongest PC1 associations 

top20_pc1 <- results %>%
  filter(exposure == "PC1_hirshorn_mor") %>%
  mutate(abs_estimate = abs(estimate)) %>%
  arrange(desc(abs_estimate)) %>%
  slice_head(n = 20) %>%
  mutate(
    estimate = estimate, 
    conf_low = conf_low, 
    conf_high = conf_high, 
    sign = ifelse(estimate >= 0, "positive", "negative"),
    exposure = "Reversed maternal dietary pattern PGS"
  )

top20_pc1 <- top20_pc1 %>%
  arrange(estimate) %>%  
  mutate(food_group = factor(food_group, levels = food_group))

diet_plot_data <- results %>%
  filter(exposure == "Diet_PC2", food_group %in% top20_pc1$food_group) %>%
  mutate(
    food_group = factor(food_group, levels = levels(top20_pc1$food_group)),  
    exposure = "Maternal Western dietary pattern PC",
    sign = ifelse(estimate >= 0, "positive", "negative")
  )

plot_pc1top <- ggplot(top20_pc1, aes(x = food_group, y = estimate, fill = sign)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.25) +
  scale_fill_manual(values = c("positive" = "#FF7F0E", "negative" = "#1F77B4")) +
  coord_flip() +
  labs(
    title = "Maternal dietary pattern PGS",
    x = "Food group",
    y = "Beta estimate (95% CI)",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 13) + theme(legend.position = "none")

plot_dietpctop <- ggplot(diet_plot_data, aes(x = food_group, y = estimate, fill = sign)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.25) +
  scale_fill_manual(values = c("positive" = "#FF7F0E", "negative" = "#1F77B4")) +
  coord_flip() +
  labs(
    title = "Maternal Western dietary pattern PC",
    x = "Food group",
    y = "Beta estimate (95% CI)",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 13) + theme(axis.text.y = element_blank(),axis.title.y = element_blank())

p_hthpc1top <- plot_pc1top + plot_dietpctop
p_hthpc1top



############################################
## Plot mean dietPC2 for intervals of PC1 ##
############################################

library(dplyr)
library(ggplot2)

pref_binned <- pref_data %>%
  mutate(
    pc1_bin = cut(
      PC1_hirshorn_mor,
      breaks = c(-Inf, -4, -2, -1, 0, 1, 2, Inf),
      labels = c("<-4", "-4,-2", "-2,-1", "-1,0", "0,1", "1,2", "+2>"),
      include.lowest = TRUE
    )
  )


summary_binned <- pref_binned %>%
  group_by(pc1_bin) %>%
  summarise(
    mean_diet_pc2 = mean(Diet_PC2, na.rm = TRUE),
    sd_diet_pc2 = sd(Diet_PC2, na.rm = TRUE),
    n = sum(!is.na(Diet_PC2)),
    se = sd_diet_pc2 / sqrt(n),
    ci_low = mean_diet_pc2 - 1.96 * se,
    ci_high = mean_diet_pc2 + 1.96 * se,
    .groups = "drop"
  ) %>%
  filter(!is.na(pc1_bin))

ggplot(summary_binned, aes(x = pc1_bin, y = mean_diet_pc2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, color = "steelblue") +
  ylim(-1.5, 1.5) +
  labs(
    title = "Mean COPSAC maternal Western dietary pattern PC by dietary pattern PGS",
    x = "Maternal dietary pattern PGS, intervals of 1 SD",
    y = "Mean maternal Western dietary pattern PC (± 95% CI)"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



################################## Figures - forest plot trio model ###########################


##Diagnoses

##Child ADHD

adhd_c <- glm(adhd ~ PC1_hirshorn_child + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_c <- adhd_c %>% subset(term =="PC1_hirshorn_child")
adhd_c <- adhd_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD")

adhd_c_trio <- glm(adhd ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_c_trio <- adhd_c_trio %>% subset(term =="PC1_hirshorn_child")
adhd_c_trio <- adhd_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD")

##Mother ADHD

adhd_m <- glm(adhd ~ PC1_hirshorn_mor + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_m <- adhd_m %>% subset(term =="PC1_hirshorn_mor")
adhd_m <- adhd_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD")

adhd_m_trio <- glm(adhd ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_m_trio <- adhd_m_trio %>% subset(term =="PC1_hirshorn_mor")
adhd_m_trio <- adhd_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD")

##Father ADHD

adhd_f <- glm(adhd ~ PC1_hirshorn_far + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_f <- adhd_f %>% subset(term =="PC1_hirshorn_far")
adhd_f <- adhd_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD")

adhd_f_trio <- glm(adhd ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_f_trio <- adhd_f_trio %>% subset(term =="PC1_hirshorn_far")
adhd_f_trio <- adhd_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD")

##Child ADD

add_c <- glm(adhd_inattention ~ PC1_hirshorn_child + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_c <- add_c %>% subset(term =="PC1_hirshorn_child")
add_c <- add_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_c_trio <- glm(adhd_inattention ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_c_trio <- add_c_trio %>% subset(term =="PC1_hirshorn_child")
add_c_trio <- add_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Mother ADD

add_m <- glm(adhd_inattention ~ PC1_hirshorn_mor + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_m <- add_m %>% subset(term =="PC1_hirshorn_mor")
add_m <- add_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_m_trio <- glm(adhd_inattention ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_m_trio <- add_m_trio %>% subset(term =="PC1_hirshorn_mor")
add_m_trio <- add_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Father ADD

add_f <- glm(adhd_inattention ~ PC1_hirshorn_far + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_f <- add_f %>% subset(term =="PC1_hirshorn_far")
add_f <- add_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_f_trio <- glm(adhd_inattention ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_f_trio <- add_f_trio %>% subset(term =="PC1_hirshorn_far")
add_f_trio <- add_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")
 
##Child ADHD - Combined Presentation

adhd1_c <- glm(adhd_ih ~ PC1_hirshorn_child + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_c <- adhd1_c %>% subset(term =="PC1_hirshorn_child")
adhd1_c <- adhd1_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_c_trio <- glm(adhd_ih ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_c_trio <- adhd1_c_trio %>% subset(term =="PC1_hirshorn_child")
adhd1_c_trio <- adhd1_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

##Mother ADHD - Combined Presentation

adhd1_m <- glm(adhd_ih ~ PC1_hirshorn_mor + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_m <- adhd1_m %>% subset(term =="PC1_hirshorn_mor")
adhd1_m <- adhd1_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_m_trio <- glm(adhd_ih ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_m_trio <- adhd1_m_trio %>% subset(term =="PC1_hirshorn_mor")
adhd1_m_trio <- adhd1_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

##Father ADHD - Combined Presentation

adhd1_f <- glm(adhd_ih ~ PC1_hirshorn_far + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_f <- adhd1_f %>% subset(term =="PC1_hirshorn_far")
adhd1_f <- adhd1_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_f_trio <- glm(adhd_ih ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_f_trio <- adhd1_f_trio %>% subset(term =="PC1_hirshorn_far")
adhd1_f_trio <- adhd1_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

forestadhd <- rbind(adhd_c, adhd_c_trio, adhd_m, adhd_m_trio, adhd_f, adhd_f_trio, add_c, add_c_trio, add_m, add_m_trio, add_f, add_f_trio, adhd1_c, adhd1_c_trio, adhd1_m, adhd1_m_trio, adhd1_f, adhd1_f_trio)
forestadhd$Condition <- factor(forestadhd$Condition, levels=c("norm", "trio"))
forestadhd$Group <- factor(forestadhd$Group, levels=c("Child direct effect", "Child PGS", "Paternal indirect effect", "Paternal PGS", "Maternal indirect effect", "Maternal PGS"))
forestadhd$Copsych <- factor(forestadhd$Copsych, levels=c("ADHD", "ADHD  Inattentive Presentation", "ADHD  Combined Presentation"))
forestadhd <- forestadhd %>% mutate(significant = ifelse(p.value<0.05,1,0))

theme_set(theme_bw())
redPalette <- c("#E69F00", "#009E73")

p_adhd_trio = ggplot(data=forestadhd, aes(x = Group,y = estimate, ymin = conf.low, ymax = conf.high, na.rm = TRUE))+
  geom_point(size = 2, shape = 19, aes(col=Condition)) +
  geom_hline(aes(fill=Group),yintercept =1, linetype=2)+
  xlab("") + ylab("Odds Ratio (95% Confidence Interval)") +
  coord_flip() + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high,col=Condition),width=0.2,cex=0.5)+ 
  scale_y_continuous(limits = c(0.3,1.8))+ 
  facet_wrap(~Copsych,strip.position="top",nrow=1,scales = "free_x", labeller = label_wrap_gen(width=8) ) +
  scale_colour_manual(values=redPalette) +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_text(size = 14, face="bold"),
        axis.text.x=element_text(size = 16, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 16, face="bold", margin = margin(0.6,0,0.6,0, "cm")), legend.position="None", strip.background = element_rect(fill = "white")) +
  geom_text(aes(label = ifelse(significant == 1, "*", ""), col=Condition), vjust = -0.1, size = 6) 


## ADHD RS total score, Inattention and IH traits

##Total score
##Child ADHD

adhd_total_symp_c <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_child + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_c <- adhd_total_symp_c %>% subset(term =="PC1_hirshorn_child")
adhd_total_symp_c <- adhd_total_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS Total traits")

adhd_total_symp_c_trio <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_c_trio <- adhd_total_symp_c_trio %>% subset(term =="PC1_hirshorn_child")
adhd_total_symp_c_trio <- adhd_total_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model", Copsych = "ADHD-RS Total traits")

##Mother ADHD

adhd_total_symp_m <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_mor + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_m <- adhd_total_symp_m %>% subset(term =="PC1_hirshorn_mor")
adhd_total_symp_m <- adhd_total_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS Total traits")

adhd_total_symp_m_trio <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_m_trio <- adhd_total_symp_m_trio %>% subset(term =="PC1_hirshorn_mor")
adhd_total_symp_m_trio <- adhd_total_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model", Copsych = "ADHD-RS Total traits")

##Father ADHD

adhd_total_symp_f <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_far + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_f <- adhd_total_symp_f %>% subset(term =="PC1_hirshorn_far")
adhd_total_symp_f <- adhd_total_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS Total traits")

adhd_total_symp_f_trio <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_f_trio <- adhd_total_symp_f_trio %>% subset(term =="PC1_hirshorn_far")
adhd_total_symp_f_trio <- adhd_total_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model", Copsych = "ADHD-RS Total traits")

##Inattention
##Child ADHD

adhd_att_symp_c <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_child + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_c <- adhd_att_symp_c %>% subset(term =="PC1_hirshorn_child")
adhd_att_symp_c <- adhd_att_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS Inattention traits")

adhd_att_symp_c_trio <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_c_trio <- adhd_att_symp_c_trio %>% subset(term =="PC1_hirshorn_child")
adhd_att_symp_c_trio <- adhd_att_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model", Copsych = "ADHD-RS Inattention traits")

##Mother ADHD

adhd_att_symp_m <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_mor + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_m <- adhd_att_symp_m %>% subset(term =="PC1_hirshorn_mor")
adhd_att_symp_m <- adhd_att_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS Inattention traits")

adhd_att_symp_m_trio <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_m_trio <- adhd_att_symp_m_trio %>% subset(term =="PC1_hirshorn_mor")
adhd_att_symp_m_trio <- adhd_att_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model", Copsych = "ADHD-RS Inattention traits")

##Father ADHD

adhd_att_symp_f <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_far + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_f <- adhd_att_symp_f %>% subset(term =="PC1_hirshorn_far")
adhd_att_symp_f <- adhd_att_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS Inattention traits")

adhd_att_symp_f_trio <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_f_trio <- adhd_att_symp_f_trio %>% subset(term =="PC1_hirshorn_far")
adhd_att_symp_f_trio <- adhd_att_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model", Copsych = "ADHD-RS Inattention traits")

##Impulsivity/Hyperactivity
##Child ADHD

adhd_ih_symp_c <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_child + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_c <- adhd_ih_symp_c %>% subset(term =="PC1_hirshorn_child")
adhd_ih_symp_c <- adhd_ih_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

adhd_ih_symp_c_trio <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_c_trio <- adhd_ih_symp_c_trio %>% subset(term =="PC1_hirshorn_child")
adhd_ih_symp_c_trio <- adhd_ih_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

##Mother ADHD

adhd_ih_symp_m <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_mor + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_m <- adhd_ih_symp_m %>% subset(term =="PC1_hirshorn_mor")
adhd_ih_symp_m <- adhd_ih_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

adhd_ih_symp_m_trio <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_m_trio <- adhd_ih_symp_m_trio %>% subset(term =="PC1_hirshorn_mor")
adhd_ih_symp_m_trio <- adhd_ih_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

##Father ADHD

adhd_ih_symp_f <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_far + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_f <- adhd_ih_symp_f %>% subset(term =="PC1_hirshorn_far")
adhd_ih_symp_f <- adhd_ih_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

adhd_ih_symp_f_trio <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_f_trio <- adhd_ih_symp_f_trio %>% subset(term =="PC1_hirshorn_far")
adhd_ih_symp_f_trio <- adhd_ih_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")


##Forest plot ADHD RS traits 

forestadhdrssub <- rbind(adhd_total_symp_c, adhd_total_symp_c_trio, adhd_total_symp_m, adhd_total_symp_m_trio, adhd_total_symp_f, adhd_total_symp_f_trio, adhd_att_symp_c, adhd_att_symp_c_trio, adhd_att_symp_m, adhd_att_symp_m_trio, adhd_att_symp_f, adhd_att_symp_f_trio, adhd_ih_symp_c, adhd_ih_symp_c_trio, adhd_ih_symp_m, adhd_ih_symp_m_trio, adhd_ih_symp_f, adhd_ih_symp_f_trio)

forestadhdrssub$Condition <- factor(forestadhdrssub$Condition, levels=c("Adjusted for child sex only", "Trio model"))
forestadhdrssub$Group <- factor(forestadhdrssub$Group, levels=c("Child direct effect", "Child PGS", "Paternal indirect effect", "Paternal PGS", "Maternal indirect effect", "Maternal PGS"))
forestadhdrssub$Copsych <- factor(forestadhdrssub$Copsych, levels=c("ADHD-RS Total traits","ADHD-RS Inattention traits","ADHD-RS  Hyperactivity/Impulsivity traits"))
forestadhdrssub <- forestadhdrssub %>% mutate(significant = ifelse(p.value<0.05,1,0))

theme_set(theme_bw())
redPalette <- c("#E69F00", "#009E73")

p_adhd_rs_sub_trio = ggplot(data=forestadhdrssub, aes(x = Group,y = estimate, ymin = conf.low, ymax = conf.high, na.rm = TRUE))+
  geom_point(size = 2, shape = 19, aes(col=Condition)) +
  geom_hline(aes(fill=Group),yintercept =0, linetype=2)+
  xlab("") + ylab("Beta Estimate (95% Confidence Interval)") +
  coord_flip() + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high,col=Condition),width=0.2,cex=0.5)+ 
  scale_y_continuous(limits = c(-2.1,1))+  # Set x-axis limits
  facet_wrap(~Copsych,strip.position="top",nrow=1,scales = "free_x", labeller = label_wrap_gen(width=15) ) +
  scale_colour_manual(values=redPalette) +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_text(size = 14, face="bold"),
        axis.text.x=element_text(size = 16, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 16, face="bold", margin = margin(0.6,0,0.6,0, "cm")), strip.background = element_rect(fill = "white"),legend.title = element_blank(), legend.text = element_text(size = 14)) +
  geom_text(aes(label = ifelse(significant == 1, "*", ""), col=Condition), vjust = -0.1, size = 6, show.legend=FALSE) 

##Final plot

library(patchwork)
p_fin <- (p_adhd_trio / p_adhd_rs_sub_trio)



################ Forest plot excluding individuals with ADHD ###################


## ADHD RS total score, Inattention and IH traits

##Total score
##Child ADHD

adhd_total_symp_c <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_child + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_c <- adhd_total_symp_c %>% subset(term =="PC1_hirshorn_child")
adhd_total_symp_c <- adhd_total_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS Total traits")

adhd_total_symp_c_trio <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_c_trio <- adhd_total_symp_c_trio %>% subset(term =="PC1_hirshorn_child")
adhd_total_symp_c_trio <- adhd_total_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model", Copsych = "ADHD-RS Total traits")

##Mother ADHD

adhd_total_symp_m <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_mor + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_m <- adhd_total_symp_m %>% subset(term =="PC1_hirshorn_mor")
adhd_total_symp_m <- adhd_total_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS Total traits")

adhd_total_symp_m_trio <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_m_trio <- adhd_total_symp_m_trio %>% subset(term =="PC1_hirshorn_mor")
adhd_total_symp_m_trio <- adhd_total_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model", Copsych = "ADHD-RS Total traits")

##Father ADHD

adhd_total_symp_f <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_far + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_f <- adhd_total_symp_f %>% subset(term =="PC1_hirshorn_far")
adhd_total_symp_f <- adhd_total_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS Total traits")

adhd_total_symp_f_trio <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_f_trio <- adhd_total_symp_f_trio %>% subset(term =="PC1_hirshorn_far")
adhd_total_symp_f_trio <- adhd_total_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model", Copsych = "ADHD-RS Total traits")

##Inattention 
##Child ADHD

adhd_att_symp_c <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_child + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_c <- adhd_att_symp_c %>% subset(term =="PC1_hirshorn_child")
adhd_att_symp_c <- adhd_att_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS Inattention traits")

adhd_att_symp_c_trio <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_c_trio <- adhd_att_symp_c_trio %>% subset(term =="PC1_hirshorn_child")
adhd_att_symp_c_trio <- adhd_att_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model", Copsych = "ADHD-RS Inattention traits")

##Mother ADHD

adhd_att_symp_m <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_mor + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_m <- adhd_att_symp_m %>% subset(term =="PC1_hirshorn_mor")
adhd_att_symp_m <- adhd_att_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS Inattention traits")

adhd_att_symp_m_trio <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_m_trio <- adhd_att_symp_m_trio %>% subset(term =="PC1_hirshorn_mor")
adhd_att_symp_m_trio <- adhd_att_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model", Copsych = "ADHD-RS Inattention traits")

##Father ADHD

adhd_att_symp_f <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_far + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_f <- adhd_att_symp_f %>% subset(term =="PC1_hirshorn_far")
adhd_att_symp_f <- adhd_att_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS Inattention traits")

adhd_att_symp_f_trio <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_f_trio <- adhd_att_symp_f_trio %>% subset(term =="PC1_hirshorn_far")
adhd_att_symp_f_trio <- adhd_att_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model", Copsych = "ADHD-RS Inattention traits")

##Impulsivity/Hyperactivity
##Child ADHD

adhd_ih_symp_c <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_child + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_c <- adhd_ih_symp_c %>% subset(term =="PC1_hirshorn_child")
adhd_ih_symp_c <- adhd_ih_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

adhd_ih_symp_c_trio <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_c_trio <- adhd_ih_symp_c_trio %>% subset(term =="PC1_hirshorn_child")
adhd_ih_symp_c_trio <- adhd_ih_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

##Mother ADHD

adhd_ih_symp_m <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_mor + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_m <- adhd_ih_symp_m %>% subset(term =="PC1_hirshorn_mor")
adhd_ih_symp_m <- adhd_ih_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

adhd_ih_symp_m_trio <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_m_trio <- adhd_ih_symp_m_trio %>% subset(term =="PC1_hirshorn_mor")
adhd_ih_symp_m_trio <- adhd_ih_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

##Father ADHD

adhd_ih_symp_f <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_far + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_f <- adhd_ih_symp_f %>% subset(term =="PC1_hirshorn_far")
adhd_ih_symp_f <- adhd_ih_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for child sex only", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

adhd_ih_symp_f_trio <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX, data = filter(pref_data, adhd==0)) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_f_trio <- adhd_ih_symp_f_trio %>% subset(term =="PC1_hirshorn_far")
adhd_ih_symp_f_trio <- adhd_ih_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")


##Forest plot ADHD RS traits 

forestadhdrssub <- rbind(adhd_total_symp_c, adhd_total_symp_c_trio, adhd_total_symp_m, adhd_total_symp_m_trio, adhd_total_symp_f, adhd_total_symp_f_trio, adhd_att_symp_c, adhd_att_symp_c_trio, adhd_att_symp_m, adhd_att_symp_m_trio, adhd_att_symp_f, adhd_att_symp_f_trio, adhd_ih_symp_c, adhd_ih_symp_c_trio, adhd_ih_symp_m, adhd_ih_symp_m_trio, adhd_ih_symp_f, adhd_ih_symp_f_trio)

forestadhdrssub$Condition <- factor(forestadhdrssub$Condition, levels=c("Adjusted for child sex only", "Trio model"))
forestadhdrssub$Group <- factor(forestadhdrssub$Group, levels=c("Child direct effect", "Child PGS", "Paternal indirect effect", "Paternal PGS", "Maternal indirect effect", "Maternal PGS"))
forestadhdrssub$Copsych <- factor(forestadhdrssub$Copsych, levels=c("ADHD-RS Total traits","ADHD-RS Inattention traits","ADHD-RS  Hyperactivity/Impulsivity traits"))
forestadhdrssub <- forestadhdrssub %>% mutate(significant = ifelse(p.value<0.05,1,0))

theme_set(theme_bw())
redPalette <- c("#E69F00", "#009E73")

p_adhd_rs_sub_trio_no_adhd = ggplot(data=forestadhdrssub, aes(x = Group,y = estimate, ymin = conf.low, ymax = conf.high, na.rm = TRUE))+
  geom_point(size = 2, shape = 19, aes(col=Condition)) +
  geom_hline(aes(fill=Group),yintercept =0, linetype=2)+
  xlab("") + ylab("Beta Estimate (95% Confidence Interval)") +
  coord_flip() + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high,col=Condition),width=0.2,cex=0.5)+ 
  scale_y_continuous(limits = c(-2,1.5))+  # Set x-axis limits
  facet_wrap(~Copsych,strip.position="top",nrow=1,scales = "free_x", labeller = label_wrap_gen(width=15) ) +
  scale_colour_manual(values=redPalette) +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_text(size = 14, face="bold"),
        axis.text.x=element_text(size = 16, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 16, face="bold", margin = margin(0.6,0,0.6,0, "cm")), strip.background = element_rect(fill = "white"),legend.title = element_blank(), legend.text = element_text(size = 14)) +
  geom_text(aes(label = ifelse(significant == 1, "*", ""), col=Condition), vjust = -0.1, size = 6, show.legend=FALSE) 


############################ Figures - forest plot, maternal ADHD PRS, adjusted for all ADHD PRS's ########################


##Diagnoses

##Child ADHD

adhd_c <- glm(adhd ~ PC1_hirshorn_child + SEX + adhd_prs, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_c <- adhd_c %>% subset(term =="PC1_hirshorn_child")
adhd_c <- adhd_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD")

adhd_c_trio <- glm(adhd ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_c_trio <- adhd_c_trio %>% subset(term =="PC1_hirshorn_child")
adhd_c_trio <- adhd_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD")

##Mother ADHD

adhd_m <- glm(adhd ~ PC1_hirshorn_mor + SEX + adhd_prs_mor, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_m <- adhd_m %>% subset(term =="PC1_hirshorn_mor")
adhd_m <- adhd_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD")

adhd_m_trio <- glm(adhd ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_m_trio <- adhd_m_trio %>% subset(term =="PC1_hirshorn_mor")
adhd_m_trio <- adhd_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD")

##Father ADHD

adhd_f <- glm(adhd ~ PC1_hirshorn_far + SEX + adhd_prs_far, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_f <- adhd_f %>% subset(term =="PC1_hirshorn_far")
adhd_f <- adhd_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD")

adhd_f_trio <- glm(adhd ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_f_trio <- adhd_f_trio %>% subset(term =="PC1_hirshorn_far")
adhd_f_trio <- adhd_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD")

##Child ADD

add_c <- glm(adhd_inattention ~ PC1_hirshorn_child + SEX + adhd_prs, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_c <- add_c %>% subset(term =="PC1_hirshorn_child")
add_c <- add_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_c_trio <- glm(adhd_inattention ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_c_trio <- add_c_trio %>% subset(term =="PC1_hirshorn_child")
add_c_trio <- add_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Mother ADD

add_m <- glm(adhd_inattention ~ PC1_hirshorn_mor + SEX + adhd_prs_mor, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_m <- add_m %>% subset(term =="PC1_hirshorn_mor")
add_m <- add_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_m_trio <- glm(adhd_inattention ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_m_trio <- add_m_trio %>% subset(term =="PC1_hirshorn_mor")
add_m_trio <- add_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Father ADD

add_f <- glm(adhd_inattention ~ PC1_hirshorn_far + SEX + adhd_prs_far, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_f <- add_f %>% subset(term =="PC1_hirshorn_far")
add_f <- add_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_f_trio <- glm(adhd_inattention ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_f_trio <- add_f_trio %>% subset(term =="PC1_hirshorn_far")
add_f_trio <- add_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Child ADHD - Combined Presentation

adhd1_c <- glm(adhd_ih ~ PC1_hirshorn_child + SEX + adhd_prs, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_c <- adhd1_c %>% subset(term =="PC1_hirshorn_child")
adhd1_c <- adhd1_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_c_trio <- glm(adhd_ih ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_c_trio <- adhd1_c_trio %>% subset(term =="PC1_hirshorn_child")
adhd1_c_trio <- adhd1_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

##Mother ADHD - Combined Presentation

adhd1_m <- glm(adhd_ih ~ PC1_hirshorn_mor + SEX + adhd_prs_mor, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_m <- adhd1_m %>% subset(term =="PC1_hirshorn_mor")
adhd1_m <- adhd1_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_m_trio <- glm(adhd_ih ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_m_trio <- adhd1_m_trio %>% subset(term =="PC1_hirshorn_mor")
adhd1_m_trio <- adhd1_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

##Father ADHD - Combined Presentation

adhd1_f <- glm(adhd_ih ~ PC1_hirshorn_far + SEX + adhd_prs_far, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_f <- adhd1_f %>% subset(term =="PC1_hirshorn_far")
adhd1_f <- adhd1_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_f_trio <- glm(adhd_ih ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data, family = binomial) %>% tidy(conf.int = TRUE, exp = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_f_trio <- adhd1_f_trio %>% subset(term =="PC1_hirshorn_far")
adhd1_f_trio <- adhd1_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

forestadhd <- rbind(adhd_c, adhd_c_trio, adhd_m, adhd_m_trio, adhd_f, adhd_f_trio, add_c, add_c_trio, add_m, add_m_trio, add_f, add_f_trio, adhd1_c, adhd1_c_trio, adhd1_m, adhd1_m_trio, adhd1_f, adhd1_f_trio)
forestadhd$Condition <- factor(forestadhd$Condition, levels=c("norm", "trio"))
forestadhd$Group <- factor(forestadhd$Group, levels=c("Child direct effect", "Child PGS", "Paternal indirect effect", "Paternal PGS", "Maternal indirect effect", "Maternal PGS"))
forestadhd$Copsych <- factor(forestadhd$Copsych, levels=c("ADHD", "ADHD  Inattentive Presentation", "ADHD  Combined Presentation"))
forestadhd <- forestadhd %>% mutate(significant = ifelse(p.value<0.05,1,0))

theme_set(theme_bw())
redPalette <- c("#E69F00", "#009E73")

p_adhd_trio_adhd_prs = ggplot(data=forestadhd, aes(x = Group,y = estimate, ymin = conf.low, ymax = conf.high, na.rm = TRUE))+
  geom_point(size = 2, shape = 19, aes(col=Condition)) +
  geom_hline(aes(fill=Group),yintercept =1, linetype=2)+
  xlab("") + ylab("Odds Ratio (95% Confidence Interval)") +
  coord_flip() + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high,col=Condition),width=0.2,cex=0.5)+ 
  scale_y_continuous(limits = c(0.3,1.9))+ 
  facet_wrap(~Copsych,strip.position="top",nrow=1,scales = "free_x", labeller = label_wrap_gen(width=8) ) +
  scale_colour_manual(values=redPalette) +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_text(size = 14, face="bold"),
        axis.text.x=element_text(size = 16, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 16, face="bold", margin = margin(0.6,0,0.6,0, "cm")), legend.position="None", strip.background = element_rect(fill = "white")) +
  geom_text(aes(label = ifelse(significant == 1, "*", ""), col=Condition), vjust = -0.1, size = 6, show.legend=FALSE) 



## ADHD RS total score, Inattention and IH traits

##Total score
##Child ADHD

adhd_total_symp_c <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_child + SEX + adhd_prs, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_c <- adhd_total_symp_c %>% subset(term =="PC1_hirshorn_child")
adhd_total_symp_c <- adhd_total_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for child sex and either\nmaternal, paternal, or child ADHD PRS", Copsych = "ADHD-RS Total traits")

adhd_total_symp_c_trio <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_c_trio <- adhd_total_symp_c_trio %>% subset(term =="PC1_hirshorn_child")
adhd_total_symp_c_trio <- adhd_total_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model with adjustment for\nmaternal, paternal, and child ADHD PRS", Copsych = "ADHD-RS Total traits")

##Mother ADHD

adhd_total_symp_m <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_mor + SEX + adhd_prs_mor, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_m <- adhd_total_symp_m %>% subset(term =="PC1_hirshorn_mor")
adhd_total_symp_m <- adhd_total_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for child sex and either\nmaternal, paternal, or child ADHD PRS", Copsych = "ADHD-RS Total traits")

adhd_total_symp_m_trio <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_m_trio <- adhd_total_symp_m_trio %>% subset(term =="PC1_hirshorn_mor")
adhd_total_symp_m_trio <- adhd_total_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model with adjustment for\nmaternal, paternal, and child ADHD PRS", Copsych = "ADHD-RS Total traits")

##Father ADHD

adhd_total_symp_f <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_far + SEX + adhd_prs_far, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_f <- adhd_total_symp_f %>% subset(term =="PC1_hirshorn_far")
adhd_total_symp_f <- adhd_total_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for child sex and either\nmaternal, paternal, or child ADHD PRS", Copsych = "ADHD-RS Total traits")

adhd_total_symp_f_trio <- glm(adhd_rs_copsych_1_18 ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_f_trio <- adhd_total_symp_f_trio %>% subset(term =="PC1_hirshorn_far")
adhd_total_symp_f_trio <- adhd_total_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model with adjustment for\nmaternal, paternal, and child ADHD PRS", Copsych = "ADHD-RS Total traits")

##Inattention
##Child ADHD

adhd_att_symp_c <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_child + SEX + adhd_prs, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_c <- adhd_att_symp_c %>% subset(term =="PC1_hirshorn_child")
adhd_att_symp_c <- adhd_att_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for child sex and either\nmaternal, paternal, or child ADHD PRS", Copsych = "ADHD-RS Inattention traits")

adhd_att_symp_c_trio <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_c_trio <- adhd_att_symp_c_trio %>% subset(term =="PC1_hirshorn_child")
adhd_att_symp_c_trio <- adhd_att_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model with adjustment for\nmaternal, paternal, and child ADHD PRS", Copsych = "ADHD-RS Inattention traits")

##Mother ADHD

adhd_att_symp_m <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_mor + SEX + adhd_prs_mor, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_m <- adhd_att_symp_m %>% subset(term =="PC1_hirshorn_mor")
adhd_att_symp_m <- adhd_att_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for child sex and either\nmaternal, paternal, or child ADHD PRS", Copsych = "ADHD-RS Inattention traits")

adhd_att_symp_m_trio <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_m_trio <- adhd_att_symp_m_trio %>% subset(term =="PC1_hirshorn_mor")
adhd_att_symp_m_trio <- adhd_att_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model with adjustment for\nmaternal, paternal, and child ADHD PRS", Copsych = "ADHD-RS Inattention traits")

##Father ADHD

adhd_att_symp_f <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_far + SEX + adhd_prs_far, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_f <- adhd_att_symp_f %>% subset(term =="PC1_hirshorn_far")
adhd_att_symp_f <- adhd_att_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for child sex and either\nmaternal, paternal, or child ADHD PRS", Copsych = "ADHD-RS Inattention traits")

adhd_att_symp_f_trio <- glm(adhd_rs_copsych_1_9 ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_f_trio <- adhd_att_symp_f_trio %>% subset(term =="PC1_hirshorn_far")
adhd_att_symp_f_trio <- adhd_att_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model with adjustment for\nmaternal, paternal, and child ADHD PRS", Copsych = "ADHD-RS Inattention traits")

##Impulsivity/Hyperactivity
##Child ADHD

adhd_ih_symp_c <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_child + SEX + adhd_prs, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_c <- adhd_ih_symp_c %>% subset(term =="PC1_hirshorn_child")
adhd_ih_symp_c <- adhd_ih_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for child sex and either\nmaternal, paternal, or child ADHD PRS", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

adhd_ih_symp_c_trio <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_child + PC1_hirshorn_mor + PC1_hirshorn_far + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_c_trio <- adhd_ih_symp_c_trio %>% subset(term =="PC1_hirshorn_child")
adhd_ih_symp_c_trio <- adhd_ih_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model with adjustment for\nmaternal, paternal, and child ADHD PRS", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

##Mother ADHD

adhd_ih_symp_m <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_mor + SEX + adhd_prs_mor, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_m <- adhd_ih_symp_m %>% subset(term =="PC1_hirshorn_mor")
adhd_ih_symp_m <- adhd_ih_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for child sex and either\nmaternal, paternal, or child ADHD PRS", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

adhd_ih_symp_m_trio <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_mor + PC1_hirshorn_child + PC1_hirshorn_far + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_m_trio <- adhd_ih_symp_m_trio %>% subset(term =="PC1_hirshorn_mor")
adhd_ih_symp_m_trio <- adhd_ih_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model with adjustment for\nmaternal, paternal, and child ADHD PRS", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

##Father ADHD

adhd_ih_symp_f <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_far + SEX + adhd_prs_far, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_f <- adhd_ih_symp_f %>% subset(term =="PC1_hirshorn_far")
adhd_ih_symp_f <- adhd_ih_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for child sex and either\nmaternal, paternal, or child ADHD PRS", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")

adhd_ih_symp_f_trio <- glm(adhd_rs_copsych_10_18 ~ PC1_hirshorn_far + PC1_hirshorn_child + PC1_hirshorn_mor + SEX + adhd_prs + adhd_prs_mor + adhd_prs_far, data = pref_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_f_trio <- adhd_ih_symp_f_trio %>% subset(term =="PC1_hirshorn_far")
adhd_ih_symp_f_trio <- adhd_ih_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model with adjustment for\nmaternal, paternal, and child ADHD PRS", Copsych = "ADHD-RS  Hyperactivity/Impulsivity traits")


##Forest plot ADHD RS traits 

forestadhdrssub <- rbind(adhd_total_symp_c, adhd_total_symp_c_trio, adhd_total_symp_m, adhd_total_symp_m_trio, adhd_total_symp_f, adhd_total_symp_f_trio, adhd_att_symp_c, adhd_att_symp_c_trio, adhd_att_symp_m, adhd_att_symp_m_trio, adhd_att_symp_f, adhd_att_symp_f_trio, adhd_ih_symp_c, adhd_ih_symp_c_trio, adhd_ih_symp_m, adhd_ih_symp_m_trio, adhd_ih_symp_f, adhd_ih_symp_f_trio)

forestadhdrssub$Condition <- factor(forestadhdrssub$Condition, levels=c("Adjusted for child sex and either\nmaternal, paternal, or child ADHD PRS", "Trio model with adjustment for\nmaternal, paternal, and child ADHD PRS"))
forestadhdrssub$Group <- factor(forestadhdrssub$Group, levels=c("Child direct effect", "Child PGS", "Paternal indirect effect", "Paternal PGS", "Maternal indirect effect", "Maternal PGS"))
forestadhdrssub$Copsych <- factor(forestadhdrssub$Copsych, levels=c("ADHD-RS Total traits","ADHD-RS Inattention traits","ADHD-RS  Hyperactivity/Impulsivity traits"))
forestadhdrssub <- forestadhdrssub %>% mutate(significant = ifelse(p.value<0.05,1,0))

theme_set(theme_bw())
redPalette <- c("#E69F00", "#009E73")

p_adhd_rs_sub_trio_adhd_prs = ggplot(data=forestadhdrssub, aes(x = Group,y = estimate, ymin = conf.low, ymax = conf.high, na.rm = TRUE))+
  geom_point(size = 2, shape = 19, aes(col=Condition)) +
  geom_hline(aes(fill=Group),yintercept =0, linetype=2)+
  xlab("") + ylab("Beta Estimate (95% Confidence Interval)") +
  coord_flip() + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high,col=Condition),width=0.2,cex=0.5)+ 
  scale_y_continuous(limits = c(-2,1.2))+  # Set x-axis limits
  facet_wrap(~Copsych,strip.position="top",nrow=1,scales = "free_x", labeller = label_wrap_gen(width=15) ) +
  scale_colour_manual(values=redPalette) +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_text(size = 14, face="bold"),
        axis.text.x=element_text(size = 16, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 16, face="bold", margin = margin(0.6,0,0.6,0, "cm")), strip.background = element_rect(fill = "white"),legend.title = element_blank(), legend.text = element_text(size = 14)) +
  geom_text(aes(label = ifelse(significant == 1, "*", ""), col=Condition), vjust = -0.1, size = 6, show.legend=FALSE) 


##Final plot

library(patchwork)
p_fin <- (p_adhd_trio_adhd_prs / p_adhd_rs_sub_trio_adhd_prs)




#################################### Table 1 ######################################

##Including all individuals in copsych evaluation 

summary(pref_data$PC1_hirshorn_child)
summary(pref_data$PC1_hirshorn_far)
summary(pref_data$PC1_hirshorn_mor)

count <- pref_data %>% filter(PC1_hirshorn_child!="NA" & PC1_hirshorn_far!="NA" & PC1_hirshorn_mor!="NA" & adhd!="NA")
count <- pref_data %>% filter(PC1_hirshorn_child!="NA" & PC1_hirshorn_far!="NA" & PC1_hirshorn_mor!="NA" & adhd_rs_copsych_1_18!="NA")

tableone_diet <- pref_data %>% mutate(pufa = relevel(factor(pufa), ref="Placebo"))
tableone_diet <- tableone_diet %>% mutate(copsych_participation_yn = 1)
tableone_diet <- tableone_diet %>% mutate(Diet_PC2 = as.numeric(Diet_PC2))


tableone_diet <- tableone_diet %>%
  mutate(copsych_participation_yn = case_when(
    copsych_participation == "copsych visit" ~ "Yes",
    copsych_participation == "copsych questionnaire" ~ "Yes",
    copsych_participation == "no copsych data" ~ "No"
  ))

tableone_diet <- tableone_diet %>% mutate(pc1_bi= ifelse(PC1_hirshorn_mor>0,1,0))

summary(tableone_diet$copsych_participation_yn)
table(tableone_diet$copsych_participation_yn)

tableone_diet <- tableone_diet  %>% mutate(pufa = relevel(factor(pufa), ref="Placebo"))
tableone_diet <- tableone_diet %>% mutate(copsych_parti = ifelse(is.na(asd),0,1))
summary(tableone_diet$copsych_parti)
table(tableone_diet$copsych_parti)

### Table 1 - stratified according to outcome
diet_table <- CreateTableOne(vars = c("PC1_hirshorn_child","PC1_hirshorn_mor","PC1_hirshorn_far", "SEX", "AGE_DAYS", "bw", "mat_edu", "pat_edu","income","alc_preg_YN_2010","smoking_preg","parity","Mother_BMI","Diet_PC2","Motherage","Fatherage","dvit", "pufa") , strata = "adhd" , data = subset(tableone_diet, (PC1_hirshorn_child!="NA"& PC1_hirshorn_far!="NA" & PC1_hirshorn_mor!="NA" & adhd!="NA")), factorVars = c("pufa","parity","alc_preg_YN_2010","smoking_preg","SEX","mat_edu", "pat_edu", "income", "dvit"))
print(diet_table, quote = TRUE, noSpaces = TRUE, addOverall = TRUE)


#################################### Drop out analysis ######################################

##Drop out analysis comparing full cohort to complete trios

summary(pref_data$PC1_hirshorn_child)
summary(pref_data$PC1_hirshorn_far)
summary(pref_data$PC1_hirshorn_mor)

#Number of complete trios (and information on ADHD diagnosis)
count <- pref_data %>% filter(PC1_hirshorn_child!="NA" & PC1_hirshorn_far!="NA" & PC1_hirshorn_mor!="NA" & copsych_participation!="no copsych data"& adhd!="NA")
#Number of individuals with at least one PRS (and information on ADHD diagnosis)
count <- pref_data %>% filter((PC1_hirshorn_child!="NA"| PC1_hirshorn_far!="NA" | PC1_hirshorn_mor!="NA") & copsych_participation!="no copsych data")

drop_out_df <- pref_data %>% 
  mutate(copsych_participation = as.character(copsych_participation)) %>%
  mutate(drop_out = ifelse(!is.na(PC1_hirshorn_child) & !is.na(PC1_hirshorn_far) & 
                             !is.na(PC1_hirshorn_mor) & 
                             !is.na(adhd), 1, 0))
drop_out_df <- drop_out_df %>% mutate(Diet_PC2 = as.numeric(Diet_PC2))
table(drop_out_df$drop_out)

### Drop out analysis
drop_out_table <- CreateTableOne(vars = c("PC1_hirshorn_child","PC1_hirshorn_mor","PC1_hirshorn_far", "SEX", "AGE_DAYS", "bw", "mat_edu", "pat_edu","income","alc_preg_YN_2010","smoking_preg","parity","Mother_BMI","Diet_PC2", "Motherage","Fatherage","dvit", "pufa") , strata = "drop_out" , data = drop_out_df, factorVars = c("pufa","parity","alc_preg_YN_2010","smoking_preg","SEX","mat_edu", "pat_edu", "income", "dvit"), addOverall = TRUE)
print(drop_out_table, quote = TRUE, noSpaces = TRUE, addoverall = TRUE)


table(pref_data$copsych_participation)

