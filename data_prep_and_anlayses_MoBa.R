

################## MoBa replication - diet ADHD trio analyses ##################

npr_full <- preload_npr(
  npr_data_root_dir = "//ess01/P471/data/durable/data/NPR/processed/",
  npr_filename = "npr2024.sav")

available_variables(source="npr")

adhd_nprcodes <- c("F900","F988", "F901", "F908", "F909")
adhd_att_nprcodes <- c("F988")
adhd_comb_nprcodes <- c("F900", "F901", "F908", "F909")

npr_groups <- c(paste0("adhddiagnoses = ",paste0(adhd_nprcodes,collapse=",")),
                paste0("adhd_comb_diagnoses = ",paste0(adhd_comb_nprcodes,collapse=",")),
                paste0("adhd_att_diagnosess = ",paste0(adhd_att_nprcodes,collapse=",")))

pheno_data <- curate_dataset(variables_required=list(moba = c("KJONN","rsdbd_adhd_c_8yr","rsdbd_adhd_c_14m","rsdbd_ina_c_8yr","rsdbd_ina_c_14m","rsdbd_hyp_c_8yr","rsdbd_hyp_c_14m"),
                                                 npr = npr_groups),
                         return_items = F, out_format="merged_df",
                         exclusions=NULL, 
                         recursive=TRUE,
                         group_all=TRUE,  
                         dx_groupname=NULL,
                         dx_owner="child",
                         npr_full=npr_full) 

#Generating child ID based on pregnancy ID and twin indicator
pheno_data <- pheno_data |> 
   mutate(child_id = paste0(preg_id, BARN_NR))

####### Loading diet data ######

load(file="N:/durable/projects/diet_adhd_trio/data/101_foodgroups.Rdata")

####### Reading in PGS ######

#Check that the pgs is available
available_pgs() %>% 
  pgs_search("diet")

#Retrieve the pgs (and an educational attainment PGS for validation)
available_pgs()

geno_dat <- fetch_pgs("diet_pc1")

geno_dat_validation <- fetch_pgs("ea3",
                                 pgs_directory = "//ess01/P471/data/durable/common/pgs_directory/pgs/",
                                 maf="0.01",clump="250_1_0.1",pgs_software="prsice2")

geno_dat_prcd <- process_pgs(geno_dat, return_covs = F)
geno_dat_val_prcd <- process_pgs(geno_dat_validation, return_covs = F)
geno_dat_val_prcd <- pgs_pca(geno_dat_val_prcd, indid= "IID", pgs_var_stem = "ea3")

# Join 
geno_dat_kids <- geno_dat_prcd |> 
  dplyr::select(-IID, -FID, -`diet-pc1_pgs`, -m_id, -f_id) |> 
  filter(!is.na(preg_id)) |> 
  left_join(geno_dat_val_prcd |> 
              dplyr::select(preg_id,BARN_NR,ea3.pgs.pc)) %>% 
  dplyr::rename("diet_pgs_barn" = `diet-pc1_pgs_res`, "ea_pgs_barn" = ea3.pgs.pc)

geno_dat_kids$BARN_NR <- as.numeric(geno_dat_kids$BARN_NR)

geno_dat_mums <- geno_dat_prcd |> 
  dplyr::select(-IID, -FID, -`diet-pc1_pgs`, -preg_id, -BARN_NR, -f_id) |> 
  filter(!is.na(m_id)) |> 
  left_join(geno_dat_val_prcd |> 
              dplyr::select(m_id,ea3.pgs.pc))|> 
  dplyr::rename("diet_pgs_mor" = `diet-pc1_pgs_res`, "ea_pgs_mor" = ea3.pgs.pc)

geno_dat_dads <- geno_dat_prcd |> 
  dplyr::select(-IID, -FID, -`diet-pc1_pgs`, -m_id, -preg_id, -BARN_NR) |> 
  filter(!is.na(f_id)) |> 
  left_join(geno_dat_val_prcd |> 
              dplyr::select(f_id,ea3.pgs.pc))|> 
  dplyr::rename("diet_pgs_far" = `diet-pc1_pgs_res`, "ea_pgs_far" = ea3.pgs.pc)
  
#Join everything to phenotypic data
merged_data <- pheno_data |> 
  left_join(geno_dat_kids) |> 
  left_join(geno_dat_mums, by = "m_id") |> 
  left_join(geno_dat_dads, by = "f_id") |> 
  left_join(moba_100_fg, by = c("preg_id"="PREG.ID.2306"))

#Perform exclusions: multiple pregnancies (Keep on child from each twin pair)
merged_data <- merged_data %>% dplyr::filter(BARN_NR==1)

#Generating sex variable
table(merged_data$KJONN_raw)
str(merged_data$KJONN_raw)
merged_data <- merged_data %>% mutate(sex = ifelse(KJONN_raw==1, "male", ifelse(KJONN_raw==2, "female", NA)))
table(merged_data$sex)

#Recode ADHD NPR diagnosis
#ADHD
table(merged_data$received_dx_adhddiagnoses_npr)
merged_data <- merged_data %>% mutate(adhd = ifelse(received_dx_adhddiagnoses_npr=="yes", 1, ifelse(received_dx_adhddiagnoses_npr=="no",0, NA)))
table(merged_data$adhd)
merged_data$adhd <- as.factor(merged_data$adhd)
#Combined type
table(merged_data$received_dx_adhd_comb_diagnoses_npr)
merged_data <- merged_data %>% mutate(adhd_comb = ifelse(received_dx_adhd_comb_diagnoses_npr=="yes", 1, ifelse(received_dx_adhd_comb_diagnoses_npr=="no",0, NA)))
table(merged_data$adhd_comb)
merged_data$adhd_comb <- as.factor(merged_data$adhd_comb)
#Inattentive type
table(merged_data$received_dx_adhd_att_diagnosess_npr)
merged_data <- merged_data %>% mutate(adhd_att = ifelse(received_dx_adhd_att_diagnosess_npr=="yes", 1, ifelse(received_dx_adhd_att_diagnosess_npr=="no",0, NA)))
table(merged_data$adhd_att)
merged_data$adhd_att <- as.factor(merged_data$adhd_att)

#Generating extra variable for ADHD diagnoses groups recieved minimum twice
#ADHD
table(merged_data$received_dx_2x_adhddiagnoses_npr)
merged_data <- merged_data %>% mutate(adhd2 = ifelse(received_dx_2x_adhddiagnoses_npr=="yes", 1, ifelse(received_dx_2x_adhddiagnoses_npr=="no",0, NA)))
table(merged_data$adhd2)
merged_data$adhd2 <- as.factor(merged_data$adhd2)
#Combined type
table(merged_data$received_dx_2x_adhd_comb_diagnoses_npr)
merged_data <- merged_data %>% mutate(adhd_comb2 = ifelse(received_dx_2x_adhd_comb_diagnoses_npr=="yes", 1, ifelse(received_dx_2x_adhd_comb_diagnoses_npr=="no",0, NA)))
table(merged_data$adhd_comb2)
merged_data$adhd_comb2 <- as.factor(merged_data$adhd_comb2)
#Inattentive type
table(merged_data$received_dx_2x_adhd_att_diagnosess_npr)
merged_data <- merged_data %>% mutate(adhd_att2 = ifelse(received_dx_2x_adhd_att_diagnosess_npr=="yes", 1, ifelse(received_dx_2x_adhd_att_diagnosess_npr=="no",0, NA)))
table(merged_data$adhd_att2)
merged_data$adhd_att2 <- as.factor(merged_data$adhd_att2)

##Testing for assortative mating
cor.test(merged_data$diet_pgs_mor, merged_data$diet_pgs_far)


################### Plotting against maternal diet ####################

##Approximately 100 food groups 
load(file="N:/durable/projects/diet_adhd_trio/data/101_foodgroups.Rdata")
str(moba_100_fg)
moba_diet <- merged_data %>% dplyr::select(White_bread:Sweet_biscuits, diet_pgs_mor, KCAL)
moba_diet %>% sjPlot::view_df()

#Cleaning names 
moba_diet <- data.frame(lapply(moba_diet, function(x) {
  attr(x, "label") <- NULL
  return(x)
}))

names(moba_diet) <- gsub("_", " ", names(moba_diet))

#Histogram of all diet categories
moba_diet <- moba_diet %>%
  mutate(across(everything(), as.numeric))

# Convert the data to long format
moba_diet_long <- moba_diet %>% select(`White bread`:`Sweet biscuits`) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

# Create histograms for all variables in one plot 
ggplot(moba_diet_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.5) +
  facet_wrap(~ variable, scales = "free") +  # Separate histogram per variable
  labs(
    title = "Histograms for Multiple Variables",
    x = "Value",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14)

#Scaling values
scale2 <- function(x, na.rm = FALSE) (x %>% scale())
moba_diet <- moba_diet %>% dplyr::mutate(across(c(`White bread`:`Sweet biscuits`), scale2)) 

glimpse(moba_diet)

allfoodgroups_morprs <- moba_diet

getmixedSTAT2 <- function(x){
  
  mdl <- glm(data = x, foodgroup_result ~ food_group_value, family= "gaussian")
  res <- mdl %>% broom::tidy(., exp=F, conf.int=T) %>% filter(term != "(Intercept)")
}

summary(allfoodgroups_morprs)
names(allfoodgroups_morprs)

output3<- allfoodgroups_morprs   %>% 
  gather(foodgroup_factor,foodgroup_result, `White bread`:`Sweet biscuits`) %>% 
  gather(food_group,food_group_value,`diet pgs mor`) %>% 
  filter(!is.na(food_group_value)) %>% 
  filter(!is.infinite(food_group_value)) %>% 
  group_by(foodgroup_factor, food_group) %>% do(getmixedSTAT2(x = .)) %>% ungroup()

output3$p.valuefdr <- p.adjust(output3$p.value,method="fdr")
output3<- output3 %>% mutate(FILL=ifelse(estimate >= 0.041, "High", ifelse(estimate < 0.041 & estimate > -0.03052, "mid", ifelse(estimate <-0.03052, "low", NA)))) %>% mutate(FILL=as.factor(FILL))
output3<- output3 %>% mutate(p.value=ifelse(p.value<0.001,"<0.001", round(p.value,3)))
output3<- output3 %>% mutate(p.valuedfrbi=ifelse(p.valuefdr<0.05,1,0))
output3<-output3[order(output3$estimate),]  

order <- output3$foodgroup_factor

output3$foodgroup_factor<- factor(output3$foodgroup_factor, levels= output3$foodgroup_factor)

saved<-output3$foodgroup_factor


a <- ggplot(output3 %>% filter(FILL != "mid"), aes(y = foodgroup_factor, 
                                                   x = as.numeric(as.character((estimate))), 
                                                   xmin = as.numeric(as.character((conf.low))), 
                                                   xmax = as.numeric(as.character((conf.high))), 
                                                   fill = FILL)) + 
  geom_bar(stat = "identity", position = "dodge", alpha = 1) + 
  geom_errorbarh(height = 0.3) +
  geom_hline(yintercept = seq(-0.3, 0.45, by = 0.1), color = "gray80", linetype = "dotted") +
  xlim(c(-0.12, 0.12)) + 
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
  label = paste0("p = ", formatC(p.value, format = "e", digits = 2), # Format p-value in scientific notation
                 ifelse(p.valuefdr < 0.05, " *", "")), # Add stars for significance
  x = 0.1), # Adjust positioning
  hjust = 0.0, # Adjust horizontal alignment
  size = 3.5) # Adjust text size

a

## Performing PCA analysis to test the ability of maternal PGS to predict PC1 ##

#Merge tot_energy to diet data
pca_fg_moba <- na.omit(allfoodgroups_morprs)

pca_fg_moba$diet_pgs_mor <- pca_fg_moba$`diet pgs mor`

#Calibrate for total energy before PCA
vars_to_adjust <- names(pca_fg_moba)[1:98]

for (v in vars_to_adjust) {
  pca_fg_moba[[paste0(v, "_adj")]] <- residuals(lm(pca_fg_moba[[v]] ~ pca_fg_moba$`KCAL`)) + mean(pca_fg_moba[[v]], na.rm = TRUE)
}

#Remove non adjusted values
pca_fg_moba <- pca_fg_moba[c("diet_pgs_mor", grep("_adj$", names(pca_fg_moba), value = TRUE))]

# Step 1: Separate the diet Pgs column
dietpgsmor <- pca_fg_moba$diet_pgs_mor

# Step 2: Remove the diet_pgs_mor column for PCA analysis
pca_data_fg_moba <- pca_fg_moba[, !names(pca_fg_moba) %in% 'diet_pgs_mor']

# Step 3: Perform the PCA
data.pca_fg_moba <- princomp(pca_data_fg_moba)

# Step 4: Summarize the PCA
summary(data.pca_fg_moba)

# Step 5: Keep track of the PCA scores with the corresponding ABCNO
pca_scores_with_IDs_moba <- data.frame(diet_pgs_mor = dietpgsmor, data.pca_fg_moba$scores)

# View the PCA scores with the corresponding IDs
head(pca_scores_with_IDs_moba)

#install.packages("factoextra")
library(factoextra)

data.pca_fg_moba$loadings[, 1:10]
fviz_eig(data.pca_fg_moba, addlabels = TRUE)
fviz_cos2(data.pca_fg_moba, choice = "var", axes = 1)
fviz_cos2(data.pca_fg_moba, choice = "var", axes = 2)

fviz_pca_var(data.pca_fg_moba, col.var = "cos2", axes = c(1,2),
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

cor.test(pca_scores_with_IDs_moba$Comp.1, pca_scores_with_IDs_moba$diet_pgs_mor)
cor.test(pca_scores_with_IDs_moba$Comp.1, pca_scores_with_IDs_moba$diet_pgs_mor)$estimate^2 
cor.test(pca_scores_with_IDs_moba$Comp.2, pca_scores_with_IDs_moba$diet_pgs_mor)
cor.test(pca_scores_with_IDs_moba$Comp.2, pca_scores_with_IDs_moba$diet_pgs_mor)$estimate^2 

fviz_pca_var(data.pca_fg_moba, col.var = "cos2", axes = c(2,3),
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

fviz_pca_var(data.pca_fg_moba, col.var = "cos2", axes = c(3,4),
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

cor.test(pca_scores_with_IDs_moba$Comp.3, pca_scores_with_IDs_moba$diet_pgs_mor)
cor.test(pca_scores_with_IDs_moba$Comp.4, pca_scores_with_IDs_moba$diet_pgs_mor)






##### Direct associations #####

#ADHD diagnosis
model <- glm(adhd ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
broom::glance(model)

model <- glm(adhd ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
broom::glance(model)

model <- glm(adhd ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
broom::glance(model)

glm(adhd ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp=T)

glm(adhd ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp =T)

glm(adhd ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp =T)


glm(adhd_comb ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp=T)

glm(adhd_comb ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp =T)

glm(adhd_comb ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp =T)


glm(adhd_att ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp=T)

glm(adhd_att ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp =T)

glm(adhd_att ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp =T)


#ADHD traits
model <- lm(rsdbd_adhd_c_8yr ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data)
broom::glance(model)

model <- lm(rsdbd_adhd_c_8yr ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data)
broom::glance(model)

model <- lm(rsdbd_adhd_c_8yr ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), data = merged_data)
broom::glance(model)

lm(rsdbd_adhd_c_8yr ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% broom::tidy(conf.int=T)

lm(rsdbd_adhd_c_8yr ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% broom::tidy(conf.int=T)

lm(rsdbd_adhd_c_8yr ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% broom::tidy(conf.int=T)

#Inattention subscale
lm(rsdbd_ina_c_8yr ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% broom::tidy(conf.int=T)

lm(rsdbd_ina_c_8yr ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% broom::tidy(conf.int=T)

lm(rsdbd_ina_c_8yr ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% broom::tidy(conf.int=T)

#Hyperactivity subscale
lm(rsdbd_hyp_c_8yr ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% broom::tidy(conf.int=T)

lm(rsdbd_hyp_c_8yr ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% broom::tidy(conf.int=T)

lm(rsdbd_hyp_c_8yr ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% broom::tidy(conf.int=T)





##### Trio analyses #####


#-----------------------------------------------------------------------------------------------------------------------------------------------------------

model <- glm(adhd ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
glance(model)
summary(model)

model <- glm(adhd ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
glance(model)

model <- glm(adhd ~ diet_pgs_barn + diet_pgs_far + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
glance(model)

glm(adhd ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp=T)

glm(adhd ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp =T)

glm(adhd ~ diet_pgs_barn + diet_pgs_far + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp =T)


glm(adhd_comb ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp=T)

glm(adhd_comb ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp =T)

glm(adhd_comb ~ diet_pgs_barn + diet_pgs_far + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp =T)


glm(adhd_att ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp=T)

glm(adhd_att ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp =T)

glm(adhd_att ~ diet_pgs_barn + diet_pgs_far + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp =T)


#ADHD traits
model <- lm(rsdbd_adhd_c_8yr ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
broom::glance(model)

model <- lm(rsdbd_ina_c_8yr ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
broom::glance(model)

model <- lm(rsdbd_hyp_c_8yr ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
broom::glance(model)

lm(rsdbd_adhd_c_8yr ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int=T)

lm(rsdbd_adhd_c_8yr ~ diet_pgs_far + diet_pgs_mor + diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int=T)

lm(rsdbd_adhd_c_8yr ~ diet_pgs_barn + diet_pgs_mor + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int=T)

#Inattention subscale
lm(rsdbd_ina_c_8yr ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int=T)

lm(rsdbd_ina_c_8yr ~ diet_pgs_far + diet_pgs_mor + diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int=T)

lm(rsdbd_ina_c_8yr ~ diet_pgs_barn + diet_pgs_mor + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int=T)

#Hyperactivity subscale
lm(rsdbd_hyp_c_8yr ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int=T)

lm(rsdbd_hyp_c_8yr ~ diet_pgs_far + diet_pgs_mor + diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int=T)

lm(rsdbd_hyp_c_8yr ~ diet_pgs_barn + diet_pgs_mor + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int=T)





################################## Figures - forest plot trio model ###########################
 
##Diagnoses

##Child ADHD

adhd_c <- glm(adhd ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp =T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_c <- adhd_c %>% subset(term =="diet_pgs_barn")
adhd_c <- adhd_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD")

adhd_c_trio <- glm(adhd ~ diet_pgs_barn + diet_pgs_far + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp =T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_c_trio <- adhd_c_trio %>% subset(term =="diet_pgs_barn")
adhd_c_trio <- adhd_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD")

##Mother ADHD

adhd_m <- glm(adhd ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_m <- adhd_m %>% subset(term =="diet_pgs_mor")
adhd_m <- adhd_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD")

adhd_m_trio <- glm(adhd ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_m_trio <- adhd_m_trio %>% subset(term =="diet_pgs_mor")
adhd_m_trio <- adhd_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD")

##Father ADHD

adhd_f <- glm(adhd ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp =T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_f <- adhd_f %>% subset(term =="diet_pgs_far")
adhd_f <- adhd_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD")

adhd_f_trio <- glm(adhd ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp =T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_f_trio <- adhd_f_trio %>% subset(term =="diet_pgs_far")
adhd_f_trio <- adhd_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD")

##Child ADD

add_c <- glm(adhd_att ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp =T) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_c <- add_c %>% subset(term =="diet_pgs_barn")
add_c <- add_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_c_trio <- glm(adhd_att ~ diet_pgs_barn + diet_pgs_far + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp =T) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_c_trio <- add_c_trio %>% subset(term =="diet_pgs_barn")
add_c_trio <- add_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Mother ADD

add_m <- glm(adhd_att ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_m <- add_m %>% subset(term =="diet_pgs_mor")
add_m <- add_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_m_trio <- glm(adhd_att ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_m_trio <- add_m_trio %>% subset(term =="diet_pgs_mor")
add_m_trio <- add_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Father ADD

add_f <- glm(adhd_att ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp =T) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_f <- add_f %>% subset(term =="diet_pgs_far")
add_f <- add_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_f_trio <- glm(adhd_att ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp =T) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_f_trio <- add_f_trio %>% subset(term =="diet_pgs_far")
add_f_trio <- add_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Child ADHD - Combined Presentation

adhd1_c <- glm(adhd_comb ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp =T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_c <- adhd1_c %>% subset(term =="diet_pgs_barn")
adhd1_c <- adhd1_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_c_trio <- glm(adhd_comb ~ diet_pgs_barn + diet_pgs_far + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp =T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_c_trio <- adhd1_c_trio %>% subset(term =="diet_pgs_barn")
adhd1_c_trio <- adhd1_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

##Mother ADHD - Combined Presentation

adhd1_m <- glm(adhd_comb ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_m <- adhd1_m %>% subset(term =="diet_pgs_mor")
adhd1_m <- adhd1_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_m_trio <- glm(adhd_comb ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_m_trio <- adhd1_m_trio %>% subset(term =="diet_pgs_mor")
adhd1_m_trio <- adhd1_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

##Father ADHD - Combined Presentation

adhd1_f <- glm(adhd_comb ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% broom::tidy(conf.int=T, exp =T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_f <- adhd1_f %>% subset(term =="diet_pgs_far")
adhd1_f <- adhd1_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_f_trio <- glm(adhd_comb ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) %>% tidy(conf.int=T, exp =T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_f_trio <- adhd1_f_trio %>% subset(term =="diet_pgs_far")
adhd1_f_trio <- adhd1_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

forestadhd <- rbind(adhd_c, adhd_c_trio, adhd_m, adhd_m_trio, adhd_f, adhd_f_trio, add_c, add_c_trio, add_m, add_m_trio, add_f, add_f_trio, adhd1_c, adhd1_c_trio, adhd1_m, adhd1_m_trio, adhd1_f, adhd1_f_trio)
forestadhd$Condition <- factor(forestadhd$Condition, levels=c("norm", "trio"))
forestadhd$Group <- factor(forestadhd$Group, levels=c("Child direct effect", "Child PGS", "Paternal indirect effect", "Paternal PGS", "Maternal indirect effect", "Maternal PGS"))
forestadhd$Copsych <- factor(forestadhd$Copsych, levels=c("ADHD", "ADHD  Inattentive Presentation", "ADHD  Combined Presentation"))
forestadhd <- forestadhd %>% mutate(significant = ifelse(p.value<0.05,1,0))

theme_set(theme_bw())
redPalette <- c("#E69F00", "#009E73")

p_adhd_trio = ggplot(data=forestadhd, aes(x = Group,y = estimate, ymin = conf.low, ymax = conf.high, na.rm = TRUE))+
  geom_point(size = 2.0, shape = 19, aes(col=Condition)) +
  geom_hline(aes(fill=Group),yintercept =1, linetype=2)+
  xlab("") + ylab("Odds Ratio (95% Confidence Interval)") +
  coord_flip() + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high,col=Condition),width=0.2,cex=0.5)+ 
  scale_y_continuous(limits = c(0.75,1.3))+ 
  facet_wrap(~Copsych,strip.position="top",nrow=1,scales = "free_x", labeller = label_wrap_gen(width=8) ) +
  scale_colour_manual(values=redPalette) +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_text(size = 14, face="bold"),
        axis.text.x=element_text(size = 16, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 16, face="bold", margin = margin(0.6,0,0.6,0, "cm")), legend.position="None", strip.background = element_rect(fill = "white")) +
  geom_text(aes(label = ifelse(significant == 1, "*", ""), col=Condition), vjust = -0.1, size = 6) 


## ADHD RS-DBD Total traits, Inattention and IH traits

##Total traits
##Child ADHD

adhd_total_symp_c <- glm(rsdbd_adhd_c_8yr ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_c <- adhd_total_symp_c %>% subset(term =="diet_pgs_barn")
adhd_total_symp_c <- adhd_total_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Total traits")

adhd_total_symp_c_trio <- glm(rsdbd_adhd_c_8yr ~ diet_pgs_barn + diet_pgs_mor + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_c_trio <- adhd_total_symp_c_trio %>% subset(term =="diet_pgs_barn")
adhd_total_symp_c_trio <- adhd_total_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Total traits")

##Mother ADHD

adhd_total_symp_m <- glm(rsdbd_adhd_c_8yr ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_m <- adhd_total_symp_m %>% subset(term =="diet_pgs_mor")
adhd_total_symp_m <- adhd_total_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Total traits")

adhd_total_symp_m_trio <- glm(rsdbd_adhd_c_8yr ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_m_trio <- adhd_total_symp_m_trio %>% subset(term =="diet_pgs_mor")
adhd_total_symp_m_trio <- adhd_total_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Total traits")

##Father ADHD

adhd_total_symp_f <- glm(rsdbd_adhd_c_8yr ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_f <- adhd_total_symp_f %>% subset(term =="diet_pgs_far")
adhd_total_symp_f <- adhd_total_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Total traits")

adhd_total_symp_f_trio <- glm(rsdbd_adhd_c_8yr ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_f_trio <- adhd_total_symp_f_trio %>% subset(term =="diet_pgs_far")
adhd_total_symp_f_trio <- adhd_total_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Total traits")

##Inattention
##Child ADHD

adhd_att_symp_c <- glm(rsdbd_ina_c_8yr ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_c <- adhd_att_symp_c %>% subset(term =="diet_pgs_barn")
adhd_att_symp_c <- adhd_att_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Inattention traits")

adhd_att_symp_c_trio <- glm(rsdbd_ina_c_8yr ~ diet_pgs_barn + diet_pgs_mor + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_c_trio <- adhd_att_symp_c_trio %>% subset(term =="diet_pgs_barn")
adhd_att_symp_c_trio <- adhd_att_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Inattention traits")

##Mother ADHD

adhd_att_symp_m <- glm(rsdbd_ina_c_8yr ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_m <- adhd_att_symp_m %>% subset(term =="diet_pgs_mor")
adhd_att_symp_m <- adhd_att_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Inattention traits")

adhd_att_symp_m_trio <- glm(rsdbd_ina_c_8yr ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_m_trio <- adhd_att_symp_m_trio %>% subset(term =="diet_pgs_mor")
adhd_att_symp_m_trio <- adhd_att_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Inattention traits")

##Father ADHD

adhd_att_symp_f <- glm(rsdbd_ina_c_8yr ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_f <- adhd_att_symp_f %>% subset(term =="diet_pgs_far")
adhd_att_symp_f <- adhd_att_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Inattention traits")

adhd_att_symp_f_trio <- glm(rsdbd_ina_c_8yr ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_f_trio <- adhd_att_symp_f_trio %>% subset(term =="diet_pgs_far")
adhd_att_symp_f_trio <- adhd_att_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Inattention traits")

##Hyperactivity/Impulsivity
##Child ADHD

adhd_ih_symp_c <- glm(rsdbd_hyp_c_8yr ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_c <- adhd_ih_symp_c %>% subset(term =="diet_pgs_barn")
adhd_ih_symp_c <- adhd_ih_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Hyperactivity/Impulsivity traits")

adhd_ih_symp_c_trio <- glm(rsdbd_hyp_c_8yr ~ diet_pgs_barn + diet_pgs_mor + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_c_trio <- adhd_ih_symp_c_trio %>% subset(term =="diet_pgs_barn")
adhd_ih_symp_c_trio <- adhd_ih_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Hyperactivity/Impulsivity traits")

##Mother ADHD

adhd_ih_symp_m <- glm(rsdbd_hyp_c_8yr ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_m <- adhd_ih_symp_m %>% subset(term =="diet_pgs_mor")
adhd_ih_symp_m <- adhd_ih_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Hyperactivity/Impulsivity traits")

adhd_ih_symp_m_trio <- glm(rsdbd_hyp_c_8yr ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_m_trio <- adhd_ih_symp_m_trio %>% subset(term =="diet_pgs_mor")
adhd_ih_symp_m_trio <- adhd_ih_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Hyperactivity/Impulsivity traits")

##Father ADHD

adhd_ih_symp_f <- glm(rsdbd_hyp_c_8yr ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_f <- adhd_ih_symp_f %>% subset(term =="diet_pgs_far")
adhd_ih_symp_f <- adhd_ih_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Hyperactivity/Impulsivity traits")

adhd_ih_symp_f_trio <- glm(rsdbd_hyp_c_8yr ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) %>% tidy(conf.int = TRUE) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_f_trio <- adhd_ih_symp_f_trio %>% subset(term =="diet_pgs_far")
adhd_ih_symp_f_trio <- adhd_ih_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Hyperactivity/Impulsivity traits")


##Forest plot ADHD RS-DBD traits 

forestadhdrssub <- rbind(adhd_total_symp_c, adhd_total_symp_c_trio, adhd_total_symp_m, adhd_total_symp_m_trio, adhd_total_symp_f, adhd_total_symp_f_trio, adhd_att_symp_c, adhd_att_symp_c_trio, adhd_att_symp_m, adhd_att_symp_m_trio, adhd_att_symp_f, adhd_att_symp_f_trio, adhd_ih_symp_c, adhd_ih_symp_c_trio, adhd_ih_symp_m, adhd_ih_symp_m_trio, adhd_ih_symp_f, adhd_ih_symp_f_trio)

forestadhdrssub$Condition <- factor(forestadhdrssub$Condition, levels=c("Adjusted for sex and birth year", "Trio model"))
forestadhdrssub$Group <- factor(forestadhdrssub$Group, levels=c("Child direct effect", "Child PGS", "Paternal indirect effect", "Paternal PGS", "Maternal indirect effect", "Maternal PGS"))
forestadhdrssub$Copsych <- factor(forestadhdrssub$Copsych, levels=c("ADHD RS-DBD Total traits","ADHD RS-DBD Inattention traits","ADHD RS-DBD Hyperactivity/Impulsivity traits"))
forestadhdrssub <- forestadhdrssub %>% mutate(significant = ifelse(p.value<0.05,1,0))

theme_set(theme_bw())
redPalette <- c("#E69F00", "#009E73")

p_adhd_rs_sub_trio = ggplot(data=forestadhdrssub, aes(x = Group,y = estimate, ymin = conf.low, ymax = conf.high, na.rm = TRUE))+
  geom_point(size = 2.0, shape = 19, aes(col=Condition)) +
  geom_hline(aes(fill=Group),yintercept =0, linetype=2)+
  xlab("") + ylab("Beta Estimate (95% Confidence Interval)") +
  coord_flip() + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high,col=Condition),width=0.2,cex=0.5)+ 
  scale_y_continuous(limits = c(-0.5,0.5))+  # Set x-axis limits
  facet_wrap(~Copsych,strip.position="top",nrow=1,scales = "free_x", labeller = label_wrap_gen (width=15) ) +
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
























########### Performing all analyses with clustered robust SE ##############



################################## Figures - forest plot trio model - clustered standard error ###########################


##Diagnoses

##Child ADHD

adhd_c <- glm(adhd ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
clustered_se <- vcovCL(adhd_c, cluster = ~ m_id)
summary_cl <- coeftest(adhd_c, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd_c <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_c <- adhd_c %>% subset(term =="diet_pgs_barn")
adhd_c <- adhd_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD")

adhd_c_trio <- glm(adhd ~ diet_pgs_barn + diet_pgs_far + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd_c_trio, cluster = ~ m_id)
summary_cl <- coeftest(adhd_c_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd_c_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_c_trio <- adhd_c_trio %>% subset(term =="diet_pgs_barn")
adhd_c_trio <- adhd_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD")

##Mother ADHD

adhd_m <- glm(adhd ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd_m, cluster = ~ m_id)
summary_cl <- coeftest(adhd_m, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd_m <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_m <- adhd_m %>% subset(term =="diet_pgs_mor")
adhd_m <- adhd_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD")

adhd_m_trio <- glm(adhd ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd_m_trio, cluster = ~ m_id)
summary_cl <- coeftest(adhd_m_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd_m_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_m_trio <- adhd_m_trio %>% subset(term =="diet_pgs_mor")
adhd_m_trio <- adhd_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD")

##Father ADHD

adhd_f <- glm(adhd ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
clustered_se <- vcovCL(adhd_f, cluster = ~ m_id)
summary_cl <- coeftest(adhd_f, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd_f <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_f <- adhd_f %>% subset(term =="diet_pgs_far")
adhd_f <- adhd_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD")

adhd_f_trio <- glm(adhd ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
clustered_se <- vcovCL(adhd_f_trio, cluster = ~ m_id)
summary_cl <- coeftest(adhd_f_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd_f_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_f_trio <- adhd_f_trio %>% subset(term =="diet_pgs_far")
adhd_f_trio <- adhd_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD")

##Child ADD

add_c <- glm(adhd_att ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
clustered_se <- vcovCL(add_c, cluster = ~ m_id)
summary_cl <- coeftest(add_c, vcov. = clustered_se) %>% tidy(conf.int=T)
add_c <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_c <- add_c %>% subset(term =="diet_pgs_barn")
add_c <- add_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_c_trio <- glm(adhd_att ~ diet_pgs_barn + diet_pgs_far + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(add_c_trio, cluster = ~ m_id)
summary_cl <- coeftest(add_c_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
add_c_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_c_trio <- add_c_trio %>% subset(term =="diet_pgs_barn")
add_c_trio <- add_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Mother ADD

add_m <- glm(adhd_att ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(add_m, cluster = ~ m_id)
summary_cl <- coeftest(add_m, vcov. = clustered_se) %>% tidy(conf.int=T)
add_m <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_m <- add_m %>% subset(term =="diet_pgs_mor")
add_m <- add_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_m_trio <- glm(adhd_att ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
clustered_se <- vcovCL(add_m_trio, cluster = ~ m_id)
summary_cl <- coeftest(add_m_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
add_m_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_m_trio <- add_m_trio %>% subset(term =="diet_pgs_mor")
add_m_trio <- add_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Father ADD

add_f <- glm(adhd_att ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
clustered_se <- vcovCL(add_f, cluster = ~ m_id)
summary_cl <- coeftest(add_f, vcov. = clustered_se) %>% tidy(conf.int=T)
add_f <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_f <- add_f %>% subset(term =="diet_pgs_far")
add_f <- add_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_f_trio <- glm(adhd_att ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(add_f_trio, cluster = ~ m_id)
summary_cl <- coeftest(add_f_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
add_f_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_f_trio <- add_f_trio %>% subset(term =="diet_pgs_far")
add_f_trio <- add_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Child ADHD - Combined Presentation

adhd1_c <- glm(adhd_comb ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd1_c, cluster = ~ m_id)
summary_cl <- coeftest(adhd1_c, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd1_c <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_c <- adhd1_c %>% subset(term =="diet_pgs_barn")
adhd1_c <- adhd1_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_c_trio <- glm(adhd_comb ~ diet_pgs_barn + diet_pgs_far + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd1_c_trio, cluster = ~ m_id)
summary_cl <- coeftest(adhd1_c_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd1_c_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_c_trio <- adhd1_c_trio %>% subset(term =="diet_pgs_barn")
adhd1_c_trio <- adhd1_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

##Mother ADHD - Combined Presentation

adhd1_m <- glm(adhd_comb ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd1_m, cluster = ~ m_id)
summary_cl <- coeftest(adhd1_m, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd1_m <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_m <- adhd1_m %>% subset(term =="diet_pgs_mor")
adhd1_m <- adhd1_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_m_trio <- glm(adhd_comb ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd1_m_trio, cluster = ~ m_id)
summary_cl <- coeftest(adhd1_m_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd1_m_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_m_trio <- adhd1_m_trio %>% subset(term =="diet_pgs_mor")
adhd1_m_trio <- adhd1_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

##Father ADHD - Combined Presentation

adhd1_f <- glm(adhd_comb ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd1_f, cluster = ~ m_id)
summary_cl <- coeftest(adhd1_f, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd1_f <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_f <- adhd1_f %>% subset(term =="diet_pgs_far")
adhd1_f <- adhd1_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_f_trio <- glm(adhd_comb ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd1_f_trio, cluster = ~ m_id)
summary_cl <- coeftest(adhd1_f_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd1_f_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_f_trio <- adhd1_f_trio %>% subset(term =="diet_pgs_far")
adhd1_f_trio <- adhd1_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

forestadhd <- rbind(adhd_c, adhd_c_trio, adhd_m, adhd_m_trio, adhd_f, adhd_f_trio, add_c, add_c_trio, add_m, add_m_trio, add_f, add_f_trio, adhd1_c, adhd1_c_trio, adhd1_m, adhd1_m_trio, adhd1_f, adhd1_f_trio)
forestadhd$Condition <- factor(forestadhd$Condition, levels=c("norm", "trio"))
forestadhd$Group <- factor(forestadhd$Group, levels=c("Child direct effect", "Child PGS", "Paternal indirect effect", "Paternal PGS", "Maternal indirect effect", "Maternal PGS"))
forestadhd$Copsych <- factor(forestadhd$Copsych, levels=c("ADHD", "ADHD  Inattentive Presentation", "ADHD  Combined Presentation"))
forestadhd <- forestadhd %>% mutate(significant = ifelse(p.value<0.05,1,0))

theme_set(theme_bw())
redPalette <- c("#E69F00", "#009E73")

p_adhd_trio = ggplot(data=forestadhd, aes(x = Group,y = estimate, ymin = conf.low, ymax = conf.high, na.rm = TRUE))+
  geom_point(size = 2.0, shape = 19, aes(col=Condition)) +
  geom_hline(aes(fill=Group),yintercept =1, linetype=2)+
  xlab("") + ylab("Odds Ratio (95% Confidence Interval)") +
  coord_flip() + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high,col=Condition),width=0.2,cex=0.5)+ 
  scale_y_continuous(limits = c(0.75,1.3))+ 
  facet_wrap(~Copsych,strip.position="top",nrow=1,scales = "free_x", labeller = label_wrap_gen(width=8) ) +
  scale_colour_manual(values=redPalette) +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_text(size = 14, face="bold"),
        axis.text.x=element_text(size = 16, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 16, face="bold", margin = margin(0.6,0,0.6,0, "cm")), legend.position="None", strip.background = element_rect(fill = "white")) +
  geom_text(aes(label = ifelse(significant == 1, "*", ""), col=Condition), vjust = -0.1, size = 6) 


## ADHD RS-DBD Total traits, Inattention and IH traits

##Total traits
##Child ADHD

adhd_total_symp_c <- glm(rsdbd_adhd_c_8yr ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), data = merged_data)
clustered_se <- vcovCL(adhd_total_symp_c, cluster = ~ m_id)
adhd_total_symp_c <- coeftest(adhd_total_symp_c, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_c <- adhd_total_symp_c %>% subset(term =="diet_pgs_barn")
adhd_total_symp_c <- adhd_total_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Total traits")

adhd_total_symp_c_trio <- glm(rsdbd_adhd_c_8yr ~ diet_pgs_barn + diet_pgs_mor + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data)
clustered_se <- vcovCL(adhd_total_symp_c_trio, cluster = ~ m_id)
adhd_total_symp_c_trio <- coeftest(adhd_total_symp_c_trio, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_c_trio <- adhd_total_symp_c_trio %>% subset(term =="diet_pgs_barn")
adhd_total_symp_c_trio <- adhd_total_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Total traits")

##Mother ADHD

adhd_total_symp_m <- glm(rsdbd_adhd_c_8yr ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
clustered_se <- vcovCL(adhd_total_symp_m, cluster = ~ m_id)
adhd_total_symp_m <- coeftest(adhd_total_symp_m, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_m <- adhd_total_symp_m %>% subset(term =="diet_pgs_mor")
adhd_total_symp_m <- adhd_total_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Total traits")

adhd_total_symp_m_trio <- glm(rsdbd_adhd_c_8yr ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data)
clustered_se <- vcovCL(adhd_total_symp_m_trio, cluster = ~ m_id)
adhd_total_symp_m_trio <- coeftest(adhd_total_symp_m_trio, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_m_trio <- adhd_total_symp_m_trio %>% subset(term =="diet_pgs_mor")
adhd_total_symp_m_trio <- adhd_total_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Total traits")

##Father ADHD

adhd_total_symp_f <- glm(rsdbd_adhd_c_8yr ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
clustered_se <- vcovCL(adhd_total_symp_f, cluster = ~ m_id)
adhd_total_symp_f <- coeftest(adhd_total_symp_f, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_f <- adhd_total_symp_f %>% subset(term =="diet_pgs_far")
adhd_total_symp_f <- adhd_total_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Total traits")

adhd_total_symp_f_trio <- glm(rsdbd_adhd_c_8yr ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
clustered_se <- vcovCL(adhd_total_symp_f_trio, cluster = ~ m_id)
adhd_total_symp_f_trio <- coeftest(adhd_total_symp_f_trio, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_total_symp_f_trio <- adhd_total_symp_f_trio %>% subset(term =="diet_pgs_far")
adhd_total_symp_f_trio <- adhd_total_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Total traits")

##Inattention
##Child ADHD

adhd_att_symp_c <- glm(rsdbd_ina_c_8yr ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
clustered_se <- vcovCL(adhd_att_symp_c, cluster = ~ m_id)
adhd_att_symp_c <- coeftest(adhd_att_symp_c, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_c <- adhd_att_symp_c %>% subset(term =="diet_pgs_barn")
adhd_att_symp_c <- adhd_att_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Inattention traits")

adhd_att_symp_c_trio <- glm(rsdbd_ina_c_8yr ~ diet_pgs_barn + diet_pgs_mor + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
clustered_se <- vcovCL(adhd_att_symp_c_trio, cluster = ~ m_id)
adhd_att_symp_c_trio <- coeftest(adhd_att_symp_c_trio, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_c_trio <- adhd_att_symp_c_trio %>% subset(term =="diet_pgs_barn")
adhd_att_symp_c_trio <- adhd_att_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Inattention traits")

##Mother ADHD

adhd_att_symp_m <- glm(rsdbd_ina_c_8yr ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
clustered_se <- vcovCL(adhd_att_symp_m, cluster = ~ m_id)
adhd_att_symp_m <- coeftest(adhd_att_symp_m, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_m <- adhd_att_symp_m %>% subset(term =="diet_pgs_mor")
adhd_att_symp_m <- adhd_att_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Inattention traits")

adhd_att_symp_m_trio <- glm(rsdbd_ina_c_8yr ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
clustered_se <- vcovCL(adhd_att_symp_m_trio, cluster = ~ m_id)
adhd_att_symp_m_trio <- coeftest(adhd_att_symp_m_trio, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_m_trio <- adhd_att_symp_m_trio %>% subset(term =="diet_pgs_mor")
adhd_att_symp_m_trio <- adhd_att_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Inattention traits")

##Father ADHD

adhd_att_symp_f <- glm(rsdbd_ina_c_8yr ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
clustered_se <- vcovCL(adhd_att_symp_f, cluster = ~ m_id)
adhd_att_symp_f <- coeftest(adhd_att_symp_f, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_f <- adhd_att_symp_f %>% subset(term =="diet_pgs_far")
adhd_att_symp_f <- adhd_att_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Inattention traits")

adhd_att_symp_f_trio <- glm(rsdbd_ina_c_8yr ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data)
clustered_se <- vcovCL(adhd_att_symp_f_trio, cluster = ~ m_id)
adhd_att_symp_f_trio <- coeftest(adhd_att_symp_f_trio, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_att_symp_f_trio <- adhd_att_symp_f_trio %>% subset(term =="diet_pgs_far")
adhd_att_symp_f_trio <- adhd_att_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Inattention traits")

##Hyperactivity/Impulsivity
##Child ADHD

adhd_ih_symp_c <- glm(rsdbd_hyp_c_8yr ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
clustered_se <- vcovCL(adhd_ih_symp_c, cluster = ~ m_id)
adhd_ih_symp_c <- coeftest(adhd_ih_symp_c, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_c <- adhd_ih_symp_c %>% subset(term =="diet_pgs_barn")
adhd_ih_symp_c <- adhd_ih_symp_c %>% mutate(Group = "Child PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Hyperactivity/Impulsivity traits")

adhd_ih_symp_c_trio <- glm(rsdbd_hyp_c_8yr ~ diet_pgs_barn + diet_pgs_mor + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data)
clustered_se <- vcovCL(adhd_ih_symp_c_trio, cluster = ~ m_id)
adhd_ih_symp_c_trio <- coeftest(adhd_ih_symp_c_trio, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_c_trio <- adhd_ih_symp_c_trio %>% subset(term =="diet_pgs_barn")
adhd_ih_symp_c_trio <- adhd_ih_symp_c_trio %>% mutate(Group = "Child direct effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Hyperactivity/Impulsivity traits")

##Mother ADHD

adhd_ih_symp_m <- glm(rsdbd_hyp_c_8yr ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
clustered_se <- vcovCL(adhd_ih_symp_m, cluster = ~ m_id)
adhd_ih_symp_m <- coeftest(adhd_ih_symp_m, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_m <- adhd_ih_symp_m %>% subset(term =="diet_pgs_mor")
adhd_ih_symp_m <- adhd_ih_symp_m %>% mutate(Group = "Maternal PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Hyperactivity/Impulsivity traits")

adhd_ih_symp_m_trio <- glm(rsdbd_hyp_c_8yr ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data)
clustered_se <- vcovCL(adhd_ih_symp_m_trio, cluster = ~ m_id)
adhd_ih_symp_m_trio <- coeftest(adhd_ih_symp_m_trio, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_m_trio <- adhd_ih_symp_m_trio %>% subset(term =="diet_pgs_mor")
adhd_ih_symp_m_trio <- adhd_ih_symp_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Hyperactivity/Impulsivity traits")

##Father ADHD

adhd_ih_symp_f <- glm(rsdbd_hyp_c_8yr ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
clustered_se <- vcovCL(adhd_ih_symp_f, cluster = ~ m_id)
adhd_ih_symp_f <- coeftest(adhd_ih_symp_f, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_f <- adhd_ih_symp_f %>% subset(term =="diet_pgs_far")
adhd_ih_symp_f <- adhd_ih_symp_f %>% mutate(Group = "Paternal PGS", Condition = "Adjusted for sex and birth year", Copsych = "ADHD RS-DBD Hyperactivity/Impulsivity traits")

adhd_ih_symp_f_trio <- glm(rsdbd_hyp_c_8yr ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), data = merged_data) 
clustered_se <- vcovCL(adhd_ih_symp_f_trio, cluster = ~ m_id)
adhd_ih_symp_f_trio <- coeftest(adhd_ih_symp_f_trio, vcov. = clustered_se) %>% broom::tidy(conf.int=T) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_ih_symp_f_trio <- adhd_ih_symp_f_trio %>% subset(term =="diet_pgs_far")
adhd_ih_symp_f_trio <- adhd_ih_symp_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "Trio model", Copsych = "ADHD RS-DBD Hyperactivity/Impulsivity traits")


##Forest plot ADHD RS traits 

forestadhdrssub <- rbind(adhd_total_symp_c, adhd_total_symp_c_trio, adhd_total_symp_m, adhd_total_symp_m_trio, adhd_total_symp_f, adhd_total_symp_f_trio, adhd_att_symp_c, adhd_att_symp_c_trio, adhd_att_symp_m, adhd_att_symp_m_trio, adhd_att_symp_f, adhd_att_symp_f_trio, adhd_ih_symp_c, adhd_ih_symp_c_trio, adhd_ih_symp_m, adhd_ih_symp_m_trio, adhd_ih_symp_f, adhd_ih_symp_f_trio)

forestadhdrssub$Condition <- factor(forestadhdrssub$Condition, levels=c("Adjusted for sex and birth year", "Trio model"))
forestadhdrssub$Group <- factor(forestadhdrssub$Group, levels=c("Child direct effect", "Child PGS", "Paternal indirect effect", "Paternal PGS", "Maternal indirect effect", "Maternal PGS"))
forestadhdrssub$Copsych <- factor(forestadhdrssub$Copsych, levels=c("ADHD RS-DBD Total traits","ADHD RS-DBD Inattention traits","ADHD RS-DBD Hyperactivity/Impulsivity traits"))
forestadhdrssub <- forestadhdrssub %>% mutate(significant = ifelse(p.value<0.05,1,0))

theme_set(theme_bw())
redPalette <- c("#E69F00", "#009E73")

p_adhd_rs_sub_trio = ggplot(data=forestadhdrssub, aes(x = Group,y = estimate, ymin = conf.low, ymax = conf.high, na.rm = TRUE))+
  geom_point(size = 2.0, shape = 19, aes(col=Condition)) +
  geom_hline(aes(fill=Group),yintercept =0, linetype=2)+
  xlab("") + ylab("Beta Estimate (95% Confidence Interval)") +
  coord_flip() + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high,col=Condition),width=0.2,cex=0.5)+ 
  scale_y_continuous(limits = c(-0.5,0.5))+  # Set x-axis limits
  facet_wrap(~Copsych,strip.position="top",nrow=1,scales = "free_x", labeller = label_wrap_gen (width=15) ) +
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









################################ Results for ADHD with a minimum of two registered diagnoses ##############################



################################## Figures - forest plot trio model ###########################

##Diagnoses

##Child ADHD

adhd_c <- glm(adhd2 ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
clustered_se <- vcovCL(adhd_c, cluster = ~ m_id)
summary_cl <- coeftest(adhd_c, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd_c <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_c <- adhd_c %>% subset(term =="diet_pgs_barn")
adhd_c <- adhd_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD")

adhd_c_trio <- glm(adhd2 ~ diet_pgs_barn + diet_pgs_far + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd_c_trio, cluster = ~ m_id)
summary_cl <- coeftest(adhd_c_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd_c_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_c_trio <- adhd_c_trio %>% subset(term =="diet_pgs_barn")
adhd_c_trio <- adhd_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD")

##Mother ADHD

adhd_m <- glm(adhd2 ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd_m, cluster = ~ m_id)
summary_cl <- coeftest(adhd_m, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd_m <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_m <- adhd_m %>% subset(term =="diet_pgs_mor")
adhd_m <- adhd_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD")

adhd_m_trio <- glm(adhd2 ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd_m_trio, cluster = ~ m_id)
summary_cl <- coeftest(adhd_m_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd_m_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_m_trio <- adhd_m_trio %>% subset(term =="diet_pgs_mor")
adhd_m_trio <- adhd_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD")

##Father ADHD

adhd_f <- glm(adhd2 ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
clustered_se <- vcovCL(adhd_f, cluster = ~ m_id)
summary_cl <- coeftest(adhd_f, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd_f <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_f <- adhd_f %>% subset(term =="diet_pgs_far")
adhd_f <- adhd_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD")

adhd_f_trio <- glm(adhd2 ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
clustered_se <- vcovCL(adhd_f_trio, cluster = ~ m_id)
summary_cl <- coeftest(adhd_f_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd_f_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd_f_trio <- adhd_f_trio %>% subset(term =="diet_pgs_far")
adhd_f_trio <- adhd_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD")

##Child ADD

add_c <- glm(adhd_att2 ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
clustered_se <- vcovCL(add_c, cluster = ~ m_id)
summary_cl <- coeftest(add_c, vcov. = clustered_se) %>% tidy(conf.int=T)
add_c <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_c <- add_c %>% subset(term =="diet_pgs_barn")
add_c <- add_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_c_trio <- glm(adhd_att2 ~ diet_pgs_barn + diet_pgs_far + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(add_c_trio, cluster = ~ m_id)
summary_cl <- coeftest(add_c_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
add_c_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_c_trio <- add_c_trio %>% subset(term =="diet_pgs_barn")
add_c_trio <- add_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Mother ADD

add_m <- glm(adhd_att2 ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(add_m, cluster = ~ m_id)
summary_cl <- coeftest(add_m, vcov. = clustered_se) %>% tidy(conf.int=T)
add_m <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_m <- add_m %>% subset(term =="diet_pgs_mor")
add_m <- add_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_m_trio <- glm(adhd_att2 ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
clustered_se <- vcovCL(add_m_trio, cluster = ~ m_id)
summary_cl <- coeftest(add_m_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
add_m_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_m_trio <- add_m_trio %>% subset(term =="diet_pgs_mor")
add_m_trio <- add_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Father ADD

add_f <- glm(adhd_att2 ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data) 
clustered_se <- vcovCL(add_f, cluster = ~ m_id)
summary_cl <- coeftest(add_f, vcov. = clustered_se) %>% tidy(conf.int=T)
add_f <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_f <- add_f %>% subset(term =="diet_pgs_far")
add_f <- add_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD  Inattentive Presentation")

add_f_trio <- glm(adhd_att2 ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(add_f_trio, cluster = ~ m_id)
summary_cl <- coeftest(add_f_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
add_f_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
add_f_trio <- add_f_trio %>% subset(term =="diet_pgs_far")
add_f_trio <- add_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD  Inattentive Presentation")

##Child ADHD - Combined Presentation

adhd1_c <- glm(adhd_comb2 ~ diet_pgs_barn + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd1_c, cluster = ~ m_id)
summary_cl <- coeftest(adhd1_c, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd1_c <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_c <- adhd1_c %>% subset(term =="diet_pgs_barn")
adhd1_c <- adhd1_c %>% mutate(Group = "Child PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_c_trio <- glm(adhd_comb2 ~ diet_pgs_barn + diet_pgs_far + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd1_c_trio, cluster = ~ m_id)
summary_cl <- coeftest(adhd1_c_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd1_c_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_c_trio <- adhd1_c_trio %>% subset(term =="diet_pgs_barn")
adhd1_c_trio <- adhd1_c_trio %>% mutate(Group = "Child direct effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

##Mother ADHD - Combined Presentation

adhd1_m <- glm(adhd_comb2 ~ diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd1_m, cluster = ~ m_id)
summary_cl <- coeftest(adhd1_m, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd1_m <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_m <- adhd1_m %>% subset(term =="diet_pgs_mor")
adhd1_m <- adhd1_m %>% mutate(Group = "Maternal PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_m_trio <- glm(adhd_comb2 ~ diet_pgs_mor + diet_pgs_barn + diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd1_m_trio, cluster = ~ m_id)
summary_cl <- coeftest(adhd1_m_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd1_m_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_m_trio <- adhd1_m_trio %>% subset(term =="diet_pgs_mor")
adhd1_m_trio <- adhd1_m_trio %>% mutate(Group = "Maternal indirect effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

##Father ADHD - Combined Presentation

adhd1_f <- glm(adhd_comb2 ~ diet_pgs_far + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd1_f, cluster = ~ m_id)
summary_cl <- coeftest(adhd1_f, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd1_f <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_f <- adhd1_f %>% subset(term =="diet_pgs_far")
adhd1_f <- adhd1_f %>% mutate(Group = "Paternal PGS", Condition = "norm", Copsych = "ADHD  Combined Presentation")

adhd1_f_trio <- glm(adhd_comb2 ~ diet_pgs_far + diet_pgs_barn + diet_pgs_mor + as.numeric(birth_yr) + as.factor(sex), family = "binomial", data = merged_data)
clustered_se <- vcovCL(adhd1_f_trio, cluster = ~ m_id)
summary_cl <- coeftest(adhd1_f_trio, vcov. = clustered_se) %>% tidy(conf.int=T)
adhd1_f_trio <- summary_cl %>% mutate(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% select(estimate, p.value, conf.low, conf.high, term) 
adhd1_f_trio <- adhd1_f_trio %>% subset(term =="diet_pgs_far")
adhd1_f_trio <- adhd1_f_trio %>% mutate(Group = "Paternal indirect effect", Condition = "trio", Copsych = "ADHD  Combined Presentation")

forestadhd <- rbind(adhd_c, adhd_c_trio, adhd_m, adhd_m_trio, adhd_f, adhd_f_trio, add_c, add_c_trio, add_m, add_m_trio, add_f, add_f_trio, adhd1_c, adhd1_c_trio, adhd1_m, adhd1_m_trio, adhd1_f, adhd1_f_trio)
forestadhd$Condition <- factor(forestadhd$Condition, levels=c("norm", "trio"))
forestadhd$Group <- factor(forestadhd$Group, levels=c("Child direct effect", "Child PGS", "Paternal indirect effect", "Paternal PGS", "Maternal indirect effect", "Maternal PGS"))
forestadhd$Copsych <- factor(forestadhd$Copsych, levels=c("ADHD", "ADHD  Inattentive Presentation", "ADHD  Combined Presentation"))
forestadhd <- forestadhd %>% mutate(significant = ifelse(p.value<0.05,1,0))

theme_set(theme_bw())
redPalette <- c("#E69F00", "#009E73")

p_adhd_trio_2diag = ggplot(data=forestadhd, aes(x = Group,y = estimate, ymin = conf.low, ymax = conf.high, na.rm = TRUE))+
  geom_point(size = 2.0, shape = 19, aes(col=Condition)) +
  geom_hline(aes(fill=Group),yintercept =1, linetype=2)+
  xlab("") + ylab("Odds Ratio (95% Confidence Interval)") +
  coord_flip() + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high,col=Condition),width=0.2,cex=0.5)+ 
  scale_y_continuous(limits = c(0.75,1.3))+ 
  facet_wrap(~Copsych,strip.position="top",nrow=1,scales = "free_x", labeller = label_wrap_gen(width=8) ) +
  scale_colour_manual(values=redPalette) +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_text(size = 14, face="bold"),
        axis.text.x=element_text(size = 16, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 16, face="bold", margin = margin(0.6,0,0.6,0, "cm")), legend.position="None", strip.background = element_rect(fill = "white")) +
  geom_text(aes(label = ifelse(significant == 1, "*", ""), col=Condition), vjust = -0.1, size = 6) 

