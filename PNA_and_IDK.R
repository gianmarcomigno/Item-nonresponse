# import full phenotypes dataset
d <- fread("ukb31063.phenotypes.20191008.csv", header=T)
colnames(d) <- paste0("f.",gsub("x","",gsub("_",".",colnames(d))))
d1 <- d %>% select(grep("20400", colnames(d)), grep("f.eid", colnames(d)),  # 110003 has only NAs, 110005 not found
                   grep("12188", colnames(d)),                             
                   grep("110001", colnames(d)), grep("110005", colnames(d)),
                   grep("2000", colnames(d))) 
saveRDS(d1, file = "recontact_variables.rds")

# read in all files
alldata <- readRDS("alldata.rds")
covariates <- fread("ukb31063.neale_gwas_covariates.both_sexes.csv", skip = 1)
covariates <- covariates[covariates$s %in% alldata$f.userId,] # subset
covariates <- covariates[match(alldata$f.userId, covariates$s),] # sort for userID

pna_eba <- data.frame(readRDS(file = "factscores_PNA.rds"))['g']
idk_eba <- data.frame(readRDS(file = "factscores_IDK.rds"))['g']

# COVARIATES: dataframe with all the  covariates in the model
covariate_fin <- as.data.frame(
  cbind(covariates, # female, age, age2, age*female, age2*female, PC1--20
        alldata$f.54.0.0, alldata$f.2178.0.0, alldata$f.6138.0.0, # f.54: assessment center, f.6138: education, f.2138: health
        scale(idk_eba), # scaled idk factorscores phenotypes
        scale(pna_eba)) # scaled pna factorscores phenotypes
)
colnames(covariate_fin)[c(1,27:31)] <- c("userID","assessm_center","health","eduYears", "scaled_idk", "scaled_pna")

# recode
covariate_fin$health <- car::recode(covariate_fin$health,   # overall health (1218 idk, 113 pna)
                                    "-1=NA ; -3=NA ")
covariate_fin$isFemale <- ifelse(covariate_fin$isFemale == T, 1,0) #1 is female
covariate_fin$eduYears <- car::recode(covariate_fin$eduYears,  # impute f.6138 as eduyears in paper ea gwas (james lee et al.)
                                      "1=20 ; 2=15 ; 3=13 ; 4=12 ;
                                      5=19 ; 6=17 ; -7=6 ; -3=NA ")

#--------------------------------------------------------
#--------------------------------------------------------
recontact_variables <- readRDS("recontact_variables.rds")
covariate_fin <- readRDS("covariate_fin.rds")

# Field 110001-1.0: Invitation to complete online 24-hour recall dietary questionnaire - acceptance (1st istance)
acceptance <- recontact_variables %>% select("f.eid", "f.110001.1.0") %>% 
  mutate(f.110001.1.0 = if_else(f.110001.1.0 == 0, 1, 0)) %>% #recode (0 non reponse, 2 completed)
  filter(f.eid %in% covariate_fin$userID) %>% #subset
  left_join(covariate_fin, by = c("f.eid" = "userID")) %>%
  rename(nonresponded_diet_questionn_1visit = f.110001.1.0) %>% 
  na.omit 

#logistic model
model1 <- glm(nonresponded_diet_questionn_1visit ~ isFemale + age + age_squared + age_isFemale + age_squared_isFemale + PC1 + PC2 + PC3 + PC4 +
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +  
                health +
                eduYears +
                scaled_idk +   #0.0003
                scaled_pna,     #0.0010
              data = acceptance, family = "binomial")

summary(model1)

blr_rsq_mckelvey_zavoina(model1)

# Field 110001-1.0: Invitation to complete online 24-hour recall dietary questionnaire - acceptance (1st istance) - 
# NON-RESPONDED ALL VISITS (1) vs COMPLETED ALL VISITS (0)
none_vs_all <- recontact_variables %>% select("f.eid", "f.110001.1.0", "f.110001.2.0", "f.110001.3.0", "f.110001.4.0") %>% 
  mutate(none = case_when(
    f.110001.1.0+f.110001.2.0+f.110001.3.0+f.110001.4.0 == 8 ~ 0, 
    f.110001.1.0+f.110001.2.0+f.110001.3.0+f.110001.4.0 == 0 ~ 1, 
    TRUE ~ NA_real_)) %>% #recode (1 non-reponded to all visits, 0 completed all visits, NA: completed some visits)
  filter(f.eid %in% covariate_fin$userID) %>% #subset
  left_join(covariate_fin, by = c("f.eid" = "userID")) %>%
  rename(nonresponded_all_visits = none) %>% 
  select(-c("f.110001.1.0", "f.110001.2.0", "f.110001.3.0", "f.110001.4.0")) %>%
  na.omit()

 model0 <- glm(nonresponded_all_visits ~ isFemale + age + age_squared + age_isFemale + age_squared_isFemale + PC1 + PC2 + PC3 + PC4 +
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +  
                health +
                eduYears +
               scaled_idk +   #0.0012
               scaled_pna,     #0.0027
              data = none_vs_all, family = "binomial")
 
summary(model0)

 
blr_rsq_mckelvey_zavoina(model1) #https://www.slideshare.net/MarcoDAlessandro11/pseudor-quadro for different methods

###########################
# ROC curve logistic regression
#####################"######
install.packages("pROC")
library(pROC)

mydata <- acceptance
model1 <- glm(nonresponded_diet_questionn_1visit ~ isFemale + age + age_squared + age_isFemale + age_squared_isFemale + PC1 + PC2 + PC3 + PC4 +
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + 
                eduYears + assessm_center + health + scaled_pna + scaled_idk, 
              data = acceptance, family = "binomial")
acceptance$prob <- predict(model1,type=c("response"))
plot(roc(nonresponded_diet_questionn_1visit ~ prob, data = acceptance))
